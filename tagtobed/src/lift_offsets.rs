use flate2::read::GzDecoder;
use rayon::prelude::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Seek, SeekFrom};
use std::sync::mpsc::sync_channel;

mod gfa;

fn main() {
    let cpgs_file = std::env::args().nth(1).expect("Missing cpgs file");
    let gfa_file = std::env::args().nth(2).expect("Missing GFA file");

    // Build nuc_index from gzipped cpgs file
    let mut nuc_index: HashMap<usize, Vec<(usize, String, String)>> = HashMap::new();
    let gz = GzDecoder::new(File::open(cpgs_file).unwrap());
    let reader = BufReader::new(gz);

    for line in reader.lines() {
        let line = line.unwrap();
        let fields: Vec<&str> = line.trim().split('\t').collect();
        if fields.len() != 4 {
            continue;
        }
        let node = fields[0].parse::<usize>().unwrap();
        let offset = fields[1].parse::<usize>().unwrap();
        let strand = fields[2];
        let nucset = fields[3];
        if nucset.contains(' ') {
            continue;
        }
        if strand == "-" {
            continue;
        }
        nuc_index.entry(node).or_insert_with(Vec::new).push((
            offset,
            strand.to_string(),
            nucset.to_string(),
        ));
    }

    // Build node_lengths from GFA
    let mut gfa = File::open(&gfa_file).unwrap();
    let reader = BufReader::new(&gfa);
    let mut node_lengths: HashMap<usize, usize> = HashMap::new();

    for line in reader.lines() {
        let line = line.unwrap();
        if line.starts_with('S') {
            let fields: Vec<&str> = line.trim().split_whitespace().collect();
            if fields.len() < 3 {
                continue;
            }
            let node_id = fields[1].parse::<usize>().unwrap();
            let node_seq = fields[2];
            node_lengths.insert(node_id, node_seq.len());
        }
    }

    // Rewind GFA file and process paths and walks in parallel
    gfa.seek(SeekFrom::Start(0)).unwrap();
    let reader = BufReader::new(gfa);

    let nuc_index = &nuc_index;
    let node_lengths = &node_lengths;

    let (tx, rx) = sync_channel::<String>(100000);

    let writer = std::thread::spawn(move || {
        for line in rx {
            println!("{}", line);
        }
    });

    reader
        .lines()
        .par_bridge()
        .map(Result::unwrap)
        .for_each(|line| {
            let Some(path) = gfa::parse_path(&line).expect("Invalid GFA path or walk") else {
                return;
            };
            let tx_ = tx.clone();
            let mut i = 0usize;
            for (node_name, strand) in path.nodes {
                if let Some(cpgs) = nuc_index.get(&node_name) {
                    for cpg in cpgs {
                        let mut offset = cpg.0;
                        if strand == "-" {
                            offset = node_lengths.get(&node_name).unwrap() - offset - 1;
                        }
                        tx_.send(format!(
                            "{}\t{}\t{}\t{}",
                            path.name,
                            path.start + i + offset,
                            path.start + i + offset + 2,
                            cpg.2
                        ))
                        .unwrap();
                    }
                }
                i += node_lengths.get(&node_name).unwrap();
            }
        });
    drop(tx);
    writer.join().unwrap();
}
