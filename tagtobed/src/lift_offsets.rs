use std::fs::File;
use std::io::{BufRead, BufReader, Seek, SeekFrom};
use std::collections::HashMap;
use flate2::read::GzDecoder;
use std::sync::mpsc::sync_channel;
use rayon::prelude::*;

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
        if fields.len() != 4 { continue; }
        let node = fields[0].parse::<usize>().unwrap();
        let offset = fields[1].parse::<usize>().unwrap();
        let strand = fields[2];
        let nucset = fields[3];
        if nucset.contains(' ') { continue; }
        if strand == "-" { continue; }
        nuc_index.entry(node)
            .or_insert_with(Vec::new)
            .push((offset, strand.to_string(), nucset.to_string()));
    }

    // Build node_lengths from GFA
    let mut gfa = File::open(&gfa_file).unwrap();
    let reader = BufReader::new(&gfa);
    let mut node_lengths: HashMap<usize, usize> = HashMap::new();

    for line in reader.lines() {
        let line = line.unwrap();
        if line.starts_with('S') {
            let fields: Vec<&str> = line.trim().split_whitespace().collect();
            if fields.len() < 3 { continue; }
            let node_id = fields[1].parse::<usize>().unwrap();
            let node_seq = fields[2];
            node_lengths.insert(node_id, node_seq.len());
        }
    }

    // Rewind GFA file and process P lines in parallel
    gfa.seek(SeekFrom::Start(0)).unwrap();
    let reader = BufReader::new(gfa);

    let nuc_index = &nuc_index;
    let node_lengths = &node_lengths;

    let (tx, rx) = sync_channel::<String>(10000);

    std::thread::spawn(move || {
        for line in rx {
            println!("{}", line);
        }
    });

    reader.lines()
        .par_bridge()
        .filter_map(Result::ok)
        .filter(|l| l.starts_with('P'))
        .for_each(|line| {
            let tx_ = tx.clone();
            let fields: Vec<&str> = line.trim().split_whitespace().collect();
            if fields.len() < 4 { return; }
            let p_name = fields[1];
            let hap = fields[2];
            let mut hap_name = p_name.to_string();
            let mut hap_start = 0usize;

            let hap_seq: Vec<&str> = hap.split(',').collect();

            if p_name.contains('[') {
                let parts: Vec<&str> = p_name.split('[').collect();
                hap_name = parts[0].to_string();
                let start = parts[1].split('-').next().unwrap();
                hap_start = start.parse::<usize>().unwrap();
            }

            let mut i = 0usize;
            for node in hap_seq {
                let node_name = node[..node.len()-1].parse::<usize>().unwrap();
                let strand = &node[node.len()-1..];

                if let Some(cpgs) = nuc_index.get(&node_name) {
                    for cpg in cpgs {
                        let mut offset = cpg.0;
                        if strand == "-" {
                            offset = node_lengths.get(&node_name).unwrap() - offset - 1;
                        }
                        tx_.send(format!(
                            "{}\t{}\t{}\t{}",
                            hap_name,
                            hap_start + i + offset,
                            hap_start + i + offset + 2,
                            cpg.2
                        ));
                    }
                }
                i += node_lengths.get(&node_name).unwrap();
            }
        });
    drop(tx);
}
