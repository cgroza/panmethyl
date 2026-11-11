use flate2::read::GzDecoder;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, Seek, SeekFrom};
use std::sync::mpsc::sync_channel;

mod gfa;

fn flip(s: &str) -> &str {
    match s {
        "+" => "-",
        "-" => "+",
        _ => s,
    }
}

fn complement(edge: &(String, usize, String, usize)) -> (String, usize, String, usize) {
    (
        flip(&edge.2).to_string(),
        edge.3,
        flip(&edge.0).to_string(),
        edge.1,
    )
}

fn main() {
    let strands = HashMap::from([(">", "+"), ("<", "-")]);
    let strands_ = HashMap::from([("+", ">"), ("-", "<")]);

    // Read edge index from gzipped file
    let edge_index_file = std::env::args().nth(1).expect("Missing edge index file");
    let gfa_file = std::env::args().nth(2).expect("Missing GFA file");

    let mut edge_index = HashSet::new();
    let gz = GzDecoder::new(File::open(edge_index_file).unwrap());
    let reader = BufReader::new(gz);

    for line in reader.lines() {
        let line = line.unwrap();
        let fields: Vec<&str> = line.trim().split('\t').collect();
        if fields.len() != 4 {
            continue;
        }
        let nucset = fields[3];
        if !nucset.contains(' ') {
            continue;
        }
        let parts: Vec<&str> = nucset.split(' ').collect();
        let left = parts[0];
        let right = parts[1];
        let lstrand = strands[left.chars().next().unwrap().to_string().as_str()].to_string();
        let rstrand = strands[right.chars().next().unwrap().to_string().as_str()].to_string();
        let lnode = left[1..].parse::<usize>().unwrap();
        let rnode = right[1..].parse::<usize>().unwrap();
        edge_index.insert((lstrand, lnode, rstrand, rnode));
    }

    // Build node_lengths from GFA
    let mut gfa = File::open(&gfa_file).unwrap();
    let reader = BufReader::new(&gfa);
    let mut node_lengths = HashMap::new();

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

    let edge_index = &edge_index;
    let node_lengths = &node_lengths;
    let strands_ = &strands_;

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
            for nodes in path.nodes.windows(2) {
                let edge = (
                    nodes[0].1.to_string(),
                    nodes[0].0,
                    nodes[1].1.to_string(),
                    nodes[1].0,
                );
                let cedge = complement(&edge);
                let mut edge_name = format!(
                    "{}{} {}{}",
                    strands_[&edge.0.as_str()],
                    edge.1,
                    strands_[&edge.2.as_str()],
                    edge.3
                );
                if edge_index.contains(&cedge) {
                    edge_name = format!(
                        "{}{} {}{}",
                        strands_[&cedge.0.as_str()],
                        cedge.1,
                        strands_[&cedge.2.as_str()],
                        cedge.3
                    );
                }
                if edge_index.contains(&edge) || edge_index.contains(&cedge) {
                    let node_len = *node_lengths.get(&edge.1).unwrap();
                    let start = path.start + i + node_len - 1;
                    let end = path.start + i + node_len + 1;
                    tx_.send(format!("{}\t{}\t{}\t{}", path.name, start, end, edge_name))
                        .unwrap();
                }
                i += *node_lengths.get(&edge.1).unwrap();
            }
        });

    drop(tx);
    writer.join().unwrap();
}
