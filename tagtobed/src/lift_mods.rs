// lift.rs

use std::collections::HashMap;
use std::env;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};
use regex::Regex;
use csv::Writer;

// CIGAR parsing from cs:Z
fn cs2cg(cs: &str) -> String {
    let re = Regex::new(r"(:[0-9]+)|(\+[ACTGNactg]+)|(-[ACTGNactg]+)|(\*[ACTGNactg]{2})").unwrap();
    let mut cigar_ops = Vec::new();
    let mut n_mismatch = 0;
    for cap in re.captures_iter(cs) {
        if n_mismatch > 0 && cap.get(4).is_none() {
            cigar_ops.push(format!("{}X", n_mismatch));
            n_mismatch = 0;
        }
        if let Some(m) = cap.get(1) {
            cigar_ops.push(format!("{}=", &m.as_str()[1..]));
            continue;
        }
        if let Some(i) = cap.get(2) {
            cigar_ops.push(format!("{}I", i.as_str()[1..].len()));
            continue;
        }
        if let Some(d) = cap.get(3) {
            cigar_ops.push(format!("{}D", d.as_str()[1..].len()));
            continue;
        }
        if let Some(_x) = cap.get(4) {
            n_mismatch += 1;
            continue;
        }
    }
    cigar_ops.join("")
}

fn parse_path_re(s: &str) -> Vec<(String, String)> {
    let re = Regex::new(r"([<>])([a-zA-Z]*[0-9]+)").unwrap();
    let mut nodes = Vec::new();
    for cap in re.captures_iter(s) {
        nodes.push((cap[1].to_string(), cap[2].to_string()));
    }
    nodes
}

fn increment_count(
    d: &mut HashMap<String, HashMap<i32, (i32, i32)>>,
    node: &str,
    pos: i32,
    ml_score: i32,
) {
    if ml_score == -1 {
        return;
    }
    let score = if ml_score < 128 { 0 } else { 1 };
    let entry = d.entry(node.to_string()).or_insert_with(HashMap::new);
    entry.entry(pos).and_modify(|v| {
        v.0 += 1;
        v.1 += score;
    }).or_insert((1, score));
}

fn main() -> io::Result<()> {
    let args: Vec<String> = env::args().collect();
    if args.len() < 3 {
        eprintln!("Usage: lift <node_sizes_path> <csv_out>");
        std::process::exit(1);
    }
    let node_sizes_path = &args[1];
    let csv_out = &args[2];

    // Load node sizes
    let mut node_sizes = HashMap::new();
    let node_sizes_file = File::open(node_sizes_path)?;
    for line in BufReader::new(node_sizes_file).lines() {
        let line = line?;
        let parts: Vec<&str> = line.split_whitespace().collect();
        node_sizes.insert(parts[0].to_string(), parts[1].parse::<usize>().unwrap());
    }

    let mut node_bmod_count: HashMap<String, HashMap<i32, (i32, i32)>> = HashMap::new();
    for node in node_sizes.keys() {
        node_bmod_count.insert(node.clone(), HashMap::new());
    }

    let stdin = io::stdin();
    let mut nrecords = 0;
    for line in stdin.lock().lines() {
        let line = line?;
        if nrecords % 1000 == 0 {
            eprintln!("Parsed {} records", nrecords);
        }
        nrecords += 1;
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() < 5 || fields[4] == "*" {
            continue;
        }
        assert_eq!(fields[4], "+");
        if fields.len() < 14 {
            continue;
        }

        let mapq = fields[11].parse::<usize>().unwrap();
        if mapq < 1 {
            continue
        }

        let bmods: Vec<&str> = fields.last().unwrap().split(',').collect();
        let qlen = fields[1].parse::<usize>().unwrap();
        let qstart = fields[2].parse::<usize>().unwrap();
        let qend = fields[3].parse::<usize>().unwrap();
        let path = fields[5];
        let plen = fields[6].parse::<usize>().unwrap();
        let pstart = fields[7].parse::<usize>().unwrap();
        let pend = fields[8].parse::<usize>().unwrap();

        let cg_str = if fields[12].starts_with("cg:Z:") {
            &fields[12][5..]
        } else if fields[12].starts_with("cs:Z:") {
            &cs2cg(&fields[12][5..])
        } else {
            panic!("Unknown alignment tag");
        };

        // Parse CIGAR string
        let mut cg_items = Vec::new();
        let re_cigar = Regex::new(r"(\d+)([XID=])").unwrap();
        for cap in re_cigar.captures_iter(cg_str) {
            let len = cap[1].parse::<usize>().unwrap();
            let op = cap[2].chars().next().unwrap();
            cg_items.push((len, op));
        }

        let parsed_path = parse_path_re(path);

        let mut pmatches = Vec::new();
        let mut qmatches = Vec::new();
        let mut qi = qstart;
        let mut pi = pstart;
        let mut nmatch = 0;
        let mut nmismatch = 0;
        for op in &cg_items {
            match op.1 {
                'X' | '=' => {
                    if op.1 == 'X' {
                        nmismatch += op.0;
                    } else {
                        nmatch += op.0;
                    }
                    pmatches.push((pi as i32, (pi + op.0) as i32));
                    qmatches.push((qi as i32, (qi + op.0) as i32));
                    qi += op.0;
                    pi += op.0;
                }
                'I' => {
                    qi += op.0;
                }
                'D' => {
                    pi += op.0;
                }
                _ => {
                    panic!("Unknown CIGAR operation {}", op.1);
                }
            }
        }

        let mut bmods_path = Vec::new();
        for bmod in bmods {
            let parts: Vec<&str> = bmod.split(':').collect();
            if parts.len() != 2 { continue; }
            let b = parts[0].parse::<i32>().unwrap();
            let bscore = parts[1].parse::<i32>().unwrap();
            for (i, qmatch) in qmatches.iter().enumerate() {
                if b >= qmatch.0 && b < qmatch.1 {
                    let offset = b - qmatch.0;
                    let pmatch = pmatches[i];
                    if pmatch.0 + offset >= pmatch.1 { continue; }
                    bmods_path.push((pmatch.0 + offset, bscore));
                }
            }
        }

        let mut node_breaks = Vec::new();
        let mut j = 0;
        for node in &parsed_path {
            let node_size = node_sizes[&node.1];
            node_breaks.push((j as i32, (j + node_size) as i32, node.clone()));
            j += node_size;
        }

        for (pb, pbscore) in bmods_path {
            for nb in &node_breaks {
                if pb >= nb.0 && pb < nb.1 {
                    let node = &nb.2;
                    let node_size = node_sizes[&node.1];
                    let mut offset = pb - nb.0;
                    if node.0 == "<" {
                        offset = -(node_size as i32 - offset - 1);
                    }
                    increment_count(&mut node_bmod_count, &node.1, offset, pbscore);
                    break;
                }
            }
        }
    }

    let mut wtr = Writer::from_path(csv_out)?;
    // wtr.write_record(&["node", "position", "count", "score_sum"])?;
    for (node, pos_map) in node_bmod_count.iter() {
        for (pos, (count, score_sum)) in pos_map.iter() {
            wtr.write_record(&[
                node,
                &pos.to_string(),
                &count.to_string(),
                &score_sum.to_string(),
            ])?;
        }
    }
    wtr.flush()?;
    Ok(())
}
