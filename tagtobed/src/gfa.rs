use flate2::read::MultiGzDecoder;
use std::fs::File;
use std::io::{self, BufRead, BufReader};

pub fn open_reader(path: &str) -> io::Result<Box<dyn BufRead + Send>> {
    let mut reader = BufReader::new(File::open(path)?);
    if reader.fill_buf()?.starts_with(&[0x1f, 0x8b]) {
        Ok(Box::new(BufReader::new(MultiGzDecoder::new(reader))))
    } else {
        Ok(Box::new(reader))
    }
}

pub struct GraphPath {
    pub name: String,
    pub start: usize,
    pub nodes: Vec<(usize, &'static str)>,
}

pub fn parse_path(line: &str) -> Result<Option<GraphPath>, String> {
    let fields: Vec<&str> = line.split_whitespace().collect();
    match fields.first().copied() {
        Some("P") => {
            if fields.len() < 4 {
                return Err("P record requires at least four fields".to_string());
            }
            let (name, start) = match fields[1].split_once('[') {
                Some((name, interval)) => (
                    name.to_string(),
                    interval
                        .split('-')
                        .next()
                        .unwrap()
                        .parse::<usize>()
                        .map_err(|_| "Invalid P record start coordinate".to_string())?,
                ),
                None => (fields[1].to_string(), 0),
            };
            let nodes = fields[2]
                .split(',')
                .map(|node| {
                    let (id, strand) = if let Some(id) = node.strip_suffix('+') {
                        (id, "+")
                    } else if let Some(id) = node.strip_suffix('-') {
                        (id, "-")
                    } else {
                        return Err(format!("Invalid P record segment: {node}"));
                    };
                    Ok((
                        id.parse::<usize>()
                            .map_err(|_| format!("Invalid segment ID: {id}"))?,
                        strand,
                    ))
                })
                .collect::<Result<Vec<_>, String>>()?;
            Ok(Some(GraphPath { name, start, nodes }))
        }
        Some("W") => {
            if fields.len() < 7 {
                return Err("W record requires at least seven fields".to_string());
            }
            let name = format!("{}#{}#{}", fields[1], fields[2], fields[3]);
            // An unknown sequence start gives walk-relative coordinates.
            let start = if fields[4] == "*" {
                0
            } else {
                fields[4]
                    .parse::<usize>()
                    .map_err(|_| "Invalid W record start coordinate".to_string())?
            };
            let mut walk = fields[6];
            let mut nodes = Vec::new();
            while !walk.is_empty() {
                let strand = match walk.as_bytes()[0] {
                    b'>' => "+",
                    b'<' => "-",
                    _ => return Err(format!("Invalid W record walk: {walk}")),
                };
                walk = &walk[1..];
                let end = walk.find(['>', '<']).unwrap_or(walk.len());
                let id = &walk[..end];
                nodes.push((
                    id.parse::<usize>()
                        .map_err(|_| format!("Invalid segment ID: {id}"))?,
                    strand,
                ));
                walk = &walk[end..];
            }
            Ok(Some(GraphPath { name, start, nodes }))
        }
        _ => Ok(None),
    }
}
