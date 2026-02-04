#!/usr/bin/env python3
import sys
import re
import gzip

def parse_path_re(s):
    nodes = []
    matches = re.finditer(r"([<>])([a-zA-Z]*[0-9]+)", s)
    for m in matches:
        nodes.append((m.group(1), m.group(2)))
    return nodes

def path_to_str(path):
    return "".join([node[0] + node[1] for node in path])

gaf_path = sys.argv[1]
mods = gzip.open(sys.argv[2])
mods.readline()
sample_name = sys.argv[3]

node_mods = {}

for line in mods:
    node, pos, strand, depth, level = line.decode().split()
    if node not in node_mods:
        node_mods[node] = []
    # skip nucleotides without coverage
    if float(depth) <= 0:
        continue
    node_mods[node].append((int(pos), strand, float(depth), float(level)))

mods.close()

gaf = open(gaf_path)
# skip header
gaf.readline()

fields = ["qname", "qlen", "qstart", "qend", "strand", "pname", "plen", "pstart", "pend", "matches", "alnblen", "mapq", "cs"]

for line in gaf:
    record = dict(zip(fields, line.rstrip().split()))
    nodes = parse_path_re(record["pname"])

    pml = []
    pmd = []

    for node in nodes:
        if node[1] in node_mods:
            pmd = pmd + [n[2] for n in node_mods[node[1]]]
            pml = pml + [n[3] for n in node_mods[node[1]]]
    pmn = len(pmd)

    if pmn == 0:
        print(record['qname'], record['pname'], "PML:PMD:PMN", ".:.:.", sep = '\t')
        continue

    print(record['qname'], record['pname'], "PML:PMD:PMN", ":".join(map(str, [round(sum(pml)/pmn, 2), round(sum(pmd)/pmn, 2), pmn])), sep = '\t')

