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

node_sizes_path = sys.argv[4]

node_sizes = {}
with open(node_sizes_path) as node_sizes_file:
    for line  in node_sizes_file:
        node_name, node_size = line.split()
        node_sizes[node_name] = int(node_size)


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

    record['pstart'] = int(record['pstart'])
    record['pend'] = int(record['pend'])
    record['plen'] = int(record['plen'])

    pml = []
    pmd = []

    if len(nodes) == 1 and nodes[0][1] in node_mods:
        # flip reversed nodes
        if nodes[0][0] == "<":
            record['pstart'] = node_sizes[nodes[0][1]] - record['pend']
            record['pend'] = node_sizes[nodes[0][1]] - record['pstart']

        pmd = pmd + [n[2] for n in node_mods[nodes[0][1]] if n[0] >= record['pstart'] and n[0] <= record['pend']]
        pml = pml + [n[3] for n in node_mods[nodes[0][1]] if n[0] >= record['pstart'] and n[0] <= record['pend']]

    else:
        if nodes[0][1] in node_mods:
            if node[0][0] == ">":
                pmd = pmd + [n[2] for n in node_mods[nodes[0][1]] if n[0] >= record['pstart']]
                pml = pml + [n[3] for n in node_mods[nodes[0][1]] if n[0] >= record['pstart']]
            # flip offsets on first node
            else:
                pmd = pmd + [n[2] for n in node_mods[nodes[0][1]] if node_sizes[nodes[0][1]] - n[0] >= record['pstart']]
                pml = pml + [n[3] for n in node_mods[nodes[0][1]] if node_sizes[nodes[0][1]] - n[0] >= record['pstart']]

        if nodes[-1][1] in node_mods:
            # where does the last node on the path begin
            last_node_start = record['plen'] - node_sizes[nodes[-1][1]]
            if nodes[-1][0] == ">":
                pmd = pmd + [n[2] for n in node_mods[nodes[-1][1]] if n[0] + last_node_start <= record['pend']]
                pml = pml + [n[3] for n in node_mods[nodes[-1][1]] if n[0] + last_node_start <= record['pend']]
            # flip offsets on last node
            else:
                pmd = pmd + [n[2] for n in node_mods[nodes[-1][1]] if (node_sizes[nodes[-1][1]] - n[0]) + last_node_start <= record['pend']]
                pml = pml + [n[3] for n in node_mods[nodes[-1][1]] if (node_sizes[nodes[-1][1]] - n[0]) + last_node_start <= record['pend']]

        for node in nodes[1:-1]:
            if node[1] in node_mods:
                pmd = pmd + [n[2] for n in node_mods[node[1]]]
                pml = pml + [n[3] for n in node_mods[node[1]]]

    pmn = len(pmd)

    if pmn == 0:
        print(record['qname'], record['pname'], "PML:PMD:PMN", ".:.:.", sep = '\t')
        continue

    print(record['qname'], record['pname'], "PML:PMD:PMN", ":".join(map(str, [round(sum(pml)/pmn, 2), round(sum(pmd)/pmn, 2), pmn])), sep = '\t')

