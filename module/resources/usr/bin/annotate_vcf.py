#!/usr/bin/env python3
import sys
import vcfpy
import re
import gzip


def parse_path_re(s):
    nodes = []
    matches = re.finditer(r"([<>])([a-zA-Z]*[0-9]+)", s)
    for m in matches:
        nodes.append((m.group(1), m.group(2)))
    return nodes


vcf = sys.argv[1]
mods = gzip.open(sys.argv[2])
mods.readline()
out_vcf = sys.argv[3]

node_mods = {}

for line in mods:
    node, pos, strand, depth, level = line.decode().split()
    pos = int(pos)
    if node not in node_mods:
        node_mods[node] = []
    node_mods[node].append((pos, strand, depth, level))

mods.close()

reader = vcfpy.Reader.from_path(vcf)
writer = vcfpy.Writer.from_path(out_vcf, reader.header)

for record in reader:
    paths = {}
    genotype = record.calls[0].gt_alleles

    for g in genotype:
        if g is not None:
            paths[g] = parse_path_re(record.INFO['AT'][g])

    nodes_positions = {}
    nodes_strands = {}
    nodes_depths = {}
    nodes_levels = {}

    for g in paths:
        nodes_positions[g] = []
        nodes_strands[g] = []
        nodes_depths[g] = []
        nodes_levels[g] = []
        for node in paths[g]:
            if node in node_mods:
                nodes_positions[g].append(node + ":" + ",".join([n[0] for n in node_mods[node]]))
                nodes_strands[g].append(node + ":" + ",".join([n[1] for n in node_mods[node]]))
                nodes_depths[g].append(node + ":" + ",".join([n[2] for n in node_mods[node]]))
                nodes_levels[g].append(node + ":" + ",".join([n[3] for n in node_mods[node]]))

    record.INFO["PNP"] = "|".join([g + "=" + ";".join(nodes_positions[g]) for g in nodes_positions])
    record.INFO["PNS"] = "|".join([g + "=" + ";".join(nodes_strands[g]) for g in nodes_strands])
    record.INFO["PND"] = "|".join([g + "=" + ";".join(nodes_depths[g]) for g in nodes_depths])
    record.INFO["PNL"] = "|".join([g + "=" + ";".join(nodes_levels[g]) for g in nodes_levels])

    writer.write_record(record)

writer.close()
