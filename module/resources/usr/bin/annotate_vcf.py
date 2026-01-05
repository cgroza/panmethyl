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

def path_to_str(path):
    return "".join([node[0] + node[1] for node in path])

def drop_source_sink(path):
    return path[1:][:-1]

vcf = sys.argv[1]
mods = gzip.open(sys.argv[2])
mods.readline()
out_vcf = sys.argv[3]

node_mods = {}

print("CHROM", "POS", "ID", "allele", "path", "node_positions",
      "strands", "depths", "levels", sep='\t')

for line in mods:
    node, pos, strand, depth, level = line.decode().split()
    if node not in node_mods:
        node_mods[node] = []
    # skip nucleotides without coverage
    if float(depth) <= 0:
        continue
    node_mods[node].append((int(pos), strand, float(depth), float(level)))

mods.close()

reader = vcfpy.Reader.from_path(vcf)
reader.header.add_format_line({'ID' : 'PMN', 'source' : 'panmethyl', 'Type' : 'Float', 'Description' : 'Number of modified nucleotides with coverage on the path.', 'Number' : '.'})
reader.header.add_format_line({'ID' : 'PML', 'source' : 'panmethyl', 'Type' : 'Float', 'Description' : 'Average modification levels across the the path.', 'Number' : '.'})
reader.header.add_format_line({'ID' : 'PMD', 'source' : 'panmethyl', 'Type' : 'Float', 'Description' : 'Average depth across the modified nucleotides.', 'Number' : '.'})

writer = vcfpy.Writer.from_path(out_vcf, reader.header)

for record in reader:
    paths = {}
    genotype = record.calls[0].gt_alleles

    for g in genotype:
        if g is not None:
            paths[g] = drop_source_sink(parse_path_re(record.INFO['AT'][g]))

    pmn = {}
    pml = {}
    pmd = {}
    for g in paths:
        nodes_positions = []
        nodes_strands = []
        nodes_depths = []
        nodes_levels = []

        pmn[g] = 0
        pml[g] = 0
        pmd[g] = 0
        for node in paths[g]:
            if node[1] in node_mods:
                nodes_positions.append(node[1] + ":" + ",".join([str(n[0]) for n in node_mods[node[1]]]))
                nodes_strands.append(",".join([n[1] for n in node_mods[node[1]]]))
                nodes_depths.append(",".join([str(n[2]) for n in node_mods[node[1]]]))
                nodes_levels.append(",".join([str(n[3]) for n in node_mods[node[1]]]))
        print(record.CHROM, record.POS, record.ID[0], str(g),
              path_to_str(paths[g]),
              ";".join(nodes_positions),
              ",".join(nodes_strands),
              ",".join(nodes_depths),
              ",".join(nodes_levels),
              sep='\t')

        pmn[g] = len([float(d) for node in nodes_depths for d in node.split(',')])
        if pmn[g] != 0:
            pml[g] = sum([float(l) for node in nodes_levels for l in node.split(',')])/pmn[g]
            pmd[g] = sum([float(d) for node in nodes_depths for d in node.split(',')])/pmn[g]

    record.add_format("PMN", [str(pmn[g]) for g in pmn])
    record.add_format("PML", [str(pml[g]) for g in pml])
    record.add_format("PMD", [str(pmd[g]) for g in pmd])
    writer.write_record(record)

writer.close()
