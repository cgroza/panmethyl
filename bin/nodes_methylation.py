#!/usr/bin/env pypy3

import pickle
import sys
import gzip

nodes_list_path = sys.argv[1]
pickle_mc_in = sys.argv[2]
mc = pickle.load(open(pickle_mc_in, "rb"))
cpgs_index = gzip.open(sys.argv[3], "r")

cpgs = dict()

for line in cpgs_index:
    node, pos, strand, pair = line.decode().split()
    pos = int(pos)
    if not node in cpgs:
      cpgs[node] = dict()

    if not pos in cpgs[node]:
        if strand == "+":
            cpgs[node][pos] = 0
        elif strand == "-":
            cpgs[node][-pos] = 0

nodes_list = set()
with open(nodes_list_path, "r") as f:
    for line in f:
        nodes_list.add(line.rstrip())

for node in nodes_list:
    # check if there is cpg on this node
    if node in cpgs: #and len(mc[node]) > 0:
        for pos in cpgs[node]:
            strand = "+"
            if pos < 0:
                strand = "-"
            depth = 0
            score = 0
            if node in mc and pos in mc[node]:
                depth = mc[node][pos][0]
                score = mc[node][pos][1]/mc[node][pos][0]
            print(node, abs(pos), strand, depth, score)
