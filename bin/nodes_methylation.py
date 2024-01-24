#!/usr/bin/env pypy3

import pickle
import sys

nodes_list_path = sys.argv[1]
pickle_mc_in = sys.argv[2]
# pickle_cov_in = sys.argv[3]

mc = pickle.load(open(pickle_mc_in, "rb"))
# mcov = pickle.load(open(pickle_cov_in, "rb"))

nodes_list = set()
with open(nodes_list_path, "r") as f:
    for line in f:
        nodes_list.add(line.rstrip())

for node in nodes_list:
    # check if there is methylation signal on this node
    if node in mc and len(mc[node]) > 0:
        for pos in mc[node]:
            strand = "+"
            if pos < 0:
                strand = "-"

            print(node, abs(pos), strand, mc[node][pos][0], mc[node][pos][1]/mc[node][pos][0])
