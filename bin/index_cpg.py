#!/usr/bin/env pypy3
import pickle
import sys

cpg_index = dict()

gfa = open(sys.argv[1])

i = 0
for s in gfa:
    fields = s.split()
    # not a segment
    if fields[0] != "S":
        continue
    name = fields[1]
    seq = fields[2]
#    if i % 1000 == 0:
#      print("Processed " + str(i) + " nodes")
    cpg_index[name] = set()
    cpg_i = seq.find("CG", 0)
    while cpg_i > -1:
        cpg_index[name].add(cpg_i)
        cpg_index[name].add(-(cpg_i + 1))
        cpg_i = seq.find("CG", cpg_i + 2)
    i = i + 1

# with open("gfa/cpg_index_xaf.pickle", "wb") as p:
#     pickle.dump(cpg_index, p)


# index = pickle.load(open("gfa/cpg_index_xaf.pickle", "rb"))

for node in cpg_index:
    for pos in cpg_index[node]:
        strand = "+"
        if pos < 0:
            strand = "-"
        print(node, abs(pos), strand)
