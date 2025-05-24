#!/usr/bin/env pypy3
import sys

gfa = open(sys.argv[1])

i = 0
for s in gfa:
    fields = s.split()
    # not a segment
    if fields[0] != "S":
        continue
    name = fields[1]
    seq = fields[2]
    cpg_i = seq.find("CG", 0)
    while cpg_i > -1:
        print(name, cpg_i, "+")
        print(name, cpg_i + 1, "-")
        cpg_i = seq.find("CG", cpg_i + 2)
    i = i + 1
