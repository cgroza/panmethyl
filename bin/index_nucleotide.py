#!/usr/bin/env pypy3
import sys

gfa = open(sys.argv[1])
nuc = sys.argv[2]

i = 0

if nuc != 'CG':
    for s in gfa:
        fields = s.split()
        # not a segment
        if fields[0] != "S":
            continue

        name = fields[1]
        seq = fields[2]

        nuc_i = seq.find(nuc, 0)
        while nuc_i > -1:
            print(name, nuc_i, "+", "N" + str(i))
        i = i + 1

elif nuc == 'CG':
    node_ends = {}
    node_sizes = {}

    for s in gfa:
        fields = s.split()
        # not a segment
        if fields[0] != "S":
            continue

        name = fields[1]
        seq = fields[2]

        node_ends[name] = (seq[0], seq[-1])
        node_sizes[name] = len(seq)

        cpg_i = seq.find("CG", 0)
        while cpg_i > -1:
            print(name, cpg_i, "+", "N" + str(i))
            print(name, cpg_i + 1, "-", "N" + str(i))
            cpg_i = seq.find("CG", cpg_i + 2)
        i = i + 1

    gfa.seek(0)

    for line in gfa:
        match line.split():
            case ["L", left, "+", right, "+", *rest]:
                if node_ends[left][1] + node_ends[right][0] == "CG":
                    print(left, node_sizes[left] - 1, "+", "E" + str(i))
                    print(right, 0, "-", "E" + str(i))
            case ["L", left, "+", right, "-", *rest]:
                if node_ends[left][1] + node_ends[right][1] == "CC":
                    print(left, node_sizes[left] - 1, "+", "E" + str(i))
                    print(right, node_sizes[right] - 1, "-", "E" + str(i))
            case ["L", left, "-", right, "+", *rest]:
                if node_ends[left][0] + node_ends[right][0] == "GG":
                    print(left, 0, "+", "E" + str(i))
                    print(right, 0, "-", "E" + str(i))
            case ["L", left, "-", right, "-", *rest]:
                if node_ends[left][0] + node_ends[right][1] == "GC":
                    print(left, 0, "+", "E" + str(i))
                    print(right, node_ends[right] - 1, "-", "E" + str(i))

        i = i + 1

gfa.close()
