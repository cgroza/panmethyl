#!/opt/pypy3/bin/pypy3
import sys

gfa = open(sys.argv[1])
nuc = sys.argv[2]

i = 0

complement = {"A" : "T",
              "G" : "C",
              "T" : "A",
              "C" : "G"}


if len(nuc) == 1:
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
            nuc_i = seq.find(nuc, nuc_i + 1)
            i = i + 1

elif len(nuc) == 2:
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

        cpg_i = seq.find(nuc, 0)
        while cpg_i > -1:
            print(name, cpg_i, "+", "N" + str(i))
            print(name, cpg_i + 1, "-", "N" + str(i))
            cpg_i = seq.find(nuc, cpg_i + 2)
            i = i + 1

    gfa.seek(0)

    for line in gfa:
        match line.split():
            case ["L", left, "+", right, "+", *rest]:
                if node_ends[left][1] + node_ends[right][0] == nuc:
                    print(left, node_sizes[left] - 1, "+", "E" + str(i))
                    print(right, 0, "-", "E" + str(i))
            case ["L", left, "+", right, "-", *rest]:
                if node_ends[left][1] + node_ends[right][1] == nuc[0] + complement[nuc[1]]:
                    print(left, node_sizes[left] - 1, "+", "E" + str(i))
                    print(right, node_sizes[right] - 1, "-", "E" + str(i))
            case ["L", left, "-", right, "+", *rest]:
                if node_ends[left][0] + node_ends[right][0] == complement[nuc[0]] + nuc[1]:
                    print(left, 0, "+", "E" + str(i))
                    print(right, 0, "-", "E" + str(i))
            case ["L", left, "-", right, "-", *rest]:
                if node_ends[left][0] + node_ends[right][1] == complement[nuc[0]] + complement[nuc[1]]:
                    print(left, 0, "+", "E" + str(i))
                    print(right, node_ends[right] - 1, "-", "E" + str(i))

        i = i + 1
else:
    sys.stderr.write("Unsupported motif\n")
    sys.exit(1)
