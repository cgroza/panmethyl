#!/opt/pypy3/bin/pypy3
import sys

gfa = open(sys.argv[1])
nuc = sys.argv[2]

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
            nuc_name = name + "," + str(nuc_i) + ",+"
            print(name, nuc_i, "+", nuc_name, sep='\t')
            nuc_i = seq.find(nuc, nuc_i + 1)

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

        dinuc_i = seq.find(nuc, 0)
        while dinuc_i > -1:
            dinuc_name = name + "," + str(dinuc_i) + ",+"
            print(name, dinuc_i, "+", dinuc_name, sep='\t')
            print(name, dinuc_i + 1, "-", dinuc_name, sep='\t')
            dinuc_i = seq.find(nuc, dinuc_i + 2)

    gfa.seek(0)

    for line in gfa:
        match line.split():
            case ["L", left, "+", right, "+", *rest]:
                if node_ends[left][1] + node_ends[right][0] == nuc:
                    name = ">" + left + " >" + right
                    print(left, node_sizes[left] - 1, "+", name, sep='\t')
                    print(right, 0, "-", name, sep='\t')
            case ["L", left, "+", right, "-", *rest]:
                if node_ends[left][1] + node_ends[right][1] == nuc[0] + complement[nuc[1]]:
                    name = ">" + left + " <" + right
                    print(left, node_sizes[left] - 1, "+", name, sep='\t')
                    print(right, node_sizes[right] - 1, "-", name, sep='\t')
            case ["L", left, "-", right, "+", *rest]:
                if node_ends[left][0] + node_ends[right][0] == complement[nuc[0]] + nuc[1]:
                    name = "<" + left + " >" + right
                    print(left, 0, "+", name, sep='\t')
                    print(right, 0, "-", name, sep='\t')
            case ["L", left, "-", right, "-", *rest]:
                if node_ends[left][0] + node_ends[right][1] == complement[nuc[0]] + complement[nuc[1]]:
                    name = "<" + left + " <" + right
                    print(left, 0, "+", name, sep='\t')
                    print(right, node_ends[right] - 1, "-", name, sep='\t')
else:
    sys.stderr.write("Unsupported motif\n")
    sys.exit(1)
