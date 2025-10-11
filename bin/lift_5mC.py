#!/usr/bin/env pypy3

import cigar
import sys
import pickle

node_sizes_path = sys.argv[1]
pickle_out = sys.argv[2]

def parse_path(s):
    node_pos = 1
    nodes = []
    while node_pos > -1:
        next_node_pos = s.find("s", node_pos + 1)
        node_name = None
        if next_node_pos == -1:
            node_name = s[node_pos:]
        else:
            node_name = s[node_pos:next_node_pos - 1]
        nodes.append((s[node_pos - 1], node_name))
        node_pos = next_node_pos
    return(nodes)

def increment_count(d, node, pos, score):
    if pos in d[node]:
        d[node][pos] = (d[node][pos][0] + 1, d[node][pos][1] + score)
    else:
        d[node][pos] = (1, score)

pb_aln = sys.stdin

# load node sizes
node_sizes = {}
with open(node_sizes_path) as node_sizes_file:
    for line  in node_sizes_file:
        node_name, node_size = line.split()
        node_sizes[node_name] = int(node_size)

node_bmod_count = {}
# note, unique dict for each node
for node in node_sizes:
    node_bmod_count[node] = {}

nrecords = 0
for line in pb_aln:
    if nrecords % 1000 == 0:
        print("Parsed", nrecords, "records", file=sys.stderr)

    nrecords = nrecords + 1

    fields = line.split()

    assert fields[4] == "+"

    # this read has no modified bases
    if len(fields) < 18:
        continue

    # base:score pairs
    bmods = list(fields[-1].split(','))

    qlen = int(fields[1])
    qstart = int(fields[2])
    qend = int(fields[3])

    path = fields[5]
    plen = int(fields[6])
    pstart = int(fields[7])
    pend = int(fields[8])

    # cigar string
    assert fields[16].startswith("cg:Z:")
    cg = cigar.Cigar(fields[16][5:])

    parsed_path = parse_path(path)

    # list of intervals describing matches on path, left closed, right closed intervals
    pmatches = []
    # list of intervals describing matches on read, left closed, right closed intervals
    qmatches = []

    # record matches along path length
    # cursor to keep track of position in strings
    qi = qstart
    pi = pstart
    nmatch = 0
    nmismatch = 0
    for op in cg.items():
        match op[1]:
            case 'X' | '=':
                # count matches vs mismatches
                if op[1] == "X":
                    nmismatch = nmismatch + op[0]
                else:
                    nmatch = nmatch + op[0]
                pmatches.append((pi, pi + op[0]))
                qmatches.append((qi, qi + op[0]))
                # move the cursors
                qi = qi + op[0]
                pi = pi + op[0]
            # skip on query
            case 'I':
                qi = qi + op[0]
            # skip on path
            case 'D':
                pi = pi + op[0]
            case _:
                raise Exception("Unknown CIGAR operation " + op[1])

    # find base modifications
    bmods_path = []
    for bmod in bmods:
        [b, bscore] = bmod.split(":")
        b = int(b)
        bscore = int(bscore)

        i = 0
        for i in range(len(qmatches)):
            if b >= qmatches[i][0] and b < qmatches[i][1]:
                offset = b - qmatches[i][0]
                assert pmatches[i][0] + offset < pmatches[i][1]
                bmods_path.append((pmatches[i][0] + offset, bscore))

    node_breaks = []
    j = 0
    for node in parsed_path:
        node_breaks.append((j, j + node_sizes[node[1]], node))
        j = j + node_sizes[node[1]]

    for (pb, pbscore) in bmods_path:
        for nb in node_breaks:
            if pb >= nb[0] and pb < nb[1]:
                node = nb[2]
                offset = pb - nb[0]
                # normalize position for the forward strand, but mark with negative sign
                if node[0] == "<":
                    offset = -(node_sizes[node[1]] - offset - 1)  # since interval open on right

                increment_count(node_bmod_count, node[1], offset, pbscore)
                break           # there is only one possible matching node

with open(pickle_out, "wb") as f:
    pickle.dump(node_bmod_count, f)
