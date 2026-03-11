#!/opt/pypy3/bin/pypy3

import pickle
import sys
import gzip

pickle_mc_in = sys.argv[1]
mc = pickle.load(open(pickle_mc_in, "rb"))
cpgs_index = gzip.open(sys.argv[2], "r")

for line in cpgs_index:
    node, pos, strand, pair = line.decode().rstrip().split('\t')

    pos = int(pos)

    depth = 0
    score = 0

    if node in mc and pos in mc[node]:
        depth = mc[node][pos][0]
        score = mc[node][pos][1]/mc[node][pos][0]
    print(node, abs(pos), strand, depth, score, pair)
