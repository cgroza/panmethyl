#!/opt/pypy3/bin/pypy3
import sys
import gzip

mc_csv = gzip.open(sys.argv[1], "rt", encoding='ascii')

mc = dict()

for line in mc_csv:
    node, pos, depth, total_score = line.rstrip().split('\t')
    pos = abs(int(pos))
    depth = int(depth)
    total_score = float(total_score)

    if depth == 0:
        continue
    if node not in mc:
        mc[node] = {}

    mc[node][pos] = (depth, total_score/depth)

cpgs_index = gzip.open(sys.argv[2], "rt", encoding='ascii')

for line in cpgs_index:
    node, pos, strand, pair = line.rstrip().split('\t')

    pos = int(pos)

    depth = 0
    score = 0

    if node in mc and pos in mc[node]:
        depth = mc[node][pos][0]
        score = mc[node][pos][1]
    print(node, pos, strand, depth, score, pair, sep='\t')
