import sys


flip = {'+' : '-', '-' : '+'}
def complement(edge):
    return (flip[edge[2]], edge[3], flip[edge[0]], edge[1])

cpg_index = set()

strands = {'>' : '+', '<' : '-'}
strands_ = {'+' : '>', '-' : '<'}

with open(sys.argv[1], 'r') as cpgs_file:
    for line in cpgs_file:
        left, right = line.rstrip().split()
        lstrand = strands[left[0]]
        rstrand = strands[right[0]]

        lnode = int(left[1:])
        rnode = int(right[1:])

        edge = (lstrand, lnode, rstrand, rnode)
        cpg_index.add(edge)

# index node lengths

node_lengths = {}

gfa = open(sys.argv[2], 'r')
for line in gfa:
    if line[0] != "S":
        continue

    _, node_id, node_seq = line.rstrip().split()

    node_lengths[int(node_id)] = len(node_seq)

gfa.seek(0)

for line in gfa:
    if line[0] != "P":
        continue
    _, p_name, hap, _ = line.rstrip().split()
    hap_name = p_name
    hap_start = 0

    hap_seq =  hap.split(',')
    hap_edges = []
    for i in range(0, len(hap_seq) - 1):
        lstrand = hap_seq[i][-1]
        rstrand = hap_seq[i+1][-1]

        lnode = int(hap_seq[i][:-1])
        rnode = int(hap_seq[i+1][:-1])

        edge = (lstrand, lnode, rstrand, rnode)

        hap_edges.append(edge)

    # unless we have clipped contig
    if '[' in p_name:
        hap_name = p_name.split('[')[0]
        hap_start = int(p_name.split('[')[1].split('-')[0])

    i = 0
    for edge in hap_edges:
        cedge = complement(edge)
        edge_name = strands_[edge[0]] + str(edge[1]) + ' ' + strands_[edge[2]] + str(edge[3])

        if cedge in cpg_index:
            edge_name = strands_[cedge[0]] + str(cedge[1]) + ' ' + strands_[cedge[2]] + str(cedge[3])

        if edge in cpg_index or cedge in cpg_index:
            print(hap_name, hap_start + i + node_lengths[edge[1]] - 1, hap_start + i + node_lengths[edge[1]] + 1 , edge_name, sep = '\t')
        # move the position past the edge
        i += node_lengths[edge[1]]

gfa.close()
