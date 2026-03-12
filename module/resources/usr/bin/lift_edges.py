#!/opt/pypy3/bin/pypy3
import sys
import threading
from concurrent.futures import ThreadPoolExecutor
from queue import Queue

flip = {'+' : '-', '-' : '+'}
def complement(edge):
    return (flip[edge[2]], edge[3], flip[edge[0]], edge[1])

edge_index = set()

strands = {'>' : '+', '<' : '-'}
strands_ = {'+' : '>', '-' : '<'}

num_workers = int(sys.argv[3])

with open(sys.argv[1], 'r') as cpgs_file:
    for line in cpgs_file:
        node, offset, strand, nucset = line.rstrip().split('\t')
        if ' ' not in nucset:
            continue
        left, right = nucset.split(' ')
        lstrand = strands[left[0]]
        rstrand = strands[right[0]]

        lnode = int(left[1:])
        rnode = int(right[1:])

        edge = (lstrand, lnode, rstrand, rnode)
        edge_index.add(edge)

# index node lengths

node_lengths = {}

gfa = open(sys.argv[2], 'r')
for line in gfa:
    if line[0] != "S":
        continue

    _, node_id, node_seq = line.rstrip().split()

    node_lengths[int(node_id)] = len(node_seq)

gfa.seek(0)

def process_edge_line(line, node_lengths, edge_index, strands_, complement, output_queue):
    if line[0] != "P":
        return
    _, p_name, hap, _ = line.rstrip().split()
    hap_name = p_name
    hap_start = 0

    hap_seq = hap.split(',')
    hap_edges = []
    for i in range(0, len(hap_seq) - 1):
        lstrand = hap_seq[i][-1]
        rstrand = hap_seq[i+1][-1]

        lnode = int(hap_seq[i][:-1])
        rnode = int(hap_seq[i+1][:-1])

        edge = (lstrand, lnode, rstrand, rnode)
        hap_edges.append(edge)

    if '[' in p_name:
        hap_name = p_name.split('[')[0]
        hap_start = int(p_name.split('[')[1].split('-')[0])

    i = 0
    for edge in hap_edges:
        cedge = complement(edge)
        edge_name = strands_[edge[0]] + str(edge[1]) + ' ' + strands_[edge[2]] + str(edge[3])

        if cedge in edge_index:
            edge_name = strands_[cedge[0]] + str(cedge[1]) + ' ' + strands_[cedge[2]] + str(cedge[3])

        if edge in edge_index or cedge in edge_index:
            output_queue.put((
                hap_name,
                hap_start + i + node_lengths[edge[1]] - 1,
                hap_start + i + node_lengths[edge[1]] + 1,
                edge_name
            ))
        i += node_lengths[edge[1]]

def output_worker(output_queue):
    while True:
        item = output_queue.get()
        if item is None:
            break
        print(*item, sep='\t')
        output_queue.task_done()

output_queue = Queue()

output_thread = threading.Thread(target=output_worker, args=(output_queue,))
output_thread.start()

with ThreadPoolExecutor(max_workers=num_workers) as executor:
    for line in gfa:
        executor.submit(
            process_edge_line,
            line,
            node_lengths,
            edge_index,
            strands_,
            complement,
            output_queue
        )

output_queue.put(None)
output_thread.join()

gfa.close()
