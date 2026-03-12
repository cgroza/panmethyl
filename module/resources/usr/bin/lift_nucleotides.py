#!/opt/pypy3/bin/pypy3
import sys
import gzip
import threading
from concurrent.futures import ThreadPoolExecutor
from queue import Queue

nuc_index = dict()


num_workers = int(sys.argv[3])

with gzip.open(sys.argv[1], 'rt', encoding = 'ascii') as cpgs_file:
    for line in cpgs_file:
        node, offset, strand, nucset = line.rstrip().split('\t')
        if ' ' in nucset:
            continue
        if strand == '-':
            continue
        if int(node) not in nuc_index:
            nuc_index[int(node)] = []
        nuc_index[int(node)].append((int(offset), strand, nucset))

# index node lengths

node_lengths = {}

gfa = open(sys.argv[2], 'r')
for line in gfa:
    if line[0] != "S":
        continue

    _, node_id, node_seq = line.rstrip().split()

    node_lengths[int(node_id)] = len(node_seq)

gfa.seek(0)

def process_line(line, nuc_index, node_lengths, output_queue):
    if line[0] != "P":
        return
    _, p_name, hap, _ = line.rstrip().split()
    hap_name = p_name
    hap_start = 0

    hap_seq = hap.split(',')

    if '[' in p_name:
        hap_name = p_name.split('[')[0]
        hap_start = int(p_name.split('[')[1].split('-')[0])

    i = 0
    for node in hap_seq:
        node_name = int(node[:-1])
        strand = node[-1]

        if node_name in nuc_index:
            for cpg in nuc_index[node_name]:
                offset = cpg[0]
                if strand == '-':
                    offset = node_lengths[node_name] - offset - 1
                output_queue.put((
                    hap_name,
                    hap_start + i + offset,
                    hap_start + i + offset + 2,
                    f"{node_name},{cpg[0]},{cpg[1]}",
                    cpg[2]
                ))
        i += node_lengths[node_name]

def output_worker(output_queue):
    while True:
        item = output_queue.get()
        if item is None:
            break
        print(*item, sep='\t')
        output_queue.task_done()

output_queue = Queue()

# Start output thread
output_thread = threading.Thread(target=output_worker, args=(output_queue,))
output_thread.start()

with ThreadPoolExecutor(max_workers=num_workers) as executor:
    for line in gfa:
        executor.submit(process_line, line, nuc_index, node_lengths, output_queue)

# Signal output thread to finish
output_queue.put(None)
output_thread.join()
gfa.close()
