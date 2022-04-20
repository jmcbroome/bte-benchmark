from Bio import Phylo
import ete3
import bte
import pandas as pd
import sys
import time

subtdf = {k:[] for k in ['Package','TreeSize','SubtreeSize','SubtreeTime']}
mat = bte.MATree(sys.argv[1])
ts = int(sys.argv[1].split("/")[1].split("_")[0])
target_size = round(ts/10)
targets = mat.get_random(target_size).get_leaves_ids()
start = time.perf_counter()
mat.subtree(targets)
end = time.perf_counter()
subtdf['Package'].append("BTE")
subtdf['TreeSize'].append(ts)
subtdf['SubtreeSize'].append(target_size)
subtdf['SubtreeTime'].append(end-start)

nwk = mat.write_newick()
nwk_f = open("subtree.nwk","w+")
print(nwk,file=nwk_f)
nwk_f.close()

etet = ete3.Tree(nwk,format=1)
start = time.perf_counter()
#ete3 prune actually prunes everything BUT the input nodes.
etet.prune(targets)
end = time.perf_counter()
subtdf['Package'].append("ETE")
subtdf['TreeSize'].append(ts)
subtdf['SubtreeSize'].append(target_size)
subtdf['SubtreeTime'].append(end-start)

phyt = Phylo.read("subtree.nwk", "newick")
start = time.perf_counter()
#phylo does not have a single-line prune, so instead we have to manually prune all non-targets.
for t in phyt.get_terminals():
    if t.name not in targets:
        phyt.prune(t)
end = time.perf_counter()
subtdf['Package'].append("BioPhylo")
subtdf['TreeSize'].append(ts)
subtdf['SubtreeSize'].append(target_size)
subtdf['SubtreeTime'].append(end-start)

subtdf = pd.DataFrame(subtdf)
subtdf.to_csv("subtree_times.csv",index=False)