from Bio import Phylo
import ete3
import bte
import pandas as pd
import sys
import time

travdf = {k:[] for k in ['Package','TreeSize','PathTime']}
mat = bte.MATree(sys.argv[1])
ts = int(sys.argv[1].split("/")[1].split("_")[0])
#pick a random leaf to rsearch from on the tree.
target = mat.get_random(1).get_leaves_ids()[0]

start = time.perf_counter()
mat.rsearch(target)
end = time.perf_counter()
travdf['Package'].append("BTE")
travdf['TreeSize'].append(ts)
travdf['PathTime'].append(end-start)

nwk = mat.write_newick()
nwk_f = open("subtree.nwk","w+")
print(nwk,file=nwk_f)
nwk_f.close()

etet = ete3.Tree(nwk,format=1)
ete_target = etet.search_nodes(name=target)[0]
start = time.perf_counter()
p = []
while ete_target.up:
    ete_target = ete_target.up
    p.append(ete_target)
end = time.perf_counter()
travdf['Package'].append("ETE")
travdf['TreeSize'].append(ts)
travdf['PathTime'].append(end-start)

phyt = Phylo.read("subtree.nwk", "newick")
phyt_target = list(phyt.find_clades(name=target))[0]
start = time.perf_counter()
phyt.get_path(phyt_target)
end = time.perf_counter()
travdf['Package'].append("BioPhylo")
travdf['TreeSize'].append(ts)
travdf['PathTime'].append(end-start)

travdf = pd.DataFrame(travdf)
travdf.to_csv("rsearch_times.csv",index=False)
