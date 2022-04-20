from Bio import Phylo
import ete3
import bte
import pandas as pd
import sys
import time

travdf = {k:[] for k in ['Package','TreeSize','FindTime']}
mat = bte.MATree(sys.argv[1])
ts = int(sys.argv[1].split("/")[1].split("_")[0])
#pick a random leaf to find from the tree.
target = mat.get_random(1).get_leaves_ids()[0]

start = time.perf_counter()
mat.get_node(target)
end = time.perf_counter()
travdf['Package'].append("BTE")
travdf['TreeSize'].append(ts)
travdf['FindTime'].append(end-start)

nwk = mat.write_newick()
nwk_f = open("subtree.nwk","w+")
print(nwk,file=nwk_f)
nwk_f.close()

etet = ete3.Tree(nwk,format=1)
start = time.perf_counter()
etet.search_nodes(name=target)
end = time.perf_counter()
travdf['Package'].append("ETE")
travdf['TreeSize'].append(ts)
travdf['FindTime'].append(end-start)

phyt = Phylo.read("subtree.nwk", "newick")
start = time.perf_counter()
phyt.find_clades(name=target)
end = time.perf_counter()
travdf['Package'].append("BioPhylo")
travdf['TreeSize'].append(ts)
travdf['FindTime'].append(end-start)

travdf = pd.DataFrame(travdf)
travdf.to_csv("find_times.csv",index=False)
