from Bio import Phylo
import ete3
import bte
import pandas as pd
import sys
import time

travdf = {k:[] for k in ['Package','TreeSize','TraversalTime']}
mat = bte.MATree(sys.argv[1])
ts = int(sys.argv[1].split("/")[1].split("_")[0])

start = time.perf_counter()
mat.depth_first_expansion()
end = time.perf_counter()
travdf['Package'].append("BTE")
travdf['TreeSize'].append(ts)
travdf['TraversalTime'].append(end-start)

nwk = mat.write_newick()
nwk_f = open("subtree.nwk","w+")
print(nwk,file=nwk_f)
nwk_f.close()

etet = ete3.Tree(nwk,format=1)
start = time.perf_counter()
for n in etet.traverse('preorder'):
    pass
end = time.perf_counter()
travdf['Package'].append("ETE")
travdf['TreeSize'].append(ts)
travdf['TraversalTime'].append(end-start)

phyt = Phylo.read("subtree.nwk", "newick")
start = time.perf_counter()
for n in phyt.find_clades():
    pass
end = time.perf_counter()
travdf['Package'].append("BioPhylo")
travdf['TreeSize'].append(ts)
travdf['TraversalTime'].append(end-start)

travdf = pd.DataFrame(travdf)
travdf.to_csv("traversal_times.csv",index=False)
