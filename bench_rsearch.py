from Bio import Phylo
import ete3
import bte
import pandas as pd
import sys
import time
from profile_decor import named_profile


ts = int(sys.argv[1].split("/")[1].split("_")[0])

@named_profile("BTE\trsearch\t" + str(ts))
def btersearch(tree,target):
   tree.rsearch(target)
@named_profile("ETE\trsearch\t" + str(ts))
def etersearch(etet,target):
    ete_target = etet.search_nodes(name=target)[0]
    p = []
    while ete_target.up:
        ete_target = ete_target.up
        p.append(ete_target)
@named_profile("BioPhylo\trsearch\t" + str(ts))
def phyrsearch(phyt,target):
    phyt_target = list(phyt.find_clades(name=target))[0]
    phyt.get_path(phyt_target)

travdf = {k:[] for k in ['Package','TreeSize','PathTime']}
matpb = bte.MATree(sys.argv[1])
#remove mutations from the tree by writing a newick and reloading from the newick. This removes systematic bias against BTE, as it 
#maintains mutations in the tree structure and the other packages do not, meaning BTE has higher memory usage (because more data is stored)
#instead, if we load from a newick alone, no mutations will be included.
nwk = matpb.write_newick()
mat = bte.MATree(nwk_string = nwk)
#pick a random leaf to rsearch from on the tree.
target = mat.get_random(1).get_leaves_ids()[0]
btersearch(mat,target)
start = time.perf_counter()
end = time.perf_counter()
travdf['Package'].append("BTE")
travdf['TreeSize'].append(ts)
travdf['PathTime'].append(end-start)

nwk = mat.write_newick()
nwk_f = open("subtree.nwk","w+")
print(nwk,file=nwk_f)
nwk_f.close()

etet = ete3.Tree(nwk,format=1)
start = time.perf_counter()
etersearch(etet,target)
end = time.perf_counter()
travdf['Package'].append("ETE")
travdf['TreeSize'].append(ts)
travdf['PathTime'].append(end-start)

phyt = Phylo.read("subtree.nwk", "newick")
start = time.perf_counter()
phyrsearch(phyt,target)
end = time.perf_counter()
travdf['Package'].append("BioPhylo")
travdf['TreeSize'].append(ts)
travdf['PathTime'].append(end-start)

travdf = pd.DataFrame(travdf)
travdf.to_csv("rsearch_times.csv",index=False)
