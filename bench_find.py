from Bio import Phylo
import ete3
import bte
import pandas as pd
import sys
import time
from profile_decor import named_profile

travdf = {k:[] for k in ['Package','TreeSize','FindTime']}
matpb = bte.MATree(sys.argv[1])
#remove mutations from the tree by writing a newick and reloading from the newick. This removes systematic bias against BTE, as it 
#maintains mutations in the tree structure and the other packages do not, meaning BTE has higher memory usage (because more data is stored)
#instead, if we load from a newick alone, no mutations will be included.
nwk = matpb.write_newick()
mat = bte.MATree(nwk_string = nwk)
ts = int(sys.argv[1].split("/")[1].split("_")[0])
#pick a random leaf to find from the tree.
target = mat.get_random(1).get_leaves_ids()[0]

@named_profile("BTE\tfind\t" + str(ts))
def btefind(tree,target):
    tree.get_node(target)

@named_profile("ETE\tfind\t" + str(ts))
def etefind(tree,target):
    tree.iter_search_nodes(name=target)

@named_profile("BioPhylo\tfind\t" + str(ts))
def phyfind(tree,target):
    tree.find_clades(name=target)

start = time.perf_counter()
btefind(mat,target)
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
etefind(etet,target)
end = time.perf_counter()
travdf['Package'].append("ETE")
travdf['TreeSize'].append(ts)
travdf['FindTime'].append(end-start)

phyt = Phylo.read("subtree.nwk", "newick")
start = time.perf_counter()
phyfind(phyt,target)
end = time.perf_counter()
travdf['Package'].append("BioPhylo")
travdf['TreeSize'].append(ts)
travdf['FindTime'].append(end-start)

travdf = pd.DataFrame(travdf)
travdf.to_csv("find_times.csv",index=False)
