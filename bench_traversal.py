from Bio import Phylo
import ete3
import bte
import pandas as pd
import sys
import time
from profile_decor import named_profile

travdf = {k:[] for k in ['Package','TreeSize','TraversalTime']}
matpb = bte.MATree(sys.argv[1])
#remove mutations from the tree by writing a newick and reloading from the newick. This removes systematic bias against BTE, as it 
#maintains mutations in the tree structure and the other packages do not, meaning BTE has higher memory usage (because more data is stored)
#instead, if we load from a newick alone, no mutations will be included.
nwk = matpb.write_newick()
mat = bte.MATree(nwk_string = nwk)
ts = int(sys.argv[1].split("/")[1].split("_")[0])

@named_profile("BTE\tTraversal\t" + str(ts))
def btetraverse(tree):
    tree.depth_first_expansion()
@named_profile("ETE\tTraversal\t" + str(ts))
def etetraverse(tree):
    for n in tree.traverse('preorder'):
        pass
@named_profile("BioPhylo\tTraversal\t" + str(ts))
def phytraverse(tree):
    for n in tree.find_clades():
        pass

start = time.perf_counter()
btetraverse(mat)
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
etetraverse(etet)
end = time.perf_counter()
travdf['Package'].append("ETE")
travdf['TreeSize'].append(ts)
travdf['TraversalTime'].append(end-start)

phyt = Phylo.read("subtree.nwk", "newick")
start = time.perf_counter()
phytraverse(phyt)
end = time.perf_counter()
travdf['Package'].append("BioPhylo")
travdf['TreeSize'].append(ts)
travdf['TraversalTime'].append(end-start)

travdf = pd.DataFrame(travdf)
travdf.to_csv("traversal_times.csv",index=False)
