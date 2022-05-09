from Bio import Phylo
import ete3
import bte
import pandas as pd
import sys
import time
from profile_decor import named_profile

ts = int(sys.argv[1].split("/")[1].split("_")[0])

@named_profile("BTE\tSubtree\t" + str(ts))
def btesubtree(tree,targets):
    tree.subtree(targets)
@named_profile("ETE\tSubtree\t" + str(ts))
def etesubtree(tree,targets):
    #ete3 prune actually prunes everything BUT the input nodes.
    tree.prune(targets)
@named_profile("BioPhylo\tSubtree\t" + str(ts))
def physubtree(tree,targets):
    #phylo does not have a single-line prune, so instead we have to manually prune all non-targets.
    for t in tree.get_terminals():
        if t.name not in targets:
            phyt.prune(t)


subtdf = {k:[] for k in ['Package','TreeSize','SubtreeSize','SubtreeTime']}
matpb = bte.MATree(sys.argv[1])
#remove mutations from the tree by writing a newick and reloading from the newick. This removes systematic bias against BTE, as it 
#maintains mutations in the tree structure and the other packages do not, meaning BTE has higher memory usage (because more data is stored)
#instead, if we load from a newick alone, no mutations will be included.
nwk = matpb.write_newick()
mat = bte.MATree(nwk_string = nwk)
target_size = round(ts/10)
targets = mat.get_random(target_size).get_leaves_ids()
start = time.perf_counter()
btesubtree(mat,targets)
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
etesubtree(etet,targets)
end = time.perf_counter()
subtdf['Package'].append("ETE")
subtdf['TreeSize'].append(ts)
subtdf['SubtreeSize'].append(target_size)
subtdf['SubtreeTime'].append(end-start)

phyt = Phylo.read("subtree.nwk", "newick")
start = time.perf_counter()
physubtree(phyt,targets)
end = time.perf_counter()
subtdf['Package'].append("BioPhylo")
subtdf['TreeSize'].append(ts)
subtdf['SubtreeSize'].append(target_size)
subtdf['SubtreeTime'].append(end-start)

subtdf = pd.DataFrame(subtdf)
subtdf.to_csv("subtree_times.csv",index=False)