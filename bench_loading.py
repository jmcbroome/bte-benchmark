from Bio import Phylo
import ete3
import bte
import pandas as pd
import sys
import time
from profile_decor import named_profile

ts = int(sys.argv[1].split("/")[1].split("_")[0])

@named_profile("BTE\tload\t" + str(ts))
def bteload(tname):
    mat = bte.MATree(nwk_string = tname)
    return mat
@named_profile("ETE\tload\t" + str(ts))
def eteload(nwk):
    return ete3.Tree(nwk,format=1)
@named_profile("BioPhylo\tload\t" + str(ts))
def phyload():
    return Phylo.read("subtree.nwk", "newick")

#first, load the tree with bte and save the time taken.
loaddf = {k:[] for k in ['Package','TreeSize','LoadTime']}
matpb = bte.MATree(sys.argv[1])
#remove mutations from the tree by writing a newick and reloading from the newick. This removes systematic bias against BTE, as it 
#maintains mutations in the tree structure and the other packages do not, meaning BTE has higher memory usage (because more data is stored)
#instead, if we load from a newick alone, no mutations will be included.
nwk = matpb.write_newick()
start = time.perf_counter()
mat = bteload(nwk)
end = time.perf_counter()
loaddf['Package'].append("BTE")
loaddf['TreeSize'].append(ts)
loaddf['LoadTime'].append(end-start)
#ete3 and biopython.phylo need newicks, so write that.
nwk = mat.write_newick()
nwk_f = open("subtree.nwk","w+")
print(nwk,file=nwk_f)
nwk_f.close()
#next, load the tree with ete3
start = time.perf_counter()
eteload(nwk)
end = time.perf_counter()
loaddf['Package'].append("ETE")
loaddf['TreeSize'].append(ts)
loaddf['LoadTime'].append(end-start)
#finally, load it with biopython.phylo
start = time.perf_counter()
phyload()
end = time.perf_counter()
loaddf['Package'].append("BioPhylo")
loaddf['TreeSize'].append(ts)
loaddf['LoadTime'].append(end-start)
#save the result
loaddf = pd.DataFrame(loaddf)
loaddf.to_csv("load_times.csv",index=False)