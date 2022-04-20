from Bio import Phylo
import ete3
import bte
import pandas as pd
import sys
import time

#first, load the tree with bte and save the time taken.
loaddf = {k:[] for k in ['Package','TreeSize','LoadTime']}
start = time.perf_counter()
mat = bte.MATree(sys.argv[1])
ts = int(sys.argv[1].split("/")[1].split("_")[0])
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
ete3.Tree(nwk,format=1)
end = time.perf_counter()
loaddf['Package'].append("ETE")
loaddf['TreeSize'].append(ts)
loaddf['LoadTime'].append(end-start)
#finally, load it with biopython.phylo
start = time.perf_counter()
Phylo.read("subtree.nwk", "newick")
stop = time.perf_counter()
loaddf['Package'].append("BioPhylo")
loaddf['TreeSize'].append(ts)
loaddf['LoadTime'].append(end-start)
#save the result
loaddf = pd.DataFrame(loaddf)
loaddf.to_csv("load_times.csv",index=False)