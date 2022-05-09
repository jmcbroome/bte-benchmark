import subprocess
import pandas as pd
from functools import reduce
import sys
#each benchmark script is being ran as a separate subprocess instead of as a unified process
#because of how memory is used and freed by the subprocesses. This ensures that memory used in generating
#subtrees is always completely freed.

subprocess.check_call("python3 create_random_subsets.py public-2022-04-08.all.masked.nextclade.pangolin.pb.gz 1 6 10",shell=True)
sizes = [int(s.strip()) for s in open("subsets/sizes_saved.txt").readlines()]
alldf = []
try:
    for p in range(0):
        for s in sizes:
            print("Benchmarking subtree size {}...".format(s),file=sys.stderr)
            subprocess.check_call("python3 bench_loading.py subsets/{}_random.pb".format(s),shell=True)
            subprocess.check_call("python3 bench_traversal.py subsets/{}_random.pb".format(s),shell=True)
            #subprocess.check_call("python3 bench_subtree.py subsets/{}_random.pb".format(s),shell=True)
            subprocess.check_call("python3 bench_rsearch.py subsets/{}_random.pb".format(s),shell=True)
            #removing the find benchmark for now because everything is really fast and its less informative.
            #subprocess.check_call("python3 bench_find.py subsets/{}_random.pb".format(s),shell=True)
            #load the resulting dataframes and concatenate them
            df = pd.read_csv("traversal_times.csv")
            # df2 = pd.read_csv("subtree_times.csv")
            df3 = pd.read_csv("load_times.csv")
            df4 = pd.read_csv("rsearch_times.csv")
            df = reduce(lambda left,right: pd.merge(left,right,on=['Package','TreeSize'], how='outer'), [df,df3,df4])
            #df['Perm'] = p
            alldf.append(df)
except KeyboardInterrupt:
    pass
# df_merged = reduce(lambda left,right: pd.merge(left,right,on=['Package','TreeSize'], how='outer'), alldf)
df_merged = pd.concat(alldf,axis=0)
df_merged.to_csv("benchmarking_results_profiled_2.csv",index=False)