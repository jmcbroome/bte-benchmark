import subprocess
import pandas as pd
from functools import reduce
#each benchmark script is being ran as a separate subprocess instead of as a unified process
#because of how memory is used and freed by the subprocesses. This ensures that memory used in generating
#subtrees is always completely freed.

subprocess.check_call("python3 create_random_subsets.py public-2022-04-08.all.masked.nextclade.pangolin.pb.gz 1 5 10",shell=True)
sizes = [int(s.strip()) for s in open("subsets/sizes_saved.txt").readlines()]
alldf = []
try:
    for s in sizes:
        print("Benchmarking subtree size {}...".format(s))
        subprocess.check_call("python3 bench_loading.py subsets/{}_random.pb".format(s),shell=True)
        subprocess.check_call("python3 bench_traversal.py subsets/{}_random.pb".format(s),shell=True)
        subprocess.check_call("python3 bench_subtree.py subsets/{}_random.pb".format(s),shell=True)
        subprocess.check_call("python3 bench_rsearch.py subsets/{}_random.pb".format(s),shell=True)
        subprocess.check_call("python3 bench_find.py subsets/{}_random.pb".format(s),shell=True)
        #load the resulting dataframes and concatenate them
        df = pd.read_csv("traversal_times.csv")
        df2 = pd.read_csv("subtree_times.csv")
        df3 = pd.read_csv("load_times.csv")
        df4 = pd.read_csv("rsearch_times.csv")
        df5 = pd.read_csv("find_times.csv")
        df = reduce(lambda left,right: pd.merge(left,right,on=['Package','TreeSize'], how='outer'), [df,df2,df3,df4,df5])
        alldf.append(df)
except KeyboardInterrupt:
    pass
# df_merged = reduce(lambda left,right: pd.merge(left,right,on=['Package','TreeSize'], how='outer'), alldf)
df_merged = pd.concat(alldf,axis=0)
df_merged.to_csv("benchmarking_results.csv",index=False)