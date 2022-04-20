import sys
import bte
import numpy as np
print("Loading initial tree...")
mat = bte.MATree(sys.argv[1])
min, max, count = int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
sizes = np.logspace(start=min,stop=max,num=count)
sizes = [round(s) for s in sizes]
print("Using the following sizes:",sizes)
print("Creating subsets...")
for s in sizes:
    mat.get_random(s).save_pb("subsets/{}_random.pb".format(s))
with open("subsets/sizes_saved.txt","w+") as f:
    for s in sizes:
        print(s,file=f)