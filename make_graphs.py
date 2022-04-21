import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt

bdf = pd.read_csv("benchmarking_results.csv")
for c in ['SubtreeTime','LoadTime','PathTime','TraversalTime']:
    sns.lineplot(x=np.log10(bdf['TreeSize']),y=np.log10(bdf[c]),hue='Package',data=bdf)
    plt.ylabel("Log"+c)
    plt.xlabel("LogTreeSize")
    plt.savefig(c+"_benchmark.png")