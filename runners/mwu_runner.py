from random import shuffle
import numpy as np
num_of_modules = [4,6,100,30,70, 20, 40, 10]
from scipy.stats import mannwhitneyu

all_i = []
for i, cur in enumerate(num_of_modules):
    all_i += [i for a in range(cur)]

pvals = []
for cur in range(1000):
    print cur
    shuffle(all_i)
    print all_i
    for a in range(len(num_of_modules)):
        x = np.where(np.array(all_i) == a)[0]
        y = np.where(np.array(all_i) != a)[0]
        pvals.append(mannwhitneyu(x,y).pvalue)

pvals = np.array(pvals)
pvals.sort()
print pvals[pvals.shape[0]/20]






