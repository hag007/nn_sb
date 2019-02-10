import scipy.spatial.distance as ssd
import numpy as np
import pandas as pd
import constants
import os
import math
from scipy.cluster.hierarchy import dendrogram, linkage
from matplotlib import pyplot as plt


distance_matrix=pd.read_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "cancer_type_go_distance.tsv"),sep='\t', index_col=0).dropna(axis=0)


files=distance_matrix.index.values # np.append(distance_matrix.index.values,["coad.csv"])

for cur_x in files:
    for cur_y in files:
        if cur_y==cur_x:
            distance_matrix.loc[cur_y,cur_x]=0
            continue
        try:
            if math.isnan(distance_matrix.loc[cur_x,cur_y]):
                distance_matrix.loc[cur_x, cur_y] =distance_matrix.loc[cur_y, cur_x]
        except:
            distance_matrix.loc[cur_x, cur_y] = distance_matrix.loc[cur_y, cur_x]



# convert the redundant n*n square matrix form into a condensed nC2 array
distArray = ssd.squareform(distance_matrix[distance_matrix.index.values].values)

linked=linkage(distArray, method='single', metric='euclidean')

plt.figure(figsize=(10, 7))
dendrogram(linked,
            orientation='top',
            labels=[x.split(".")[0] for x in distance_matrix.index.values],
            distance_sort='descending',
            show_leaf_counts=True)
plt.show()