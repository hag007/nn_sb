import scipy.spatial.distance as ssd
import numpy as np
import pandas as pd
import constants
import os
import math
from scipy.cluster.hierarchy import dendrogram, linkage
from matplotlib import pyplot as plt


distance_matrix=pd.read_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "cancer_type_go_distance.tsv"),sep='\t', index_col=0)


distance_matrix[np.isnan(distance_matrix)]=1

distance_matrix+=(np.abs(np.min(distance_matrix.values))+1)
for cur in distance_matrix.index:
    distance_matrix.loc[cur,cur]=0


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

linked=linkage(distArray, method='weighted', metric='euclidean')

plt.figure(figsize=(10, 7))
dendrogram(linked,
            orientation='top',
            labels=[x.split(".")[0] for x in distance_matrix.index.values],
            distance_sort='descending',
            show_leaf_counts=True)

# label_colors = {'jactivemodules_greedy': 'red', 'jactivemodules_sa': 'green', 'hotnet2': 'blue', 'bionet': 'black', 'netbox': 'brown', 'keypathwayminer_INES_GREEDY': 'orange'}



label_colors = {'TNFa_2': 'red', 'HC12': 'green', 'SHERA': 'blue', 'ROR_1': 'black', 'SHEZH_1': 'brown', 'ERS_1': 'orange', 'IEM': 'purple'}

ax = plt.gca()
xlbls = ax.get_xmajorticklabels()
for lbl in xlbls:
    for cur_label_color in label_colors:
        if cur_label_color in lbl.get_text():
            lbl.set_color(label_colors[cur_label_color])


plt.tight_layout()
plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "dataset_hierarchical.png"))
plt.show()