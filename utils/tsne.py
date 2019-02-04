import numpy as np
from sklearn.manifold import TSNE
import constants
import time
import os
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.colors as ml_colors
import matplotlib.cm as cm

def plot_tsne(dataset, labels_assignment, meta_groups, tested_gene_list_file_name=None, n_components=3):

    groups_md = [(cur["_name"], int(cur["_label"])) for i, cur in enumerate(meta_groups[0])]
    groups_md = [("unknown", 0)]+groups_md

    label_ids_unique = np.unique(labels_assignment[0])
    labels = [x[0] for x in groups_md if x[1] in label_ids_unique]
    label_ids = [x[1] for x in groups_md if x[1] in label_ids_unique]

    X = []
    y = []
    for cur in label_ids_unique :
        for cur_sample in dataset[np.where(cur == labels_assignment[0]), :][0]:
            X.append(cur_sample)
            y.append(cur)
    X = np.array(X)
    y = np.array(y)

    X = TSNE(n_components=n_components, metric="correlation", perplexity=30.0).fit_transform(X)

    fig = plt.figure(1, figsize=(20, 20))
    plt.clf()
    if n_components == 3:
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(X[:, 0], X[:, 1], X[:, 2], c=y, cmap='jet')
    if n_components == 2:
        ax = fig.add_subplot(111)
        ax.scatter(X[:, 0], X[:, 1], c=y, cmap='jet')

    colormap=cm.jet
    colorlist_unique = [ml_colors.rgb2hex(colormap(i)) for i in
                        label_ids_unique / float(max(label_ids))]
    patches = [Line2D([0], [0], marker='o', color='gray', label=a,
                      markerfacecolor=c) for a, c in
               zip(labels, colorlist_unique)]
    ax.legend(handles=patches)

    plt.savefig(
        os.path.join(constants.BASE_PROFILE, "output", "TSNE_by_samples_{}_{}_{}_{}.png").format(constants.CANCER_TYPE,
                                                                                                tested_gene_list_file_name.split(
                                                                                                    ".")[
                                                                                                    0] if tested_gene_list_file_name is not None else "none",
                                                                                                n_components,
                                                                                                time.time()))
