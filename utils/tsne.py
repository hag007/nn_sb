import numpy as np
from sklearn.manifold import TSNE
import constants
import time
import os
import matplotlib.pyplot as plt

def plot_tsne(dataset, labels_assignment, meta_groups, tested_gene_list_file_name=None, n_components=3):

    actual_labels = list(range(1, len(meta_groups[0]) + 1))

    X = []
    y = []
    for cur in actual_labels:
        for cur_sample in dataset[np.where(cur == labels_assignment[0]), :][0]:
            X.append(cur_sample)
            y.append(cur)
    X = np.array(X)
    y = np.array(y)

    X = TSNE(n_components=n_components).fit_transform(X)

    fig = plt.figure(1, figsize=(20, 20))
    plt.clf()
    if n_components == 3:
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(X[:, 0], X[:, 1], X[:, 2], c=y, cmap='jet')
    if n_components == 2:
        ax = fig.add_subplot(111)
        ax.scatter(X[:, 0], X[:, 1], c=y, cmap='jet')
    plt.savefig(
        os.path.join(constants.BASE_PROFILE, "output", "TSNE_by_samples_{}_{}_{}_{}.png").format(constants.CANCER_TYPE,
                                                                                                tested_gene_list_file_name.split(
                                                                                                    ".")[
                                                                                                    0] if tested_gene_list_file_name is not None else "none",
                                                                                                n_components,
                                                                                                time.time()))