from sklearn.cluster import KMeans
import numpy as np
from scipy.stats import rankdata
import matplotlib.pyplot as plt
import os
from utils.kmeans import kmeanssample
import constants
import time
from utils.ensembl2gene_symbol import e2g_convertor


def find_clusters(end_k, gene_expression_top_var, gene_expression_top_var_headers_rows, start_k,
                gene_expression_top_var_headers_columns, tested_gene_list_file_name, labels_assignment=None, phenotype_heatmap=None, clustering_algorithm="euclidean", plot=True):
    clfs_results = {}
    for n_clusters in range(start_k, end_k + 1):
        clfs_results[n_clusters] = []

        centres, km_clf, dist = kmeanssample(k=n_clusters, X=gene_expression_top_var, metric=clustering_algorithm)
        ranks = []
        for i in range(n_clusters):
            ranks.append(np.average(np.delete(gene_expression_top_var, np.where(km_clf != i)[0], axis=0)))
        ranks = rankdata(ranks)
        cluster_labels = np.array(km_clf)
        for i in range(n_clusters):
            cluster_labels[np.where(km_clf == ranks[i] - 1)] = i
        if labels_assignment is not None:
            labels_assignment = [cluster_labels + 1] + labels_assignment
        else:
            labels_assignment = [cluster_labels + 1]

        for i in range(n_clusters):
            cluster_indices = np.where(cluster_labels != i)[0]
            gene_expression_cluster = np.delete(gene_expression_top_var, cluster_indices, axis=0)
            gene_headers_rows_cluster = np.delete(gene_expression_top_var_headers_rows, cluster_indices, axis=0)
            clfs_results[n_clusters].append(gene_headers_rows_cluster)
            desc = "k={} clustering cluster {} has {} patients".format(n_clusters, i, len(gene_headers_rows_cluster))

        if plot:
            plot_heatmap(gene_expression_top_var, gene_expression_top_var_headers_columns,
                         labels_assignment, gene_expression_top_var_headers_rows,
                         tested_gene_list_file_name, n_clusters, phenotype_heatmap=phenotype_heatmap)
    return clfs_results

def plot_heatmap(gene_expression_top_var, gene_expression_top_var_headers_columns, labels_assignment,
                 gene_expression_top_var_headers_rows, tested_gene_list_file_name, n_clusters = None, label_index=None ,phenotype_heatmap=None):
    cluster_labels = labels_assignment[0]
    fig = plt.figure()
    abs_left = 0.05
    abs_bottom = 0.05
    abs_height = 0.9
    main_width = 0.8
    label_width = 0.025
    axes = []

    if len(gene_expression_top_var_headers_columns) < 80:
        abs_bottom = 0.1
    if len(gene_expression_top_var_headers_rows) < 80:
        abs_left = 0.2
        main_width = 0.6

    ax1 = fig.add_axes([abs_left, abs_bottom, main_width, abs_height])
    axes.append(ax1)

    if len(gene_expression_top_var_headers_columns) < 80:
        ax1.set_xticks(np.arange(len(gene_expression_top_var_headers_columns)))
        ax1.set_xticklabels(gene_expression_top_var_headers_columns)
    if len(gene_expression_top_var_headers_rows) < 80:
        ax1.set_yticks(np.arange(len(gene_expression_top_var_headers_rows)))
        ax1.set_yticklabels(gene_expression_top_var_headers_rows[cluster_labels.argsort()])

    for i, cur in enumerate(labels_assignment + [phenotype_heatmap]):
        if cur is None: continue
        ax2 = fig.add_axes(
            [abs_left + main_width + label_width + (i * 2) * label_width, abs_bottom, label_width, abs_height])
        axes.append(ax2)
        ax2.grid(False)
        ax2.set_xticks([])
        ax2.set_yticks([])
        data = ax2.imshow(cur[cluster_labels.argsort()].reshape((len(cluster_labels), 1)), cmap='jet', aspect='auto')
    # ax = plt.subplot(111)
    for label in ax1.xaxis.get_ticklabels():
        label.set_fontsize(7)
        label.set_rotation(90)
    for label in ax1.yaxis.get_ticklabels():
        label.set_fontsize(7)
    data = ax1.imshow(gene_expression_top_var[cluster_labels.argsort(), :], cmap='jet', aspect='auto')
    cb = plt.colorbar(data, ax=axes, fraction=0.05, pad=0.04)
    interval = 0 # abs(np.min(gene_expression_top_var) - np.max(gene_expression_top_var))*0.3
    data.set_clim(np.percentile(gene_expression_top_var,5) , np.percentile(gene_expression_top_var,90) ) # - interval
    plt.savefig(os.path.join(constants.BASE_PROFILE, "output",
                             "heatmap_cluster_by_p_{}_{}_k={}_label_i={}_{}.png".format(constants.CANCER_TYPE,
                                                                                        tested_gene_list_file_name.split("/")[-1].split('.')[0],
                                                                                        n_clusters, label_index, time.time())))
    # "heatmap_cluster_by_p_{}_{}_k={}.svg".format(constants.CANCER_TYPE, tested_gene_list_file_name.split(".")[0], n_clusters)), format='svg', dpi=1200)
    # plt.show()
    cb.remove()
    plt.cla()