#!/usr/bin/env python
"""Calculate differentially expressed genes using EdgeR from bioconductor.
http://bioconductor.org/packages/2.5/bioc/html/edgeR.html
Usage:
    count_diffexp.py <count_file>
"""
import sys
sys.path.insert(0, '../')

import random
import os
import numpy as np
import pandas as pd
from pandas.errors import EmptyDataError
import subprocess
import networkx as nx
from networkx.drawing.nx_agraph import write_dot, graphviz_layout

import fastsemsim
from fastsemsim.SemSim import *
from fastsemsim.Ontology import ontologies
from fastsemsim.Ontology import AnnotationCorpus

import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D

# import rpy2.robjects.numpy2ri  as numpy2ri
# numpy2ri.activate()

from rpy2.robjects import pandas2ri
pandas2ri.activate()

import constants

from utils.r_runner import run_rscript
from utils.network import get_network_genes
from utils.network import build_all_reports
import utils.server as server
import infra

import utils.go
import utils.go_hierarcies
import DEG_runner

from utils.ensembl2entrez import ensembl2entrez_convertor
from utils.scripts import format_script
from utils.server import get_parameters
from utils.server import get_score
from utils.network import output_modules

from sklearn.decomposition import PCA
from mpl_toolkits.mplot3d import proj3d

import matplotlib.cm as cm
import matplotlib.colors as ml_colors
from matplotlib.colors import ListedColormap
import matplotlib.pylab as pl

ENABLE_GO_GRAPH = False
IS_GO_GRAPH_ONLY=False
GO_PCA=False
TERMS_SIMILARITY_TO_NUM_OF_TERMS=True
RATIO_TO_GO_TERM = True
QVAL_TH = 0.001



def rewire_edges(cur_node, active_father, vertices, edges, contained_nodes, path):
    if cur_node in contained_nodes and cur_node != active_father:
        # cur_edge = (active_father, cur_node)
        # if cur_node not in edges:
        for cur in path:
            edges.append(cur)

        active_father=cur_node
        path = []

    for cur in vertices[cur_node]["obj"].children:
        child_path = list(path)
        child_path.append((cur_node, cur.id))
        rewire_edges(cur.id, active_father, vertices, edges, contained_nodes, child_path)


def plot_go_tree(dict_result, positive_terms, algo_name, module_index):
    print "num of positive terms: {}".format(len(positive_terms))

    for cur_root, cur_tree in dict_result.iteritems():
        fig, axs = plt.subplots(1, 1, figsize=(100, 100))
        edges = []
        rewire_edges(cur_root, cur_root, cur_tree["vertices"], edges, positive_terms, [])
        G = nx.DiGraph(edges)
        pos = graphviz_layout(G, prog='dot')
        nc = nx.draw_networkx_nodes(G, pos,
                                    node_size=[5000 if x in positive_terms else 0 for x in G.nodes], with_labels=True)
        nx.draw_networkx_edges(G, pos, edge_color="black", alpha=0.4)
        nx.draw_networkx_labels(G, pos, {x: "".join(["\n" for a in range(int(random.random()*10))]) + "\n" + cur_tree["vertices"][x]["name"] for x in G.nodes if
                                         x in positive_terms or x == cur_root}, font_size=24, font_color='blue')
        props = dict(boxstyle='round', facecolor='gray', alpha=0.5)
        fig.text(0.01, 0.01,
                 "\n".join(["{}: {}".format(x, df_go_metadata.loc[x, "GO name"]) for x in G.nodes if x in positive_terms]),
                 bbox=props, fontsize=36)

        print "number of nodes in tree {}: {}".format(cur_root, len(G.nodes))
        plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "go_tree_{}_{}_{}.png".format(constants.DATASET_NAME, algo_name, cur_root)))
        pl.close(fig)
        plt.clf()


# def init_graph():
#     dict_result, go2geneids, geneids2go, entrez2ensembl = utils.go_hierarcies.build_hierarcy()
#     G = nx.DiGraph()
#     G.add_node(ROOT)
#     G = nx.DiGraph([tuple(x.split("=")) for x in dict_result.values()[0]["edges"].keys()])
#     G.add_nodes_from([(k, v) for k, v in dict_result.values()[0]["vertices"].iteritems()])
#     cmap = pl.cm.bwr
#     my_cmap = cmap(np.arange(cmap.N))
#     my_cmap[:, -1] = np.linspace(0, 1, cmap.N)
#     my_cmap = ListedColormap(my_cmap)
#     return dict_result, G, my_cmap


def plot_pca(all_hg_score_labels, df_all_hg_pval, ml_colors, algos):
    n_components = 2
    colormap = cm.rainbow
    pca = PCA(n_components=n_components)
    y = np.array(all_hg_score_labels)
    X = np.transpose(df_all_hg_pval.values)
    X = pca.fit_transform(X)
    fig = plt.figure(1, figsize=(15, 15))
    plt.clf()
    # colorlist = [ml_colors.rgb2hex(colormap(i)) for i in all_hg_score_labels]
    colorlist = [ml_colors.rgb2hex(colormap(i)) for i in
                 np.array(all_hg_score_labels) / float(np.max(all_hg_score_labels))]
    if n_components == 3:
        ax = fig.add_subplot(111, projection='3d')
        for x, c, a in zip([x for x in X], colorlist, algos):
            ax.scatter(x[0], x[1], x[2], c=c, cmap='jet', label=a)
    elif n_components == 2:
        ax = fig.add_subplot(111)
        for x, c in zip([x for x in X], colorlist):
            ax.scatter(x[0], x[1], c=c, cmap='jet')
    colorlist = [ml_colors.rgb2hex(colormap(i)) for i in
                 np.array(list(range(len(algos)))) / float(len(algos) - 1)]
    patches = [Line2D([0], [0], marker='o', color='gray', label=a,
                      markerfacecolor=c) for i, a, c in zip(list(range(len(algos))), algos, colorlist)]
    ax.legend(handles=patches)
    ax.grid(True)
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    coeff = np.transpose(pca.components_[0:2, :]) * 150
    labels = df_all_hg_pval.index.values
    top_pca_x_go_terms = labels[np.abs(coeff[:, 0]).argsort()[::-1][:5]]
    top_pca_x_details = ["{}: {}".format(x, df_go_metadata.loc[x, "GO name"]) for x in top_pca_x_go_terms]
    top_pca_y_go_terms = labels[np.abs(coeff[:, 1]).argsort()[::-1][:5]]
    top_pca_y_details = ["{}: {}".format(x, df_go_metadata.loc[x, "GO name"]) for x in top_pca_y_go_terms]
    top_pca_go_terms = np.unique(np.append(top_pca_x_go_terms, top_pca_y_go_terms))
    props = dict(boxstyle='round', facecolor='gray', alpha=0.5)
    fig.text(0.01, 0.01,
             "\n".join(["PC1 top GO terms:"] + top_pca_x_details + ["PC2 top GO terms:"] + top_pca_y_details),
             bbox=props)
    n = coeff.shape[0]
    for i in range(n):
        if labels[i] in top_pca_go_terms:
            ax.arrow(0, 0, coeff[i, 0], coeff[i, 1], color='r', alpha=0.5)
            ax.annotate(labels[i], tuple(coeff[i]), (-100 + random.random() * 200, -100 + random.random() * 200),
                        textcoords='offset points',
                        bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.5),
                        arrowprops=dict(facecolor='black', shrink=0.05, width=2, headwidth=3, headlength=2))
    for module_i, txt in enumerate(all_hg_score_modules):
        ax.annotate(str(txt), (X[:, 0][module_i], X[:, 1][module_i]))
    plt.savefig(
        os.path.join(constants.OUTPUT_GLOBAL_DIR, "PCA_by_samples_{}.png".format(constants.DATASET_NAME)))


def main(dataset_name):
    global df_go_metadata, all_hg_score_modules, df_hg_output
    colormap = cm.rainbow
    constants.update_dirs(DATASET_NAME_u=dataset_name)
    GO_RANK_CUTOFF = 150
    ##########################
    if TERMS_SIMILARITY_TO_NUM_OF_TERMS:
        ontology_type = 'GeneOntology'
        ignore_parameters = {'ignore': {}}
        source_type = 'obo'
        source = os.path.join(os.path.join(constants.GO_DIR, constants.GO_FILE_NAME))

        print "\n######################"
        print "# Loading ontology... #"
        print "######################\n"

        ontology = ontologies.load(source=source, source_type=source_type, ontology_type=ontology_type,
                                   parameters=ignore_parameters)

        print "\n######################"
        print "# Loading Annotation Corpus... #"
        print "######################\n"
        ac = AnnotationCorpus.AnnotationCorpus(ontology)
        ac.parse(os.path.join(constants.GO_DIR, "goa_human.gaf"), "gaf-2.0")
        ac.isConsistent()

        print "\n#################################"
        print "# Annotation corpus successfully loaded."
        print "#################################\n"

        semsim = GSESAMESemSim(ontology, ac)  # maxSemSim(ontology, ac) #
    #################
    if ENABLE_GO_GRAPH:
        dict_result, go2geneids, geneids2go, entrez2ensembl = utils.go_hierarcies.build_hierarcy(
            roots=['GO:0008150', 'GO:0005575', 'GO:0003674'])
    #################
    all_homogeneity = []
    all_separability = []
    agg_homogeneity = []
    agg_separability = []
    algo_go_sim_score = []
    colors = []
    df_all_hg_pval = pd.DataFrame()
    df_go_metadata = pd.DataFrame()
    all_hg_score_labels = []
    all_hg_score_modules = []
    labels_by_sample = []
    total_num_genes = []
    algos = ["keypathwayminer_INES_GREEDY", "netbox", "hotnet2", "jactivemodules_greedy", "bionet",
             "jactivemodules_sa"]  # "matisse", "reactomefi"  # "keypathwayminer_INES_OPTIMAL", "keypathwayminer_INES_ACO"
    algos_signals = []
    modules_signals = []
    df_all_hg_qval = pd.DataFrame()
    df_module2best_rows = []
    df_module2avg_rows = []
    df_algo2best_rows = []
    df_algo2avg_rows = []
    for i_algo, ALGO_NAME in enumerate(algos):
        df_all_hg_qval = pd.DataFrame()
        print "current algo: {}".format(ALGO_NAME)

        go2modules = {}
        modules2go = {}
        homogeneity = []
        separability = []
        all_go_terms = []

        try:
            total_num_genes.append(pd.read_csv(
                os.path.join(constants.OUTPUT_GLOBAL_DIR, constants.DATASET_NAME, ALGO_NAME, "all_modules_general.tsv"),
                sep="\t")["total_num_genes"][0])
        except pd.errors.EmptyDataError:
            total_num_genes.append(0)
        algo2best_go_ratio = 0
        algo2avg_go_ratio = 0
        module2best_go_ratio=[]
        module2avg_go_ratio = []
        df_modules_summary = pd.read_csv(
            os.path.join(constants.OUTPUT_GLOBAL_DIR, constants.DATASET_NAME, ALGO_NAME, "modules_summary.tsv"),
            sep='\t')
        i = -1
        for i in range(len(df_modules_summary.index)):
            hg_file_name = os.path.join(constants.OUTPUT_GLOBAL_DIR, constants.DATASET_NAME, ALGO_NAME,
                                        "module_{}_separated_modules_hg_samples.tsv".format(i))
            print "reading module: {} from file".format(i)
            if os.path.getsize(hg_file_name) < 2:
                modules2go[i] = np.array([])
                modules_signals.append(0)

            else:
                df_hg_output = pd.read_csv(hg_file_name, sep="\t")
                df_hg_output.index = df_hg_output["GO id"]
                df_go_metadata = pd.concat([df_go_metadata, df_hg_output[["GO id", "GO name"]]], axis=0)
                df_go_metadata = df_go_metadata[~df_go_metadata.index.duplicated(keep='first')]
                df_all_hg_pval = pd.concat([df_all_hg_pval, df_hg_output["pval"].apply(lambda x: -np.log10(x))],
                                           join='outer', axis=1)
                df_all_hg_qval = pd.concat([df_all_hg_qval, df_hg_output["qval"].apply(lambda x: -np.log10(x))],
                                           join='outer', axis=1)
                df_hg_output = df_hg_output.iloc[:min(len(df_hg_output.index), GO_RANK_CUTOFF), :]
                df_hg_output = df_hg_output[
                    df_hg_output["qval"] <= QVAL_TH]  # .iloc[:min(len(df_hg_output.index),5),:]  #

                score = df_hg_output['value'].iloc[0] / float(df_modules_summary.loc[i, "#_genes"]) if len(
                    df_hg_output.index) > 0 else 0
                df_module2best_rows.append({'name' : "{}_{}".format(ALGO_NAME ,i), 'algo': ALGO_NAME, 'module' : i, 'score' : score, 'num_of_genes' : df_modules_summary.loc[i, "#_genes"]})
                algo2best_go_ratio += score
                score = (df_hg_output['value'] / float(df_modules_summary.loc[i, "#_genes"]) if len(
                    df_hg_output.index) > 0 else np.array([0])).mean()
                df_module2avg_rows.append({'name' : "{}_{}".format(ALGO_NAME, i), 'algo': ALGO_NAME, 'module' : i, 'score' : score, 'num_of_genes' : df_modules_summary.loc[i, "#_genes"]})
                algo2avg_go_ratio += score
                modules_signals.append(len(df_hg_output.index))
                all_hg_score_modules.append(i)
                all_hg_score_labels.append(i_algo)
                for x in df_hg_output["GO id"]:
                    if x in go2modules:
                        go2modules[x].append(i)
                    else:
                        go2modules[x] = [i]

                modules2go[str(i)] = df_hg_output["GO id"].values
                all_go_terms = np.append(all_go_terms, df_hg_output["GO id"].values)

        i+=1
        if RATIO_TO_GO_TERM:
            df_algo2best_rows.append(
                {'name': '{}_total_avg'.format(ALGO_NAME), 'score': algo2best_go_ratio / max(i, 1)})
            df_algo2avg_rows.append({'name': '{}_total_avg'.format(ALGO_NAME), 'score': algo2avg_go_ratio / max(i, 1)})


        df_all_hg_pval[pd.isna(df_all_hg_pval)] = 0
        df_all_hg_qval[pd.isna(df_all_hg_qval)] = 0
        max_per_go_term = 0
        if df_all_hg_qval.values.size > 0:
            max_per_go_term = np.sum(np.max(df_all_hg_qval.values, axis=1) >= -np.log10(QVAL_TH))
        algos_signals.append(max_per_go_term)
        all_go_terms = list(np.unique(all_go_terms))

        print "added signal : {}".format(algos_signals[-1])
        print "all_go_terms : {}".format(len(all_go_terms))

        if ENABLE_GO_GRAPH:
            plot_go_tree(dict_result, all_go_terms, ALGO_NAME, i)

        if ENABLE_GO_GRAPH and IS_GO_GRAPH_ONLY:
            continue

        adj = np.ones((len(all_go_terms), len(all_go_terms))) * (-2)

        if TERMS_SIMILARITY_TO_NUM_OF_TERMS:
            for i_x, x in enumerate(all_go_terms):
                print "calc distance between terms {}/ {}".format(i_x, len(all_go_terms))
                for i_y, y in enumerate(all_go_terms):
                    if adj[i_x, i_y] != -2: continue
                    adj[i_x, i_y] = semsim.SemSim(x, y)  # , ResnikSemSim(ontology,ac))
                    if np.isnan(adj[i_x, i_y]):
                        adj[i_x, i_y] = -1
                    adj[i_y, i_x] = adj[i_x, i_y]

            algo_go_sim_score.append(
                [1 if np.isnan(x) else x for x in [np.sum(adj[adj != -1]) / (np.size(adj) - np.sum(adj == -1))]][0])

            for k, v in sorted([(int(k), v) for k, v in modules2go.iteritems()], key=lambda x: x[0]):
                print "calc homogeneity and seperability for module: {}".format(k)
                v_filtered = [x for x in v]  # if x in G.nodes
                labels_by_sample.append(ALGO_NAME)
                homogeneity.append(np.nan_to_num(np.sum(
                    [adj[all_go_terms.index(x), all_go_terms.index(y)] for x in v for y in v if
                     adj[all_go_terms.index(x), all_go_terms.index(y)] != -1 and x != y]) / (
                                                             len(v_filtered) * (len(v_filtered) - 1))))
                separability.append(np.nan_to_num(np.sum(
                    [adj[all_go_terms.index(x), all_go_terms.index(y)] for x in v for y in all_go_terms if
                     y not in v and adj[all_go_terms.index(x), all_go_terms.index(y)] != -1]) / (
                                                          len(v_filtered) * (len(all_go_terms) - len(v_filtered)))))

            all_separability = all_separability + separability
            all_homogeneity = all_homogeneity + homogeneity

            agg_separability.append(np.average(separability))
            agg_homogeneity.append(np.average(homogeneity))

            fig, ax = plt.subplots(figsize=(15, 15))

            ax.scatter(homogeneity, separability)
            ax.legend()
            ax.set_xlabel("Intra-Similarity")
            ax.set_ylabel("Inter-Similarity")
            ax.grid(True)
            colors = colors + [float(i_algo) / len(algos) for x in range(len(modules2go))]

            for module_i, txt in enumerate(range(i)):
                ax.annotate(str(txt), (homogeneity[module_i], separability[module_i]))
            plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR,
                                     "hs_plot_{}_{}.png".format(ALGO_NAME, constants.DATASET_NAME)))
    if RATIO_TO_GO_TERM:
        pd.DataFrame(df_module2best_rows).set_index('name').to_csv(
            os.path.join(constants.OUTPUT_GLOBAL_DIR, constants.DATASET_NAME, "GO_terms_per_module_ratio_best.tsv"),
            sep='\t')
        pd.DataFrame(df_module2avg_rows).set_index('name').to_csv(
            os.path.join(constants.OUTPUT_GLOBAL_DIR, constants.DATASET_NAME, "GO_terms_per_module_ratio_avg.tsv"),
            sep='\t')

        pd.DataFrame(df_algo2best_rows).set_index('name').to_csv(
            os.path.join(constants.OUTPUT_GLOBAL_DIR, constants.DATASET_NAME, "GO_terms_ratio_best.tsv"), sep='\t')
        pd.DataFrame(df_algo2avg_rows).set_index('name').to_csv(
            os.path.join(constants.OUTPUT_GLOBAL_DIR, constants.DATASET_NAME, "GO_terms_ratio_avg.tsv"), sep='\t')
    if TERMS_SIMILARITY_TO_NUM_OF_TERMS:
        fig, ax = plt.subplots(figsize=(15, 15))
        colorlist = [ml_colors.rgb2hex(colormap(i)) for i in colors]
        for h, s, c, a in zip(all_homogeneity, all_separability, colorlist, labels_by_sample):
            ax.scatter(h, s, s=50, c=c, vmin=0, vmax=1, cmap='jet')
            colorlist = [ml_colors.rgb2hex(colormap(i)) for i in
                         np.array(list(range(len(algos)))) / float(len(algos) - 1)]
            patches = [Line2D([0], [0], marker='o', color='gray', label=a,
                              markerfacecolor=c) for i, a, c in zip(list(range(len(algos))), algos, colorlist)]
            ax.legend(handles=patches)
            ax.set_xlabel("Intra-Similarity")
            ax.set_ylabel("Inter-Similarity")
            ax.grid(True)
        plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "hs_plot_all_{}.png".format(constants.DATASET_NAME)))

        fig, ax = plt.subplots(figsize=(10, 10))
        colorlist = [ml_colors.rgb2hex(colormap(i)) for i in colors]
        for h, s, c, a in zip(modules_signals, all_separability, colorlist, labels_by_sample):
            ax.scatter(h, s, s=50, c=c, vmin=0, vmax=1, cmap='jet')
        colorlist = [ml_colors.rgb2hex(colormap(i)) for i in
                     np.array(list(range(len(algos)))) / float(len(algos) - 1)]
        patches = [Line2D([0], [0], marker='o', color='gray', label=a,
                          markerfacecolor=c) for i, a, c in zip(list(range(len(algos))), algos, colorlist)]
        ax.legend(handles=patches)
        ax.set_xlabel("# GO terms (-log10(qval)) above threshold")
        ax.set_ylabel("Algorithm Inter-Similarity")
        ax.grid(True)
        plt.savefig(
            os.path.join(constants.OUTPUT_GLOBAL_DIR, "hs_plot_signal_all_{}.png".format(constants.DATASET_NAME)))

        colorlist = [ml_colors.rgb2hex(colormap(i)) for i in
                     np.array(list(range(len(agg_homogeneity)))) / float(len(agg_homogeneity) - 1)]
        fig, ax = plt.subplots(figsize=(10, 10))
        for h, s, c, a in zip(agg_homogeneity, agg_separability, colorlist, algos):
            ax.scatter(h, s, s=50, c=c, cmap='jet', label=a)
        ax.set_xlabel("Intra-Similarity")
        ax.set_ylabel("Inter-Similarity")
        ax.legend()
        ax.grid(True)

        for module_i, txt in enumerate(algos):
            ax.annotate(str(txt), (agg_homogeneity[module_i], agg_separability[module_i]))
        plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "hs_plot_agg_{}.png".format(constants.DATASET_NAME)))
    # colorlist = [ml_colors.rgb2hex(colormap(i)) for i in algo_go_sim_score/np.max(algo_go_sim_score)]

        colorlist = [ml_colors.rgb2hex(colormap(i)) for i in
                     np.array(list(range(len(algos)))) / float(len(algos) - 1)]

        fig, ax = plt.subplots(figsize=(10, 10))
        for h, s, c, a, gene_size in zip(algos_signals, algo_go_sim_score, colorlist, algos,
                                         total_num_genes):  # [0 for x in range(len(algo_go_sim_score))]
            print (h, s)
            ax.scatter(h, s, s=(50 + 2000 * (float(gene_size) / np.max(total_num_genes))),
                       c=c, cmap='jet', label=a)

            colorlist = [ml_colors.rgb2hex(colormap(i)) for i in
                         np.array(list(range(len(algos)))) / float(len(algos) - 1)]
            patches = [Line2D([0], [0], marker='o', color='gray', label=a,
                              markerfacecolor=c) for i, a, c in zip(list(range(len(algos))), algos, colorlist)]
            ax.set_xlabel("# GO terms (-log10(qval)) above threshold")
            ax.set_ylabel("Algorithm Inter-Similarity")
            ax.legend(handles=patches)
            ax.grid(True)
        plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR,
                                 "hs_plot_terms_signal_algo_{}.png".format(constants.DATASET_NAME)))
    if GO_PCA:
        plot_pca(all_hg_score_labels, df_all_hg_pval, ml_colors, algos)


if __name__ == "__main__":

    go_ratio_ds_summary = pd.DataFrame()
    # datasets = [name for name in os.listdir(constants.DATASETS_DIR) if
    #             os.path.isdir(os.path.join(constants.DATASETS_DIR, name)) and name.startswith("GWAS_")]
    datasets=["TNFa_2", "MCF7_2", "SOC", "HC12", "IEM", "IES"]

    for cur_ds in datasets:
        print "current dataset: {}".format(cur_ds)
        main(dataset_name=cur_ds)
        go_ratio_ds_summary=pd.concat([go_ratio_ds_summary, pd.read_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, cur_ds,"GO_terms_ratio_avg.tsv"), sep='\t', index_col=0).rename(columns={"score": cur_ds})], axis=1  )
    go_ratio_ds_summary['avg'] = go_ratio_ds_summary.mean(axis=1)
    temp = go_ratio_ds_summary['avg'].values.argsort()
    ranks = np.empty_like(temp)
    ranks[temp] = np.arange(len(temp))[ : :-1]
    go_ratio_ds_summary['rank'] = ranks
    go_ratio_ds_summary.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR,"ds_go_score_summary.tsv"),sep='\t')
    for cur_col in go_ratio_ds_summary.columns:
        temp = go_ratio_ds_summary[cur_col].values.argsort()
        ranks = np.empty_like(temp)
        ranks[temp] = np.arange(len(temp))[::-1]
        go_ratio_ds_summary[cur_col]=ranks
    go_ratio_ds_summary['avg'] = go_ratio_ds_summary.mean(axis=1)
    temp = go_ratio_ds_summary['avg'].values.argsort()
    ranks = np.empty_like(temp)
    ranks[temp] = np.arange(len(temp))
    go_ratio_ds_summary['rank'] = ranks
    go_ratio_ds_summary.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "ds_go_rank_summary.tsv"), sep='\t')


