
import sys
sys.path.insert(0, '../')

import pandas as pd
from pandas.errors import EmptyDataError

from fastsemsim.SemSim import *
from fastsemsim.Ontology import ontologies
from fastsemsim.Ontology import AnnotationCorpus

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from rpy2.robjects import pandas2ri
pandas2ri.activate()

import constants

import matplotlib.cm as cm
import matplotlib.colors as ml_colors

import argparse


def main(datasets, algos):

    colormap = cm.rainbow
    colorlist = [ml_colors.rgb2hex(colormap(i)) for i in
                 np.array(list(range(len(algos)))) / float(len(algos) - 1)]
    df_matrix = pd.DataFrame()
    df_summary = pd.DataFrame()
    for cur_ds in datasets:

        constants.update_dirs(DATASET_NAME_u=cur_ds)
        total_num_genes=[]
        avg_num_genes=[]
        std_num_genes=[]
        algos_signals=[]
        algo_go_sims = []

        for i_algo, cur_algo in enumerate(algos):
            print "current aggregation: {}, {}".format(cur_ds,cur_algo)
            try:
                total_num_genes.append(pd.read_csv(
                    os.path.join(constants.OUTPUT_GLOBAL_DIR, constants.DATASET_NAME, cur_algo, "all_modules_general.tsv"),
                    sep="\t")["total_num_genes"][0])
                avg_num_genes.append(pd.read_csv(
                    os.path.join(constants.OUTPUT_GLOBAL_DIR, constants.DATASET_NAME, cur_algo,
                                 "modules_summary.tsv"),
                    sep="\t")["#_genes"].mean())
                std_num_genes.append(pd.read_csv(
                    os.path.join(constants.OUTPUT_GLOBAL_DIR, constants.DATASET_NAME, cur_algo,
                                 "modules_summary.tsv"),
                    sep="\t")["#_genes"].std())
            except:
                print "no genes were found for: {}, {}".format(cur_ds, cur_algo)
                total_num_genes.append(0)
            algos_signals.append(float(file(os.path.join(constants.OUTPUT_GLOBAL_DIR, "emp_fdr", "ds_2_alg_scores",
                              "{}_{}_{}".format(cur_ds, cur_algo, "n_sig.txt"))).read()))
            algo_go_sims.append(float(file(os.path.join(constants.OUTPUT_GLOBAL_DIR, "emp_fdr", "ds_2_alg_scores",
                                                       "{}_{}_{}".format(cur_ds, cur_algo, "var.txt"))).read()))

        fig, ax = plt.subplots(figsize=(10, 10))

        print "all data: \n{}\n{}\n{}\n{}".format(algos_signals, algo_go_sims, algos, total_num_genes)
        for h, s, c, a, gene_size, module_mean, module_std in zip(algos_signals, algo_go_sims, colorlist, algos,
                                         total_num_genes, avg_num_genes, std_num_genes):  # [0 for x in range(len(algo_go_sim_score))]
            print (h, s)
            ax.scatter(h, s, s=(50 + 2000 * (float(gene_size) / (1+np.max(total_num_genes)))),
                       c=c, cmap='jet', label=a)
            df_series=pd.Series({"algo": a, "dataset": cur_ds, "sig_terms": h,
                       "sig_terms_rank": pd.Series(np.array(algos_signals)).rank(ascending=0).values[
                           np.where(np.array(algos_signals) == h)[0][0]], "variability": s,
                       "variability_rank": pd.Series(np.array(algo_go_sims)).rank(ascending=0).values[
                           np.where((np.array(algo_go_sims)) == s)[0][0]], "n_genes": gene_size, "module_size_mean": module_mean, "module_size_std": module_std})
            df_series.name = "{}_{}".format(cur_ds, a)
            df_summary=df_summary.append(df_series)
            df_matrix.loc[a, cur_ds]=h
            colorlist = [ml_colors.rgb2hex(colormap(i)) for i in
                         np.array(list(range(len(algos)))) / float(len(algos) - 1)]
            patches = [Line2D([0], [0], marker='o', color='gray', label=a,
                              markerfacecolor=c) for i, a, c in zip(list(range(len(algos))), algos, colorlist)]
            ax.set_xlabel("# GO terms (-log10(qval)) above threshold")
            ax.set_ylabel("GO terms variability")
            ax.legend(handles=patches)
            ax.grid(True)
        plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR,
                                 "hs_plot_terms_signal_algo_{}.png".format(constants.DATASET_NAME)))
    return df_summary, df_matrix


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='args')
    parser.add_argument('--datasets', dest='datasets', default="SOC")
    parser.add_argument('--prefix', dest='prefix', default="GE")
    parser.add_argument('--algos', dest='algos', default="jactivemodules_greedy")

    args = parser.parse_args()


    prefix = args.prefix
    datasets=["{}_{}".format(prefix,x) for x in args.datasets.split(",")]
    algos = args.algos.split(",")

    go_ratio_ds_summary = pd.DataFrame()
    ds_summary=pd.DataFrame()
    ds_summary, df_matrix=main(datasets=datasets, algos=algos)
    ds_summary.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "ds_go_rank_summary.tsv"), sep='\t')
    df_matrix.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "ds_go_matrix.tsv"), sep='\t')



