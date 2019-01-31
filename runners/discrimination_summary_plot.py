import wget
from utils.ensembl2gene_symbol import e2g_convertor
import time
import requests
import scipy.special
import matplotlib.pyplot as plt
from matplotlib import style
import matplotlib.ticker as ticker
style.use("ggplot")
import logging
sh = logging.StreamHandler()
logger = logging.getLogger("log")
logger.addHandler(sh)
from infra import *
import os
import pandas as pd
from sklearn.decomposition import PCA
from mpl_toolkits.mplot3d import proj3d
from matplotlib.lines import Line2D
import random
import matplotlib.cm as cm
import matplotlib.colors as ml_colors
from utils.omic_svm import apply_svm, svm_linear_default, DISTANCE ,make_meshgrid, plot_contours


def main():
    high_th=0.9
    algos = ["hotnet2", "jactivemodules_greedy",
             "jactivemodules_sa", "keypathwayminer_INES_GREEDY", "netbox", "bionet"] # "matisse", "reactomefi"

    datapoints = []
    discrimination_files = [os.path.join(constants.OUTPUT_GLOBAL_DIR, "discrimination", name) for name in os.listdir(os.path.join(constants.OUTPUT_GLOBAL_DIR, "discrimination")) if
                  not os.path.isdir(os.path.join(constants.OUTPUT_GLOBAL_DIR, "discrimination", name))]

    fig = plt.figure(1, figsize=(15, 15))
    ax = fig.add_subplot(111)
    colormap = cm.rainbow
    df_summary = pd.DataFrame()

    for discrimination_file in discrimination_files:
        is_countable = []
        df_discrimination = pd.read_csv(discrimination_file, sep='\t', index_col=0)
        df_discrimination = df_discrimination.loc[~df_discrimination.index.isin(["reactomefi", "matisse"]), :]
        for algo, row in df_discrimination.iterrows():
            datapoints.append((row['algo_top_sig_ratio'], row['algo_kmean_ratio'], ml_colors.rgb2hex(colormap(algos.index(algo)/float(len(algos)-1)))))
            is_countable.append(datapoints[-1][0]>high_th and datapoints[-1][1]>high_th)

        df_discrimination["is_countable"] = pd.Series(is_countable, index=df_discrimination.index)
        df_summary = pd.concat((df_summary, df_discrimination))

        for x, y, c in datapoints:
            ax.scatter(x, y, c=c)
        colorlist = [ml_colors.rgb2hex(colormap(i)) for i in
                     np.array(list(range(len(algos)))) / float(len(algos) - 1)]
    patches = [Line2D([0], [0], marker='o', color='gray', label=a,
                      markerfacecolor=c) for i, a, c in zip(list(range(len(algos))), algos, colorlist)]
    ax.legend(handles=patches)
    ax.grid(True)
    ax.set_xlabel("top_sig_ratio")
    ax.set_ylabel("kmean_ratio")

    plt.savefig(
        os.path.join(constants.OUTPUT_GLOBAL_DIR, "discrimination.png".format(constants.DATASET_NAME)))

    plt.clf()
    plt.cla()
    datapoints = []
    df_mean = df_summary.groupby(df_summary.index)['algo_top_sig_ratio', 'algo_kmean_ratio', 'is_countable'].agg(['mean', 'std', 'sum', 'count'])
    algos=list(df_mean.index.values)
    for algo, row in df_mean.iterrows():
        datapoints.append((algo, row[("algo_top_sig_ratio", "mean")], row[("algo_top_sig_ratio", "std")], row[('algo_kmean_ratio','mean')], row[('algo_kmean_ratio','std')], row[('algo_kmean_ratio','count')], row[('is_countable','sum')],
                           ml_colors.rgb2hex(colormap(algos.index(algo) / float(len(algos)-1)))))



    plt.clf()
    plt.cla()

    fig = plt.figure(2, figsize=(10, 10))
    ax = fig.add_subplot(111)
    details=[]
    for algo, x, x_std, y, y_std, count, num_of_countable, c in datapoints:
        ax.scatter(x, y, count*20, c=c)
        ax.errorbar(x, y, xerr=x_std,  yerr=y_std, linestyle='None', marker='^',c=c)
        ax.annotate("{}, {}".format(int(count), int(num_of_countable)), (x,y))
        details.append("{}: got results in {}/{} datasets. of them {} were high".format(algo, int(count), len(discrimination_files), int(num_of_countable)))
    colorlist = [ml_colors.rgb2hex(colormap(i)) for i in
                 np.array(list(range(len(algos)))) / float(len(algos) - 1)]

    props = dict(boxstyle='round', facecolor='gray', alpha=0.5)
    fig.text(0.01, 0.85,
             "\n".join(details),
             bbox=props)
    patches = [Line2D([0], [0], marker='o', color='gray', label=a,
                      markerfacecolor=c) for i, a, c in zip(list(range(len(algos))), algos, colorlist)]
    ax.legend(handles=patches)
    ax.grid(True)
    ax.set_xlabel("top_sig_ratio")
    ax.set_ylabel("kmean_ratio")


    plt.savefig(
        os.path.join(constants.OUTPUT_GLOBAL_DIR, "discrimination_avg.png".format(constants.DATASET_NAME)))

main()
