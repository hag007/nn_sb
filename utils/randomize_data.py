import json
from matplotlib import style
from pandas._libs.parsers import k

style.use("ggplot")
import seaborn as sns
sns.set(color_codes=True)
import logging
sh = logging.StreamHandler()
logger = logging.getLogger("log")
logger.addHandler(sh)
from constants import *
from infra import *
from param_builder import build_gdc_params
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import networkx as nx
import shutil
from pandas.errors import EmptyDataError
from utils.permute_network import EdgeSwapGraph
from functools import partial

import random


def get_permutation_name(prefix, dataset, algo, index):
    random_ds_name=prefix + "_random_" + dataset
    if algo is not None:
        random_ds_name += "_{}".format(algo)
    if index is not None:
        random_ds_name += "_{}".format(index)

    return random_ds_name


def permutation_output_exists(prefix, dataset, algo, index):
    # print "try to find results for {}... {}".format(
    #         os.path.join(constants.OUTPUT_GLOBAL_DIR, get_permutation_name(prefix, dataset, algo, index), algo,
    #                      "modules_summary.tsv"),
    #     os.path.exists(
    #         os.path.join(constants.OUTPUT_GLOBAL_DIR, get_permutation_name(prefix, dataset, algo, index), algo,
    #                      "modules_summary.tsv"))
    # )
    return os.path.exists(
        os.path.join(constants.OUTPUT_GLOBAL_DIR, get_permutation_name(prefix, dataset, algo, index), algo,
                     "modules_summary.tsv"))

def create_random_ds(prefix, cur_ds, index=None, algo=None):
    data_type = "score.tsv"
    if prefix=="GE":
        data_type="ge.tsv"
    cur_ds = cur_ds[len(prefix)+1:]
    random_ds_name=get_permutation_name(prefix, cur_ds, algo, index)
    root_random_dir=os.path.join(constants.DATASETS_DIR, random_ds_name)
    if os.path.exists(root_random_dir):
        shutil.rmtree(root_random_dir)
    os.makedirs(os.path.join(root_random_dir, "data"))
    os.makedirs(os.path.join(root_random_dir, "output"))
    os.makedirs(os.path.join(root_random_dir, "cache"))
    data_file_name = os.path.join(root_random_dir, "data", data_type)
    classes_file_name = os.path.join(root_random_dir, "data", "classes.tsv")
    if prefix=="GE":
        shutil.copy(os.path.join(constants.DATASETS_DIR, "{}_{}".format(prefix, cur_ds), 'data', "classes.tsv"), classes_file_name)
    data=pd.read_csv(os.path.join(constants.DATASETS_DIR, "{}_{}".format(prefix, cur_ds), 'data', data_type), sep='\t', index_col=0)
    np.random.seed(int(random.random()*10000))
    data=pd.DataFrame(data=np.random.permutation(data.values), index=data.index, columns=data.columns)
    data.to_csv(data_file_name, sep='\t', index_label="id")
    return random_ds_name

def create_permuted_network(network_file_name):
    network_file_data=file(os.path.join(constants.NETWORKS_DIR, network_file_name+".sif")).readlines()
    edges=[(x.split('\t')[0].strip(), x.split('\t')[2].strip()) for x in network_file_data[1:]]

    G=EdgeSwapGraph()
    G.add_edges_from(edges)
    G=G.randomize_by_edge_swaps(10)
    permuted_network_file_data=[[network_file_data[0].split('\t')[0] , network_file_data[0].split('\t')[2]]] + [[x[0], x[1]] for x in G.edges()]
    permuted_network_file_name=network_file_name+"_perm.sif"
    file(os.path.join(constants.NETWORKS_DIR, permuted_network_file_name),'w+').write("\n".join(["\t".join([x[0], "ppi", x[1]]) for x in permuted_network_file_data]))
    return os.path.splitext(permuted_network_file_name)[0]




if __name__ == "__main__":
    prefix="GE"
    # algos = np.array(["netbox", "hotnet2", "keypathwayminer_INES_GREEDY", "jactivemodules_sa", "jactivemodules_greedy", "bionet"])
    algos = np.array(["jactivemodules_greedy"])
    for cur_algo in algos:
        algos_filter = cur_algo

        # datasets = [name for name in os.listdir(constants.OUTPUT_GLOBAL_DIR) if
        #             os.path.isdir(os.path.join(constants.OUTPUT_GLOBAL_DIR, name)) and name.startswith("GWAS_") and not name.startswith("GWAS_random")]
        # datasets = ["TNFa_2", "MCF7_2", "SOC", "HC12", "IEM", "IES"]
        datasets = ["TNFa_2"]
        df_go = pd.DataFrame(columns=['qval', 'pval'])
        for cur_ds in datasets:
            go_results = [os.path.join(constants.OUTPUT_GLOBAL_DIR, cur_ds, cur_algo, cur_module) for cur_algo in os.listdir(os.path.join(constants.OUTPUT_GLOBAL_DIR, cur_ds))
                          if os.path.isdir(os.path.join(constants.OUTPUT_GLOBAL_DIR, cur_ds, cur_algo)) and cur_algo in algos_filter for cur_module in os.listdir(os.path.join(constants.OUTPUT_GLOBAL_DIR, cur_ds, cur_algo)) if "separated_modules" in cur_module]

            for cur in go_results:
                try:
                    df_go = pd.concat((df_go, pd.read_csv(cur,sep='\t')))
                except EmptyDataError:
                    pass

        fig, ax = plt.subplots(figsize=(12, 10))

        df_go=df_go[df_go['qval']<0.05]
        df_go = df_go[~df_go.index.duplicated(keep='first')]
        pval=-np.log10(df_go["pval"].values)
        if np.size(pval)==0:
            pval=np.array([1])
        sns.distplot(pval, kde=False)
        plt.title("pval dist. -log10(pval) go terms. algo: {}".format("".join(algos_filter)))
        real_n_term=np.size(pval)
        txt = "# terms={}".format(np.size(pval))
        fig.text(.5, .05, txt, ha='center')
        plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "go_pval_dist_{}_{}.png".format(prefix.lower(), "".join(algos_filter))))
        plt.clf()
        # datasets = [name for name in os.listdir(constants.OUTPUT_GLOBAL_DIR) if
        #             os.path.isdir(os.path.join(constants.OUTPUT_GLOBAL_DIR, name)) and name.startswith("GWAS_random")]
        datasets = ["GE_random_TNFa_2"] # , "GE_random_MCF7_2", "GE_random_SOC", "GE_random_HC12", "GE_random_IEM", "GE_random_IES"]
        df_go = pd.DataFrame(columns=['qval', 'pval'])
        for cur_ds in datasets:
            go_results = [os.path.join(constants.OUTPUT_GLOBAL_DIR, cur_ds, cur_algo, cur_module) for cur_algo in os.listdir(os.path.join(constants.OUTPUT_GLOBAL_DIR, cur_ds))
                          if os.path.isdir(os.path.join(constants.OUTPUT_GLOBAL_DIR, cur_ds, cur_algo)) and cur_algo in algos_filter for cur_module in os.listdir(os.path.join(constants.OUTPUT_GLOBAL_DIR, cur_ds, cur_algo)) if "separated_modules" in cur_module]
            for cur in go_results:
                try:
                    df_go = pd.concat((df_go, pd.read_csv(cur,sep='\t')))
                except EmptyDataError:
                    pass

        fig, ax = plt.subplots(figsize=(12, 10))
        df_go = df_go[df_go['qval'] < 0.05]
        df_go = df_go[~df_go.index.duplicated(keep='first')]
        pval=-np.log10(df_go["pval"].values)
        if np.size(pval)==0:
            pval=np.array([1])
        sns.distplot(pval, kde=False)
        plt.title("random pval dist. -log10(pval) go terms. algo: {}".format("".join(algos_filter)))
        random_n_term = np.size(pval)
        txt="# terms={}".format(np.size(pval))
        fig.text(.5, .05, txt, ha='center')
        plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "go_pval_dist_random_{}_{}.png".format(prefix.lower(), "".join(algos_filter))))
        print "{}\t{}\t{}\t{}".format(algos_filter, real_n_term, random_n_term, round(real_n_term/float(random_n_term),2))
        # for cur_ds in ["GWAS_random_2hr_glucose", "GWAS_random_adhd", "GWAS_random_alzheimers"]:
        #     pval = pd.read_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, cur_ds, "data", "score.tsv"), sep='\t')["pval"]
        #     sns.distplot(pval, kde=False)
        #     plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "pval_dist_gwas_pval.png"))

        # for cur_ds in ["TNFa_2", "MCF7_2", "SOC"]:
        #     pval = pd.read_csv(os.path.join(constants.DATASETS_DIR, cur_ds, "cache", "deg_edger.tsv"), sep='\t')["pval"]
        #     sns.distplot(pval, kde=False)
        #     plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "pval_dist_rnaseq.png"))
        #
        # for cur_ds in ["IEM", "IES"]:
        #     pval = pd.read_csv(os.path.join(constants.DATASETS_DIR, cur_ds, "cache", "deg_t.tsv"), sep='\t')["pval"]
        #     sns.distplot(pval, kde=False)
        #     plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "pval_dist_microarray.png"))




