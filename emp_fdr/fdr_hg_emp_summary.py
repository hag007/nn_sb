import json
import matplotlib
from matplotlib_venn import venn2
matplotlib.use('Agg')

from pandas._libs.parsers import k
import sys
sys.path.insert(0, '../')

import argparse
import seaborn as sns
sns.set(color_codes=True)
import logging
sh = logging.StreamHandler()
logger = logging.getLogger("log")
logger.addHandler(sh)
from constants import *
from infra import *
from utils.param_builder import build_gdc_params
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import networkx as nx
import shutil
from pandas.errors import EmptyDataError
from utils.permute_network import EdgeSwapGraph
import scipy
from scipy.optimize import least_squares
from runners.FDR_runner import run_FDR
import random
from statsmodels.sandbox.stats.multicomp import fdrcorrection0
import multiprocessing
from utils.daemon_multiprocessing import func_star
def calc_emp_pval(cur_rv, cur_dist):
    pos = np.size(cur_dist) - np.searchsorted(np.sort(cur_dist), cur_rv, side='left')


    return pos / float(np.size(cur_dist))



def main(algo_sample = None, dataset_sample = None,tsv_file_name=os.path.join(constants.OUTPUT_GLOBAL_DIR, "emp_fdr", "MAX", "emp_diff_{}_{}_md.tsv")):
    output_md = pd.read_csv(
        tsv_file_name.format(dataset_sample, algo_sample),
        sep='\t', index_col=0).dropna()
    output_md = output_md.rename(columns={"filtered_pval": "hg_pval"})
    filtered_genes=output_md.loc[np.logical_and.reduce([output_md["n_genes"].values > 5, output_md["n_genes"].values < 500]), ["GO name","hg_pval", "emp_pval", "passed_oob_permutation_test"]]


    print "total n_genes with pval:{}/{}".format(np.size(filtered_genes["hg_pval"].values), 7435)

    sorted_genes_hg = filtered_genes.sort_values(by=['hg_pval'], ascending=False)
    sig_genes_hg_pval=np.append(sorted_genes_hg["hg_pval"].values,np.zeros(7435-np.size(sorted_genes_hg["hg_pval"].values)))
    sig_genes_hg_pval = [10**(-x) for x in sig_genes_hg_pval]
    fdr_results = fdrcorrection0(sig_genes_hg_pval, alpha=0.05, method='indep', is_sorted=False)
    n_hg_true = len([cur for cur in fdr_results[0] if cur])
    sig_hg_genes= sorted_genes_hg.iloc[:n_hg_true, :] if n_hg_true > 0 else 0
    HG_CUTOFF = 10**(-sig_hg_genes.iloc[- 1]["hg_pval"])
    print "HG cutoff: {}, n={}".format(HG_CUTOFF, len(sig_hg_genes.index))

    sorted_genes_emp = filtered_genes.sort_values(by=['emp_pval'])
    sorted_genes_emp.loc[sorted_genes_emp['emp_pval']==0,'emp_pval']=1.0/1000
    sig_genes_emp_pval = sorted_genes_emp["emp_pval"].values
    fdr_results = fdrcorrection0(sig_genes_emp_pval, alpha=0.05, method='indep', is_sorted=False)
    n_emp_true = len([cur for cur in fdr_results[0] if cur])
    sig_emp_genes = sorted_genes_emp.iloc[:n_emp_true, :]
    EMP_CUTOFF = sig_emp_genes.iloc[- 1]["emp_pval"] if n_emp_true > 0 else 0
    print "EMP cutoff: {}, n={}".format(EMP_CUTOFF, len(sig_emp_genes.index))

    # genes_oob = filtered_genes.loc[filtered_genes["passed_oob_permutation_test"],:]
    # print "OOB genes: n={}".format(len(genes_oob.index))



    return sig_hg_genes, sig_emp_genes


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='args')
    parser.add_argument('--datasets', dest='datasets', default="TNFa_2,HC12,SHERA,SHEZH_1,ROR_1,ERS_1,IEM")
    parser.add_argument('--prefix', dest='prefix', default="GE")
    parser.add_argument('--algos', dest='algos', default="jactivemodules_greedy,jactivemodules_sa,hotnet2,bionet,netbox,keypathwayminer_INES_GREEDY")

    args = parser.parse_args()

    datasets=args.datasets.split(",")
    algos=args.algos.split(",")
    prefix = args.prefix

    terms_limit = 0
    sig_terms_summary=pd.DataFrame()
    full_report=pd.DataFrame()
    df_all=pd.DataFrame()
    df_matrix=pd.DataFrame()
    for cur_ds in datasets:
        df_ds=pd.DataFrame()
        for cur_alg in algos:

            sig_hg_genes, sig_emp_genes=main(cur_alg, cur_ds,tsv_file_name=os.path.join("/home/hag007/Desktop/aggregate_report/oob", "emp_diff_{}_{}_passed_oob.tsv"))
            venn2([set(sig_hg_genes.index), set(sig_emp_genes.index)], set_labels = ('HG', 'EMP'))
            plt.title("EMP/HG ratio: {}".format(round(float(len(sig_emp_genes.index))/len(sig_hg_genes.index),3)))
            plt.savefig("/home/hag007/Desktop/venn/venn_{}_{}.png".format(cur_ds, cur_alg))
            plt.clf()
            df_ds.loc["{}_{}".format(cur_ds,cur_alg), "algo"]=cur_alg
            df_ds.loc["{}_{}".format(cur_ds,cur_alg), "dataset"]=cur_ds
            df_ds.loc["{}_{}".format(cur_ds,cur_alg),"n_emp"]=len(sig_emp_genes.index)
            df_ds.loc["{}_{}".format(cur_ds,cur_alg), "n_hg"] = len(sig_hg_genes.index)
            df_ds.loc["{}_{}".format(cur_ds,cur_alg), "ratio"] = round(float(len(sig_emp_genes.index))/len(sig_hg_genes.index),3)
            df_matrix.loc[cur_alg, cur_ds]=df_ds.loc["{}_{}".format(cur_ds,cur_alg), "ratio"]
        df_ds["ratio_rank"]=df_ds["ratio"].rank(ascending=False)
        for cur_alg in algos:
            df_matrix.loc[cur_alg, cur_ds] = df_ds.loc["{}_{}".format(cur_ds,cur_alg), "ratio_rank"]
        df_all=pd.concat([df_all,df_ds])

    df_all.to_csv("/home/hag007/Desktop/venn/summary.tsv", sep='\t')
    df_matrix.to_csv("/home/hag007/Desktop/venn/rank_matrix.tsv", sep='\t')