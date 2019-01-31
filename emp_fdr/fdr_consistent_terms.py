import json
import matplotlib
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



def main(algo_sample = None, dataset_sample = None, n_dist_samples = 300, n_total_samples = None, shared_list=None, n_start_i=None, limit = 10000):
    output_md = pd.read_csv(
        os.path.join(constants.OUTPUT_GLOBAL_DIR, "emp_fdr", "MAX", "emp_diff_{}_{}_md.tsv".format(dataset_sample, algo_sample)),
        sep='\t', index_col=0).dropna()
    # constants.OUTPUT_GLOBAL_DIR, "emp_fdr", "MAX", "emp_diff_{}_{}_md.tsv"
    output_md = output_md.rename(columns={"filtered_pval": "hg_pval"})
    n_genes_pvals=output_md.loc[np.logical_and.reduce([output_md["n_genes"].values > 5, output_md["n_genes"].values < 500]), "hg_pval"].values

    print "total n_genes with pval:{}/{}".format(np.size(n_genes_pvals), 7435)
    n_genes_pvals=np.append(n_genes_pvals,np.zeros(7435-np.size(n_genes_pvals)))
    n_genes_pvals = [10**(-x) for x in n_genes_pvals]
    fdr_results = fdrcorrection0(n_genes_pvals, alpha=0.05, method='indep', is_sorted=False)
    n_hg_true = len([cur for cur in fdr_results[0] if cur == True])
    HG_CUTOFF=np.sort(n_genes_pvals)[n_hg_true-1]
    print "HG cutoff: {}".format(HG_CUTOFF)

    output_md = output_md.loc[np.logical_and.reduce(
        [output_md["n_genes"].values > 5, output_md["n_genes"].values < 500,
         output_md["hg_pval"].values > HG_CUTOFF]), :]

    output = pd.read_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "emp_fdr", "MAX", "emp_diff_{}_{}.tsv").format(dataset_sample, algo_sample),
        sep='\t', index_col=0).dropna()
    output = output.rename(columns={"filtered_pval": "hg_pval"})
    output = output.loc[output_md.index.values, :]
    counter = 0
    emp_dists = []
    emp_pvals = []

    n_total_samples=n_total_samples if n_total_samples is not None else len(output.iloc[0].loc["dist_n_samples"][1:-1].split(", "))
    if n_start_i is None:
        np.random.seed(int(random.random()*1000))        
        i_choice=np.random.choice(n_total_samples, n_dist_samples, replace=False)
        i_dist=i_choice[:n_dist_samples]
    else:
        i_dist=np.arange(n_start_i, n_start_i+n_dist_samples)



    for index, cur in output.iterrows():
        if counter == limit: break
        # print "cur iteration in real data: {}/{}".format(counter, len(output.index))
        pval = np.array([float(x) for x in cur["dist_n_samples"][1:-1].split(", ")])[i_dist]

        emp_pvals.append(calc_emp_pval(cur["hg_pval"], pval))
        emp_dists.append(pval)
        counter += 1
    dist_name = "emp"
    df_dists = pd.DataFrame(index=output.index)
    df_dists["emp"] = pd.Series(emp_dists, index=output.index[:limit])

    emp_pvals = [x if x != 0 else 1.0 / n_dist_samples for x in emp_pvals]
    zero_bool=[x<=0.004 for x in emp_pvals]
    fdr_bh_results = fdrcorrection0(emp_pvals, alpha=0.05, method='indep', is_sorted=False)[0]
    fdr_robust_results=[True if x <0.05 else False for x in run_FDR(emp_pvals)]
    mask_terms=fdr_robust_results
    go_ids_result=output.index.values[mask_terms]
    go_names_result=output["GO name"].values[mask_terms]
    n_emp_true =sum(mask_terms)
    BH_TH = np.sort(emp_pvals)[n_emp_true - 1]

    print "BH cutoff: {} # true terms passed BH cutoff: {}".format(BH_TH, n_emp_true)
    if shared_list is not None:
        shared_list.append((BH_TH, n_emp_true, HG_CUTOFF, n_hg_true, go_ids_result, go_names_result))
    return BH_TH, n_emp_true, HG_CUTOFF, n_hg_true, go_ids_result, go_names_result


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='args')
    parser.add_argument('--datasets', dest='datasets', default="SOC")
    parser.add_argument('--prefix', dest='prefix', default="GE")
    parser.add_argument('--algos', dest='algos', default="jactivemodules_greedy")
    parser.add_argument('--n_iteration', dest='n_iteration', default=100)
    parser.add_argument('--n_total_samples', help="n_total_samples", dest='n_total_samples', default=1000)
    parser.add_argument('--n_dist_samples', help="n_dist_samples", dest='n_dist_samples', default=200)
    parser.add_argument('--pf', dest='pf', help="parallelization factor", default=20)

    args = parser.parse_args()

    datasets=args.datasets.split(",")
    algos=args.algos.split(",")
    prefix = args.prefix
    n_iteration=int(args.n_iteration)
    n_total_samples=int(args.n_total_samples)
    n_dist_samples = int(args.n_dist_samples)
    pf=int(args.pf)

    terms_limit = 0
    sig_terms_summary=pd.DataFrame()
    full_report=pd.DataFrame()
    p = multiprocessing.Pool(pf)

    for cur_ds in datasets:
        for cur_alg in algos:

            l_emp_cutoff=[]
            l_n_emp_true=[]
            l_hg_cutoff=[]
            l_n_hg_true=[]

            go_ids_results=[]
            go_names_results=[]
            go_names_intersection=[]
            ys1=[]
            ys2=[]
            print "{}_{}".format(cur_alg, cur_alg)

            shared_list = multiprocessing.Manager().list()
            params=[]
            print "about to start {} intersections".format(n_iteration)
            for cur in range(n_iteration):
                params.append([main, [cur_alg, cur_ds, n_dist_samples, n_total_samples, shared_list]]) # , n_dist_samples*cur

            p.map(func_star, params)            
           
            for BH_TH, n_emp_true, HG_CUTOFF, n_hg_true, go_ids_result, go_names_result in list(shared_list):
                l_emp_cutoff.append(str(BH_TH))
                l_n_emp_true.append(str(n_emp_true))
                l_hg_cutoff.append(str(HG_CUTOFF))
                l_n_hg_true.append(str(n_hg_true))


                if len(go_names_result) >terms_limit:
                    go_ids_results.append(go_ids_result)
                    go_names_results.append(go_names_result)
                    go_names_intersection = list(set(go_names_results[0]).intersection(*go_names_results))

                    print "total intersection size : {}".format(len(go_names_intersection))
                    ys1.append(len(go_names_intersection))
                    ys2.append(len(go_names_result))

            plt.clf()

            go_ids_intersection = list(set(go_ids_results[0]).intersection(*go_ids_results)) if len(go_ids_results) > 0 else []

            output_md = pd.read_csv(
                os.path.join(constants.OUTPUT_GLOBAL_DIR, "emp_fdr", "MAX",
                             "emp_diff_{}_{}_md.tsv".format(cur_ds, cur_alg)),
                sep='\t', index_col=0)

            output_md.loc[go_ids_intersection,"passed_oob_permutation_test"]=True
            output_md.loc[~np.isin(output_md.index.values, np.array(go_ids_intersection)), "passed_oob_permutation_test"] = False

            output_md.to_csv(
                os.path.join(constants.OUTPUT_GLOBAL_DIR, "emp_fdr", "MAX",
                             "emp_diff_{}_{}_passed_oob.tsv".format(cur_ds, cur_alg)),
                sep='\t')

            sig_terms_summary.loc[cur_alg,cur_ds]=len(go_names_intersection)
            modules_summary=pd.read_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR,"GE_{}".format(cur_ds),cur_alg,"modules_summary.tsv"), sep='\t')
            full_report=full_report.append({"index" : "{}_{}".format(cur_alg, cur_ds), "dataset" : cur_ds, "algo" : cur_alg, "n_mutual_sig_terms" : len(go_names_intersection), "mean_sig_terms" : np.mean(ys2) if len(ys2) > 0 else 0, "n_genes" : modules_summary["#_genes"].sum() , "n_modules" : modules_summary["#_genes"].count(), "modules_size_mean" : modules_summary["#_genes"].mean(), "module_size_std" : modules_summary["#_genes"].sum()}, ignore_index=True)

            fig,ax = plt.subplots()
            plt.plot(np.arange(len(ys1)),ys1, label='# intersected terms')
            plt.plot(np.arange(len(ys2)), ys2, label='# sig. terms in permutation')
            plt.title("min # terms above threshold ({}):  ({}/{}). ratio: {}/{}".format(terms_limit, len(ys1), n_iteration, ys1[-1] if len(ys1) > 0 else 0,
                                                                                min([len(x) for x in go_names_results]) if len(go_names_results) > 0 else 0 ))
            ax.set_xlabel("iteration index")
            ax.set_ylabel("# terms")

            plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "intersection_function_{}_{}.png".format(cur_ds, cur_alg)))

            fdr_test_summary=pd.DataFrame()
            fdr_test_summary["emp_cutoff"]=l_emp_cutoff
            fdr_test_summary["n_emp_true"]=l_n_emp_true
            fdr_test_summary["hg_cutoff"]=l_hg_cutoff
            fdr_test_summary["n_hg_true"]=l_n_hg_true

            fdr_test_summary.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "emp_fdr", "fdr_test_summary_{}_{}.tsv".format(cur_ds, cur_alg)), sep='\t')

    sig_terms_summary.to_csv(
        os.path.join(constants.OUTPUT_GLOBAL_DIR, "emp_fdr", "sig_terms_summary.tsv"),
        sep='\t')

    rank_col=[]
    for cur_ds_i, cur_ds in enumerate(datasets):
        rank_col = np.append(rank_col,full_report["n_mutual_sig_terms"].iloc[cur_ds_i:cur_ds_i+len(algos)].rank(ascending=1))

    full_report["mutual_sig_terms_rank"]=rank_col
    full_report.to_csv(
        os.path.join(constants.OUTPUT_GLOBAL_DIR, "emp_fdr", "full_report.tsv"),
        sep='\t')
