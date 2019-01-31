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
import scipy
from scipy.optimize import least_squares
import random
from statsmodels.sandbox.stats.multicomp import fdrcorrection0

def calc_emp_pval(cur_rv, cur_dist):
    pos = np.size(cur_dist) - np.searchsorted(np.sort(cur_dist), cur_rv, side='left')


    return pos / float(np.size(cur_dist))



def main(dataset_test ="SOC2", algo_test ="hotnet2_5000", algo_sample = None, dataset_sample = None, n_dist_samples = 300, n_test_samples = 300, limit = 10000):
    single_dataset = False
    if dataset_sample is None:
        single_dataset=True
        dataset_sample=dataset_test
    if algo_sample is None:
        algo_sample = algo_test

    output_md = pd.read_csv(
        os.path.join(constants.OUTPUT_GLOBAL_DIR, "emp_fdr", "{}_MAX".format(dataset_sample), "emp_diff_{}_{}_md.tsv".format(dataset_sample, algo_sample)),
        sep='\t', index_col=0).dropna()


    n_genes_pvals=output_md.loc[np.logical_and.reduce([output_md["n_genes"].values > 5, output_md["n_genes"].values < 500]), "filtered_pval"].values

    print "total n_genes with pval:{}/{}".format(np.size(n_genes_pvals), 7435)
    n_genes_pvals=np.append(n_genes_pvals,np.zeros(7435-np.size(n_genes_pvals)))
    n_genes_pvals = [10**(-x) for x in n_genes_pvals]
    fdr_results = fdrcorrection0(n_genes_pvals, alpha=0.05, method='indep', is_sorted=False)
    n_hg_true = len([cur for cur in fdr_results[0] if cur == True])
    HG_CUTOFF=np.sort(n_genes_pvals)[n_hg_true-1]
    print "HG cutoff: {}".format(HG_CUTOFF)

    output_md = output_md.loc[np.logical_and.reduce(
        [output_md["n_genes"].values > 5, output_md["n_genes"].values < 500,
         output_md["filtered_pval"].values > HG_CUTOFF]), :]

    output = pd.read_csv(
        os.path.join(constants.OUTPUT_GLOBAL_DIR, "emp_fdr", "{}_MAX".format(dataset_sample), "emp_diff_{}_{}.tsv".format(dataset_sample, algo_sample)),
        sep='\t', index_col=0).dropna()
    output = output.loc[output_md.index.values, :]
    counter = 0
    emp_dists = []
    emp_pvals = []

    i_dist = None
    i_test = None
    if single_dataset:
        n_total_samples=len(output.iloc[0].loc["dist_n_samples"][1:-1].split(", "))
        i_choice=np.random.choice(n_total_samples, n_dist_samples + n_test_samples)
        i_dist=i_choice[:n_dist_samples]
        i_test=i_choice[n_dist_samples:]
    else:
        n_total_samples = len(output.iloc[0].loc["dist_n_samples"][1:-1].split(", "))
        i_dist = np.random.choice(n_total_samples, n_dist_samples)



    for index, cur in output.iterrows():
        if counter == limit: break
        # print "cur iteration in real data: {}/{}".format(counter, len(output.index))
        pval = np.array([float(x) for x in cur["dist_n_samples"][1:-1].split(", ")])[i_dist]

        emp_pvals.append(calc_emp_pval(cur["filtered_pval"], pval))
        emp_dists.append(pval)
        counter += 1
    dist_name = "emp"
    df_dists = pd.DataFrame(index=output.index)
    df_dists["emp"] = pd.Series(emp_dists, index=output.index[:limit])

    fdr_results = fdrcorrection0(emp_pvals, alpha=0.05, method='indep', is_sorted=False)
    go_ids_results.append(output.index.values[fdr_results[0]])
    go_names_results.append(output["GO name"].values[fdr_results[0]])
    n_emp_true = len([cur for cur in fdr_results[0] if cur == True])
    BH_TH = np.sort(emp_pvals)[n_emp_true - 1]

    output = pd.read_csv(
        os.path.join(constants.OUTPUT_GLOBAL_DIR, "emp_fdr", "{}_MAX".format(dataset_test), "emp_diff_{}_{}.tsv".format(dataset_test, algo_test)),
        sep='\t', index_col=0).dropna()
    output = output.loc[output_md.index.values, :]

    if i_test is None:
        n_total_samples = len(output.iloc[0].loc["dist_n_samples"][1:-1].split(", "))
        i_test = np.random.choice(n_total_samples, n_test_samples)





    counter = 0
    pvals = []
    for index, cur in output.iterrows():
        if counter == limit: break
        # print "cur iteration in reading random data: {}/{}".format(counter, len(output.index))
        pvals.append(np.array([float(x) for x in cur["dist_n_samples"][1:-1].split(", ")])[i_test])
        counter += 1
    pvals = np.array(pvals)

    # dist_name="chi2"
    random_above_th = []
    counter = 0
    for cur_random in pvals.T:
        if counter == limit: break
        # print "cur iteration in testing random data: {}/{}".format(counter, pvals.T.shape[0])
        random_above_th.append(sum([a < BH_TH for a in [calc_emp_pval(cur_rv, cur_dist) for cur_rv, cur_dist in
                                                        zip(cur_random, df_dists[dist_name])]]))
        counter += 1

    random_above_th = np.array(random_above_th)
    fig, ax = plt.subplots()
    ax.set_yscale('log')
    ax.set_xlabel("# sig. terms")
    ax.set_ylabel("# randomized DSs")
    plt.hist(random_above_th, bins=np.max(random_above_th) + 1)
    plt.title("# non-zero perm={}/{}, total sum={}\nqval={}, real DS' true terms: {}/{}".format(
        len(random_above_th[random_above_th != 0]), len(random_above_th),
        np.sum(random_above_th[random_above_th != 0]),
        BH_TH, n_hg_true, len(output_md.index)
        ))
    plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "fdr_dist_{}_{}.png".format(dataset_test, algo_test)))
    random_above_th = np.array(random_above_th)
    print random_above_th[random_above_th != 0]
    print "BH cutoff: {} # true terms passed BH cutoff: {}".format(BH_TH, n_emp_true)
    return BH_TH, n_emp_true, HG_CUTOFF, n_hg_true


if __name__ == "__main__":

    dataset_sample = "TNFa_22"
    algo_sample = "bionet_1000" # "hotnet2_5000"
    dataset_test = "TNFa_22"
    algo_test = "bionet_1000"

    dataset_sample = None
    algo_sample = None



    n_dist_samples = 300
    n_test_samples = 300

    l_emp_cutoff=[]
    l_n_emp_true=[]
    l_hg_cutoff=[]
    l_n_hg_true=[]

    for cur in range(10):
        BH_TH, n_emp_true, HG_CUTOFF, n_hg_true = main(dataset_test= dataset_test, algo_test= algo_test, n_dist_samples = n_dist_samples, n_test_samples = n_test_samples, dataset_sample=dataset_sample, algo_sample=algo_sample)
        l_emp_cutoff.append(str(BH_TH))
        l_n_emp_true.append(str(n_emp_true))
        l_hg_cutoff.append(str(HG_CUTOFF))
        l_n_hg_true.append(str(n_hg_true))

    fdr_test_summary=pd.DataFrame()
    fdr_test_summary["emp_cutoff"]=l_emp_cutoff
    fdr_test_summary["n_emp_true"]=l_n_emp_true
    fdr_test_summary["hg_cutoff"]=l_hg_cutoff
    fdr_test_summary["n_hg_true"]=l_n_hg_true

    fdr_test_summary.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "emp_fdr", "fdr_test_summary_{}_5000_300_300.tsv".format("bionet")), sep='\t')