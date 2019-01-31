import multiprocessing
from multiprocessing import Process
import go_hierarcies
import pandas as pd
import numpy as np
import scipy
from statsmodels.sandbox.stats.multicomp import fdrcorrection0
import utils.mixture_dist

from symfit import Parameter, Variable, Fit, sqrt, pi, Equality, Abs, GreaterThan
from sympy.functions.elementary.exponential import exp
from sympy.functions.elementary.piecewise import Piecewise
from sympy.core.add import Add
from sympy.core.mul import Mul
from sympy.core.power import Pow
from sympy.functions.special.gamma_functions import gamma
from symfit.core.objectives import LogLikelihood
import subprocess
import constants
import os
from utils.ensembl2entrez import entrez2ensembl_convertor
dict_result, go2geneids, geneids2go, entrez2ensembl = go_hierarcies.build_hierarcy(
            roots=['GO:0008150'])
vertices=dict_result.values()[0]['vertices']


def mean_difference(row, dataset_data, classes_data):
    try:
        return dataset_data.loc[entrez2ensembl_convertor(get_all_genes_for_term(vertices, row["index"], row["index"],
                                                                             True)), classes_data == 2].dropna().values.mean() - \
        dataset_data.loc[entrez2ensembl_convertor(get_all_genes_for_term(vertices, row["index"], row["index"],
                                                                             True)), classes_data == 1].dropna().values.mean()
    except:
        print "no gene were found for {}, {} (pval={})".format(row["index"], row["GO name"],
                                                                        row["filtered_pval"])



def get_all_genes_for_term(vertices, cur_root, term, in_subtree):

    in_subtree= in_subtree or term==cur_root
    all_genes = set()
    if in_subtree:
        all_genes.update(go2geneids[cur_root])
    for cur_child in vertices[cur_root]["obj"].children:
        all_genes.update(get_all_genes_for_term(vertices, cur_child.id, term, in_subtree))
    return all_genes


def main(dataset="SOC", algo="jactivemodules_sa", csv_file_name=os.path.join(constants.OUTPUT_GLOBAL_DIR, "emp_fdr","{dataset}_MAX/emp_diff_{dataset}_{algo}_md.tsv" )):

    df_hg_cutoffs=pd.read_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "emp_fdr", "hg_cutoffs_summary.tsv"), sep='\t')
    df_emp_cutoffs=pd.read_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "emp_fdr", "emp_cutoffs_summary.tsv"), sep='\t')
    df_n_terms=pd.read_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "emp_fdr", "n_terms_summary.tsv"), sep='\t')

    print "dataset: {}".format(dataset)

    dataset_data=pd.read_csv(os.path.join(constants.DATASETS_DIR, "GE_{}".format(dataset),"data", "ge.tsv"), sep='\t', index_col=0)
    classes_data=np.array(file(os.path.join(constants.DATASETS_DIR, "GE_{}".format(dataset), "data", "classes.tsv")).readlines()[0].strip().split("\t")).astype(np.int)



    csv_file_name=csv_file_name.format(dataset=dataset, algo=algo)
    df=None
    try:
        df=pd.read_csv(csv_file_name,sep='\t',index_col=0)
    except:
        return None
    df=df.dropna()


    n_genes_pvals=df.loc[np.logical_and.reduce([df["n_genes"].values > 5, df["n_genes"].values < 500]), "filtered_pval"].values

    print "total n_genes with pval:{}/{}".format(np.size(n_genes_pvals), 7435)
    n_genes_pvals=np.append(n_genes_pvals,np.zeros(7435-np.size(n_genes_pvals)))
    n_genes_pvals = [10**(-x) for x in n_genes_pvals]
    fdr_results = fdrcorrection0(n_genes_pvals, alpha=0.05, method='indep', is_sorted=False)
    true_counter = len([cur for cur in fdr_results[0] if cur == True])
    HG_CUTOFF=np.sort(n_genes_pvals)[true_counter-1]
    print "HG cutoff: {}".format(HG_CUTOFF)

    n_genes_pvals = df["emp_pval"].values
    fdr_results = fdrcorrection0(n_genes_pvals, alpha=0.05, method='indep', is_sorted=False)
    true_counter = len([cur for cur in fdr_results[0] if cur == True])
    EMP_CUTOFF = np.sort(n_genes_pvals)[true_counter-1] if true_counter>0 else 0
    print "emp cutoff: {}".format(EMP_CUTOFF)

    df_filtered_in=df.loc[np.logical_and.reduce([df["n_genes"].values > 5, df["n_genes"].values < 500, df["filtered_pval"].values > HG_CUTOFF]), :]
    df_filtered_out=df.loc[~np.logical_and.reduce([df["n_genes"].values > 5, df["n_genes"].values < 500, df["filtered_pval"].values > HG_CUTOFF]), :]

    counter=0
    oob_params_counter=0




    pvals_corrected=df_filtered_in["emp_pval"].values
    fdr_results = fdrcorrection0(pvals_corrected , alpha=0.05, method='indep', is_sorted=False)
    true_counter = len([cur for cur in fdr_results[0] if cur == True])
    emp_cutoff=np.sort(pvals_corrected)[true_counter-1] if true_counter > 0 else 0
    print "emp true hypothesis: {} (emp cutoff: {})".format(true_counter, emp_cutoff)

    # fdr_results = fdrcorrection0(df_filtered_in["chi2_pval"].values, alpha=0.05, method='indep', is_sorted=False)
    # true_counter = len([cur for cur in fdr_results[0] if cur == True])
    # print "chi2 true hypothesis: {}".format(true_counter)

    # fdr_results = fdrcorrection0(df_filtered_in["beta_pval"].values, alpha=0.05, method='indep', is_sorted=False)
    # true_counter = len([cur for cur in fdr_results[0] if cur == True])
    # print "beta true hypothesis: {}".format(true_counter)

    # fdr_results = fdrcorrection0(df_filtered_in["mixture_pval"].values, alpha=0.05, method='indep', is_sorted=False)
    # true_counter = len([cur for cur in fdr_results[0] if cur == True])
    # print "mixture true hypothesis: {}".format(true_counter)


    pd.concat((df_filtered_in, df_filtered_out), axis=0).loc[df["filtered_pval"].values > 0, :][["GO name", "filtered_pval", "emp_pval", "hg_rank", "emp_rank", "passed_fdr", "mean_difference", "n_genes", "depth", "diff", "th_pval"]].to_csv(csv_file_name[:-4]+"_x.tsv",sep='\t') #

    return len(df_filtered_in.index), true_counter, HG_CUTOFF, emp_cutoff


    # "emp_pval", "beta_pval", "chi2_pval", "mixture_pval", "hg_rank", "emp_rank", "beta_rank", "chi2_rank", "mixture_rank",


if __name__ == "__main__":
    csv_file_name=os.path.join(constants.OUTPUT_GLOBAL_DIR,"emp_fdr","SOC_MAX", "emp_diff_{dataset}_{algo}_md.tsv")
    main(dataset="SOC2", algo="jactivemodules_greedy_5000", csv_file_name=csv_file_name)