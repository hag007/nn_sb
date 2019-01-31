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

dict_result, go2geneids, geneids2go, entrez2ensembl = go_hierarcies.build_hierarcy(
            roots=['GO:0008150'])
vertices=dict_result.values()[0]['vertices']


def get_all_genes_for_term(vertices, cur_root, term, in_subtree):

    in_subtree= in_subtree or term==cur_root
    all_genes = set()
    if in_subtree:
        all_genes.update(go2geneids[cur_root])
    for cur_child in vertices[cur_root]["obj"].children:
        all_genes.update(get_all_genes_for_term(vertices, cur_child.id, term, in_subtree))
    return all_genes


def main(dataset="SOC", algo="jactivemodules_sa", csv_file_name=os.path.join(constants.OUTPUT_GLOBAL_DIR, "emp_fdr","{dataset}_MAX/emp_diff_{dataset}_{algo}.tsv" )):

    csv_file_name=csv_file_name.format(dataset=dataset, algo=algo)
    df=None
    try:
        df=pd.read_csv(csv_file_name,sep='\t',index_col=0)
    except:
        return None
    df=df.dropna()
    n_genes=[]
    depth=[]
    for i, cur_go_id in enumerate(df.index.values):
        # print "current i: {}/{}".format(i, len(df.index.values))
        n_genes.append(len(get_all_genes_for_term(vertices, cur_go_id, cur_go_id, cur_go_id==cur_go_id)))
        depth.append(dict_result.values()[0]['vertices'][cur_go_id]['D'])

    df["n_genes"]=pd.Series(n_genes, index=df.index)
    df["depth"]=pd.Series(depth, index=df.index)




    # HG_CUTOFF=1.92 # SOC: 2.109
    n_genes_pvals=df.loc[np.logical_and.reduce([df["n_genes"].values > 5, df["n_genes"].values < 500]), "filtered_pval"].values

    print "total n_genes with pval:{}/{}".format(np.size(n_genes_pvals), 7435)
    n_genes_pvals=np.append(n_genes_pvals,np.zeros(7435-np.size(n_genes_pvals)))
    n_genes_pvals = [10**(-x) for x in n_genes_pvals]
    fdr_results = fdrcorrection0(n_genes_pvals, alpha=0.05, method='indep', is_sorted=False)
    true_counter = len([cur for cur in fdr_results[0] if cur == True])
    HG_CUTOFF=-np.log10(n_genes_pvals[true_counter])
    print "cutoff: {}".format(HG_CUTOFF)

    df_filtered_in=df.loc[np.logical_and.reduce([df["n_genes"].values > 5, df["n_genes"].values < 500, df["filtered_pval"].values > HG_CUTOFF]), :]
    df_filtered_in["emp_rank"]=df_filtered_in["emp_pval"].rank(ascending=1)
    df_filtered_in["hg_rank"]=df_filtered_in["filtered_pval"].rank(ascending=0)
    df_filtered_out=df.loc[~np.logical_and.reduce([df["n_genes"].values > 5, df["n_genes"].values < 500, df["filtered_pval"].values > HG_CUTOFF]), :]

    counter=0
    oob_params_counter=0
    for index, cur in df_filtered_in.iterrows():

        # if counter <109:
        # if cur["GO name"]!="cell cycle phase transition":
        #     counter+=1
        #     continue





        print "{}/{}".format(counter, len(df_filtered_in.index))


        pval = np.array([float(x) for x in cur["dist_n_samples"][1:-1].split(", ")])

        # pval = np.append(pval, np.zeros(8000-np.size(pval)))
        #pval[pval != 0]

        # pval = pval-np.min(pval)

        if len(pval) < 1000:
            print "too few samples in {}, {} : {}".format(dataset, algo, len(pval))
            exit(1)
        else:
            pval=pval[:1000]
            pos = np.size(pval) - np.searchsorted(np.sort(pval), cur["filtered_pval"], side='left')

            df_filtered_in.loc[index, "emp_pval"]=pos / float(np.size(pval))

        ###
        # sigma2 = Parameter("sigma2", value=1, min=0.5, max=5.0)
        # mu = Parameter("mu", value=max(np.max(pval) - 1, 5.0), min=5.0, max=max(np.max(pval), 5.0))
        # V = Parameter("V", value=0.5, min=0.0001, max=1)
        # A = Parameter("A", value=1.0, min=0.0001, max=10)
        # B = Parameter("B", value=0.0, min=0.0, max=max(min(pval) + 1, 1))
        #
        # x = Variable()
        #
        # model = Add(Mul(V, Piecewise(
        #     (Mul(exp(-Mul(Add(x, -B), 1 / A)), 1 / A), GreaterThan(Add(x, -B), 0)), (1e-09, True))),
        #             Mul(Add(1, -V), Piecewise(((1 / (sqrt(2 * pi * np.abs(sigma2)))) * exp(
        #                 -(x - mu) ** 2 / (2 * np.abs(sigma2))), GreaterThan(Add(x, -B), 0)), (1e-09, True))))

        ###


        # beta_params = scipy.stats.beta.fit(pval)
        # chi2_params = scipy.stats.chi2.fit(pval)


        # param_counter = 0
        # opt_objective_value = np.inf
        # opt_mixture_params = None
        # optional_mixture_params = [{"v": 0.5}, {"v": 0.75}, {"v": 1.0}, {"v": 0.25}, {"v": 0.001}]
        # prcs = []
        # for i, cur_params in enumerate(optional_mixture_params):
        #     prcs.append(subprocess.Popen("{}/../python27/bin/python {} --data {} --v {}"
        #                                  .format(constants.dir_path, "mixture_dist_sa.py", ",".join(pval.astype(np.str)),
        #                                          cur_params["v"]), shell=True,
        #                                  stdout=subprocess.PIPE, cwd=os.path.join(constants.dir_path, "utils")))
        #
        # param_counter = 0
        # opt_objective_value = np.inf
        # optional_mixture_params = [{"v": 0.5}, {"v": 0.75}, {"v": 1.0}, {"v": 0.25}, {"v": 0.0001}]
        # mixture_result = None
        # opt_mixture_params = None
        # manager = multiprocessing.Manager()
        # return_list = manager.list()
        # prcs = []
        # for i, cur_params in enumerate(optional_mixture_params):
        #     prcs.append(subprocess.Popen("{}/../python27/bin/python {} --data {} --v {}"
        #                                  .format(constants.dir_path, "mixture_dist_sa.py", ",".join(pval.astype(np.str)),
        #                                          cur_params["v"]), shell=True,
        #                                  stdout=subprocess.PIPE, cwd=os.path.join(constants.dir_path, "utils")))
        # for i, prc in enumerate(prcs):
        #
        #     output = prc.stdout.read()
        #
        #     output = output.split("\n")[-1]
        #     if output == "None":
        #         continue
        #     objective_value = float(output.split("#")[1])
        #     mixture_params = {cur.split(":")[0]: float(cur.split(":")[1]) for cur in output.split("#")[0].split(",")}
        #
        #     # mixture_params['k']< k.min or
        #     if mixture_params['A'] < A.min or \
        #             mixture_params['B'] < B.min or mixture_params['sigma2'] < sigma2.min or \
        #             mixture_params['mu'] < mu.min:
        #         print "params out of bounds: {}".format(mixture_params)
        #         oob_params_counter += 1
        #         continue
        #
        #     print("cur params: {}, cur objective: {}".format(mixture_params, objective_value))
        #     if opt_objective_value > np.abs(objective_value):
        #         opt_objective_value = np.abs(objective_value)
        #         opt_mixture_params = mixture_params



        # param_counter=0
        # max_objective_values = -1
        # optional_mixture_params = [{"v" : 1.0},{"v" : 0.75},{"v" : 0.5}]
        # mixture_results = []
        # mixture_result=None
        # opt_mixture_params=None
        # manager = multiprocessing.Manager()
        # return_list = [] # manager.list()
        # prcs=[]
        # for i, cur_params in enumerate(optional_mixture_params):
        # #     prcs.append(Process(target=utils.mixture_dist.fit, args=[pval, cur_params['v'], return_list]))
        # #     prcs[-1].start()
        # #
        # # for prc in prcs:
        # #     prc.join()
        #     return_list.append(utils.mixture_dist.fit(pval, cur_params['v'], [])) # return_list
        # for mixture_result in return_list:
        #     if mixture_result != None:
        #         model, mixture_params, objective_value = mixture_results
        #         if max_objective_values < objective_value:
        #             max_objective_values=objective_value
        #             opt_mixture_params=mixture_params
        #

        # df_filtered_in.loc[index, "mixture_pval"]=(1-opt_mixture_params['V']) * scipy.stats.norm(opt_mixture_params['mu'], opt_mixture_params['sigma2']).sf(cur["filtered_pval"]) + opt_mixture_params['V'] * scipy.stats.expon(loc=opt_mixture_params['B'], scale=opt_mixture_params['A']).sf(cur["filtered_pval"])
        # df_filtered_in.loc[index, "beta_pval"]=1-scipy.stats.beta(*beta_params).cdf(cur["filtered_pval"])
        # df_filtered_in.loc[index, "chi2_pval"]=1-scipy.stats.chi2(*chi2_params).cdf(cur["filtered_pval"])

        print cur["GO name"]
        # print opt_mixture_params
        # print df_filtered_in.loc[index, "beta_pval"]
        # print df_filtered_in.loc[index, "chi2_pval"]
        print df_filtered_in.loc[index, "filtered_pval"]
        # print df_filtered_in.loc[index, "mixture_pval"]
        counter+=1

    # df_filtered_in["chi2_rank"]=df_filtered_in["chi2_pval"].rank(ascending=1)
    # df_filtered_in["beta_rank"]=df_filtered_in["beta_pval"].rank(ascending=1)
    # df_filtered_in["mixture_rank"]=df_filtered_in["mixture_pval"].rank(ascending=1)

    df_filtered_in["emp_rank"] = df_filtered_in["emp_pval"].rank(ascending=1)
    df_filtered_in["hg_rank"] = df_filtered_in["filtered_pval"].rank(ascending=0)


    pvals_corrected=df_filtered_in["emp_pval"].values
    fdr_results = fdrcorrection0(pvals_corrected , alpha=0.05, method='indep', is_sorted=False)
    true_counter = len([cur for cur in fdr_results[0] if cur == True])
    emp_cutoff=np.sort(pvals_corrected)[true_counter-1] if true_counter > 0 else 0
    print "emp true hypothesis: {} (emp cutoff: {}, n={})".format(true_counter, emp_cutoff,len(pval))

    # fdr_results = fdrcorrection0(df_filtered_in["chi2_pval"].values, alpha=0.05, method='indep', is_sorted=False)
    # true_counter = len([cur for cur in fdr_results[0] if cur == True])
    # print "chi2 true hypothesis: {}".format(true_counter)

    # fdr_results = fdrcorrection0(df_filtered_in["beta_pval"].values, alpha=0.05, method='indep', is_sorted=False)
    # true_counter = len([cur for cur in fdr_results[0] if cur == True])
    # print "beta true hypothesis: {}".format(true_counter)

    # fdr_results = fdrcorrection0(df_filtered_in["mixture_pval"].values, alpha=0.05, method='indep', is_sorted=False)
    # true_counter = len([cur for cur in fdr_results[0] if cur == True])
    # print "mixture true hypothesis: {}".format(true_counter)


    pd.concat((df_filtered_in, df_filtered_out), axis=0).loc[df["filtered_pval"].values > 0, :][["GO name", "filtered_pval", "emp_pval", "hg_rank", "emp_rank",  "n_genes", "depth", "diff", "th_pval"]].to_csv(csv_file_name[:-4]+"_md.tsv",sep='\t') #

    return len(df_filtered_in.index), true_counter, HG_CUTOFF, emp_cutoff


    # "emp_pval", "beta_pval", "chi2_pval", "mixture_pval", "hg_rank", "emp_rank", "beta_rank", "chi2_rank", "mixture_rank",


if __name__ == "__main__":
    csv_file_name=os.path.join(constants.OUTPUT_GLOBAL_DIR,"emp_fdr","SOC2_MAX", "emp_diff_{dataset}_{algo}.tsv")
    main(dataset="SOC2", algo="jactivemodules_greedy_5000", csv_file_name=csv_file_name)