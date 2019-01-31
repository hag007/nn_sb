import sys
sys.path.insert(0, '../')

import json
from matplotlib import style
from pandas._libs.parsers import k
import argparse

style.use("ggplot")
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
import utils.goids2gonames as goids2gonames
import shutil
from datasets_multithread_runner import run_dataset

from pandas.errors import EmptyDataError
import scipy



def plot_dist(pval ,algos_filter):



    fig, ax = plt.subplots(figsize=(12, 10))
    sns.distplot(pval, kde=False)
    plt.title("pval dist. -log10(pval) go terms. algo: {}".format("_".join(algos_filter)))
    real_n_term = np.size(pval)
    txt = "# terms={}".format(np.size(pval))
    fig.text(.5, .05, txt, ha='center')
    plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR,
                             "go_pval_dist_{}_{}.png".format(prefix.lower(), "_".join(algos_filter))))
    plt.clf()

def calc_dist(algos, datasets,is_plot=False,empirical_th=None):
    for cur_algo in algos:
        algos_filter = cur_algo

        df_go = pd.DataFrame(columns=['qval', 'pval'])
        df_go_pvals = pd.DataFrame()
        df_go_pvals.index.name="GO id"
        for cur_ds in datasets:
            go_results = [os.path.join(constants.OUTPUT_GLOBAL_DIR, cur_ds, cur_algo, cur_module) for cur_algo in
                          os.listdir(os.path.join(constants.OUTPUT_GLOBAL_DIR, cur_ds))
                          if os.path.isdir(
                    os.path.join(constants.OUTPUT_GLOBAL_DIR, cur_ds, cur_algo)) and cur_algo in algos_filter for
                          cur_module in os.listdir(os.path.join(constants.OUTPUT_GLOBAL_DIR, cur_ds, cur_algo)) if
                          "separated_modules" in cur_module]

            for cur in go_results:
                try:
                    df_go = pd.concat((df_go, pd.read_csv(cur, sep='\t')))
                    df_go_pvals = pd.concat((df_go_pvals, pd.read_csv(cur, sep='\t').set_index("GO id")['pval']), axis=1)
                except EmptyDataError:
                    pass
        df_go_pvals[df_go_pvals.isna()]=1
        df_go = df_go[df_go['qval'] < 0.05]
        if empirical_th:
            df_go = df_go[df_go['pval'].apply(lambda x:-np.log10(x)) > empirical_th]

        # df_go = df_go[~df_go.index.duplicated(keep='first')]
        pval = -np.log10(df_go["pval"].values)
        if np.size(pval) == 0:
            pval = np.array([0])

        if is_plot:
            plot_dist(pval, algos_filter=algos+datasets)

        return pval, df_go, df_go_pvals

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='args')
    parser.add_argument('--datasets', dest='datasets', default="SOC")
    parser.add_argument('--prefix', dest='prefix', default="GE")
    parser.add_argument('--algos', dest='algos', default="jactivemodules_greedy")
    parser.add_argument('--network', dest='network', default="dip.sif")
    parser.add_argument('--n_start', help="number of iterations (total n permutation is pf*(n_end-n_start))", dest='n_start', default=0)
    parser.add_argument('--n_end', help="number of iterations (total n permutation is pf*(n_end-n_start))", dest='n_end', default=1000)
    parser.add_argument('--network_permutations', dest='network_perm', default="false")
    parser.add_argument('--max_dist', help="takes max or all samples", dest='max_dist', default="true")
    parser.add_argument('--pval_measurement', help="how to calc pval. can be emp or beta", dest='pval_measurement', default="emp")
    parser.add_argument('--generate_plot', dest='generate_plot', default="false")
    parser.add_argument('--recalc_true_modules', dest='recalc_true_modules', default="false")

    args = parser.parse_args()

    datasets=args.datasets.split(",")
    algos=args.algos.split(",")
    prefix = args.prefix
    network_file_name = args.network
    n_start=int(args.n_start)
    n_end=int(args.n_end)
    max_dist=args.max_dist.lower()=="true"
    network_perm=args.network_perm.lower()=="true"
    generate_plot=args.generate_plot.lower()=="true"
    recalc_true_modules=args.recalc_true_modules.lower()=="true"
    pval_measurement=args.pval_measurement


    for dataset in datasets:


        score_method = constants.PREDEFINED_SCORE
        if prefix == "GE":
            score_method = constants.DEG_EDGER
            if dataset.startswith("IE"):
                score_method = constants.DEG_T


        for algo in algos:

            total_n_terms = []
            n_terms_filtered_in = []
            n_terms_2_oom = []
            diff_values = np.array([0])
            df_all_terms = pd.DataFrame()

            bad_iterations=[]
            for cur in range(n_start,n_end):


                random_ds = "{}_random_{}_{}_{}".format(prefix, dataset, algo,  cur)

                print "\ncur ds/n: {}/{}".format(random_ds, cur)

                try:
                    cur_pval, df_terms, df_pval_terms = calc_dist([algo], [random_ds.format(prefix, dataset)])
                    if max_dist:
                        df_pval_terms=df_pval_terms.min(axis=1).to_frame()

                    df_all_terms=pd.concat((df_all_terms, df_pval_terms), axis=1)
                except Exception:
                    print "cannot read iteration #{}".format(cur)
                    bad_iterations.append(cur)
                    pass

            print "total # bad iterations: {}".format(len(bad_iterations))
            print "bad iterations: {}".format(bad_iterations)


            # threshold difference
            df_all_terms = df_all_terms.apply(lambda x: -np.log10(x))
            df_all_terms[df_all_terms.isna()]=0
            df_all_terms_th=df_all_terms.apply(lambda r: np.percentile(r.values, 100, interpolation='higher'), axis=1)

            df_temp=df_all_terms.apply(lambda r: np.percentile(r.values, 100, interpolation='higher'), axis=1)
            df_temp=df_temp.to_frame()
            df_temp.columns=["pval"]
            df_temp["GO name"] = pd.Series(goids2gonames.get_go_names(list(df_temp.index)),
                                               index=df_temp.index)
            df_temp.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "emp_go_terms_th_values.tsv"), sep='\t', index_label="GO id")

            ## z statistic for each term
            df_agg_statistic = pd.concat((df_all_terms.mean(axis=1), df_all_terms.std(axis=1, ddof=0)), axis=1)
            df_agg_statistic.columns = ['mean', 'std']
            df_agg_statistic.loc[df_agg_statistic['std'].isna().values, 'std'] = 1



            if recalc_true_modules:
                if os.path.exists(os.path.join(constants.OUTPUT_GLOBAL_DIR,  "{}_{}".format(prefix, dataset), algo)):
                    shutil.rmtree(os.path.join(constants.OUTPUT_GLOBAL_DIR, "{}_{}".format(prefix, dataset), algo))
                run_dataset("{}_{}".format(prefix, dataset), score_method=score_method,
                            algos=[algo], network_file_name=network_file_name)
                if generate_plot:
                    plot_dist(pval, [algo, random_ds])


            ## empirical pvalues
            pval, df_go, df_agg_pval = calc_dist([algo], ["{}_{}".format(prefix, dataset)], is_plot=False, empirical_th=None)
            df_agg_pval=df_agg_pval.apply(lambda x: -np.log10(x))
            df_max_pval=df_agg_pval.max(axis=1).to_frame()
            print "total # real enriched terms: {}".format(len(df_max_pval.index))
            df_pvals = pd.DataFrame()
            df_counter=0
            df_agg_pval_filtered = df_agg_pval.loc[df_max_pval.iloc[:,0].values > 1.3, :]
            for index, row in df_agg_pval_filtered.iterrows():
                print "current dist_index: {}/{}".format(df_counter, len(df_agg_pval_filtered.index))
                df_counter+=1
                # if df_agg_pval.loc[index,"filtered_pval"] < 1.3:
                #     df_pvals.loc[index, "emp_pval"] =0
                #     df_pvals.loc[index, "dist_n_samples"] = 0
                #     df_pvals.loc[index, "sample_pos"] = 0
                row=df_all_terms.loc[index,:].values if index in df_all_terms.index  else np.array([0])
 
                # row=row[row!=0]
                if row.shape[0]==0:
                   row=np.array([0])
                print "shape: {}".format(row.shape)
                if pval_measurement=="beta":

                    bata_params = scipy.stats.chi2.fit(row)
                    real_value = df_max_pval.loc[index, :] if index in df_max_pval.index else 0
                    df_pvals.loc[index, "emp_pval"] = 1 - scipy.stats.chi2.cdf(real_value, *bata_params)
                    df_pvals.loc[index, "sample_pos"] = 0

                elif pval_measurement=="emp":
                    pos = np.size(row) - np.searchsorted(np.sort(row), df_max_pval.loc[index,:].iloc[0] if index in df_max_pval.index else 0, side='left')
                    df_pvals.loc[index,"emp_pval"]=pos/float(np.size(row))
                    df_pvals.loc[index, "sample_pos"] = pos

                df_pvals.loc[index, "dist_n_samples"] = str(list(row))


            for index, row in df_max_pval.iterrows():
                if index not in df_pvals.index:
                    df_pvals.loc[index, "emp_pval"] =0
                    df_pvals.loc[index, "dist_n_samples"] = 0
                    df_pvals.loc[index, "sample_pos"] = 0

            df_max_pval.columns=['pval']
            df_max_pval["GO name"] = pd.Series(goids2gonames.get_go_names(list(df_max_pval.index)),
                                               index=df_max_pval.index)


            ## z statistic for all terms together
            x = df_all_terms.loc[df_agg_pval.index.values]
            x[x.isna()]=0
            x=x.values.flatten()
            x = x[x != 0]
            joint_dist_mean = x.mean()
            joint_dist_std = x.std()
            if generate_plot:
                fig, ax = plt.subplots(figsize=(12, 10))
                ax.set_yscale('log')
                sns.distplot(x, kde=False)
                plt.title("joint dist real data terms")
                plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR,
                                         "joint_dist_real_terms_{}.png".format(algo)))
                plt.clf()


                x= (x - x.mean()) / x.std()
                fig, ax = plt.subplots(figsize=(12, 10))
                ax.set_yscale('log')
                sns.distplot(x, kde=False)
                plt.title("joint dist real data terms")
                plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR,
                                         "joint_dist_real_terms_normalized_{}.png".format(algo)))
                plt.clf()


            # filter terms
            print "total # of rows: {}".format(df_max_pval.shape)
            print "total # rows after 1 filter: {}".format(df_max_pval.loc[~df_max_pval.index.isin(df_all_terms_th.index.values)].shape)
            df_agg_pval_filtered=df_max_pval.loc[~df_max_pval.index.isin(df_all_terms_th.index.values) | (df_max_pval.iloc[:, 0] > df_all_terms_th.loc[df_max_pval.index]).values]
            print "total # row after both: {}".format(df_agg_pval_filtered.shape[0])

            # intermediate summary tables
            df_max_pval["filtered_in"]=pd.Series(df_max_pval.index.isin(df_agg_pval_filtered.index.values), df_max_pval.index)
            df_max_pval.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "emp_go_terms.tsv"), sep='\t', index_label="GO id")
            df_agg_pval_filtered.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "emp_filtered_go_terms.tsv"),
                                        sep='\t', index_label="GO id")

            # final summary tables
            df_diff=pd.concat((df_temp.rename(columns={"pval": "th_pval"}), df_max_pval.rename(columns={"pval": "filtered_pval"})) ,axis=1)
            df_diff.loc[df_diff["th_pval"].isna().values, "th_pval"] = 0
            df_diff.loc[df_diff["filtered_pval"].isna().values, "filtered_pval"] = 0
            df_diff["diff"]=df_diff["filtered_pval"]-df_diff["th_pval"]
            df_diff["emp_pval"]=df_pvals["emp_pval"]
            df_diff["dist_n_samples"] = df_pvals["dist_n_samples"]
            df_diff["sample_pos"] = df_pvals["sample_pos"]
            # zscore for each term separately
            # df_diff["zscore"] = (df_diff["filtered_pval"] - df_agg_statistic['mean'])/df_agg_statistic['std']
            df_diff["zscore"] = (df_diff["filtered_pval"] - joint_dist_mean) / joint_dist_std
            df_diff.loc[df_diff["zscore"].isna().values, "zscore"] = np.inf
            df_diff=df_diff.sort_values(ascending=False, by=['zscore', 'diff'])
            df_diff.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "emp_diff_{}_{}.tsv".format(dataset, algo)),
                                        sep='\t', index_label="GO id")

            if generate_plot:
                sns.distplot(df_diff["diff"].values, kde=False)
                fig, ax = plt.subplots(figsize=(12, 10))
                ax.set_yscale('log')
                ax.set_xlabel('-log10(pval)')
                ax.set_ylabel('log10(count)')
                sns.distplot(df_diff["diff"].values, kde=False)
                plt.title("diff dist")
                txt = "total # terms: {}. # terms passed filter={}. # terms passed 2 OOM: {}".format(len(df_diff.index), np.sum(df_diff["diff"].values >0), np.sum(df_diff["diff"].values>2))
                fig.text(.5, .05, txt, ha='center')
                plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR,
                                         "diff_dist_{}.png".format(algo)))
                plt.clf()

            total_n_terms.append(len(df_diff.index))
            n_terms_filtered_in.append(np.sum(df_diff["diff"].values > 0))
            diff_values=np.append(diff_values, df_diff["diff"][df_diff["diff"]>0].values)
            n_terms_2_oom.append(np.sum(df_diff["diff"].values > 2))

            print(total_n_terms)
            print(np.mean(total_n_terms))
            print(np.std(total_n_terms))
            print(n_terms_filtered_in)
            print(np.mean(n_terms_filtered_in))
            print(n_terms_filtered_in)
            print(np.mean(n_terms_filtered_in))
            print(np.std(n_terms_filtered_in))
            print(n_terms_2_oom)
            print(np.mean(n_terms_2_oom))
            print(np.std(n_terms_2_oom))



            # fig, ax = plt.subplots(figsize=(12, 10))
            # ax.set_yscale('log')
            # sns.distplot(diff_values, kde=False)
            # plt.title("diff dist")
            # txt = "total # diff: {}. 90'th percentile ={}. 95'th percentile : {}".format(np.size(diff_values), round(np.percentile(diff_values, 90, interpolation='higher'),2), round(np.percentile(diff_values, 95, interpolation='higher'),2))
            # fig.text(.5, .05, txt, ha='center')
            # plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR,
            #                          "diff_positive_agg_dist_{}_{}.png".format(dataset,algo)))
            # plt.clf()





