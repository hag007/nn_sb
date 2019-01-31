import sys
sys.path.insert(0, '../')

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
from utils.param_builder import build_gdc_params
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import utils.goids2gonames as goids2gonames
import shutil
from datasets_multithread_runner import run_dataset
from utils.go_pval_dist import create_random_ds
from utils.go_pval_dist import create_permuted_network

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
    # prefix="GE"
    # datasets = ["{}_TNFa_2".format(prefix)]
    # algos = np.array(["jactivemodules_greedy"])
    #
    # calc_dist(algos, datasets)
    network_file_name="dip"
    prefix="GE"
    datasets = ["SOC"] # , "IEM", "IEN", "HC12", "MCF7_2", "TNFa_2"   alzheimers, schizophrenia
    algos = ["jactivemodules_greedy"] # "bionet"

    summary = []
    for dataset in datasets:


        score_method = constants.PREDEFINED_SCORE
        if prefix == "GE":
            score_method = constants.DEG_EDGER
            if dataset.startswith("IE"):
                score_method = constants.DEG_T

        total_n_terms=[]
        n_terms_filtered_in=[]
        n_terms_2_oom=[]
        summary = []
        diff_values=np.array([0])
        for cur_real_from_rand in range(0, 1):
            df_all_terms = pd.DataFrame()
            cur_real_ds= "{}_{}".format(prefix, dataset) # "{}_random_{}_{}".format(prefix, dataset, cur_real_from_rand)

            for algo in algos:
                pval = np.array([])
                random_ds = "{}_random_{}".format(prefix, dataset)
                bad_iterations=[]
                for cur in range(650):
                    # if cur==cur_real_from_rand: continue

                    print "\ncur iteration: {}".format(cur)
                    random_ds=create_random_ds(prefix, "{}_{}".format(prefix, dataset),cur)
                    permuted_network_file_name="dip" # _perm
                    # if cur==0:
                    #     permuted_network_file_name=create_permuted_network(network_file_name=network_file_name)
                    # run_dataset(random_ds, score_method=score_method,
                    #             algos=[algo], network_file_name=permuted_network_file_name)
                    try:
                        cur_pval, df_terms, df_pval_terms = calc_dist([algo], [random_ds.format(prefix, dataset)])
                        pval = np.append(pval, cur_pval)
                        df_all_terms=pd.concat((df_all_terms, df_pval_terms), axis=1)
                    except Exception:
                        print "cannot read iteration #{}".format(cur)
                        bad_iterations.append(cur)
                        pass

                print "total # bad iterations: {}".format(len(bad_iterations))
                print "bad iterations: {}".format(bad_iterations)
        #
                df_all_terms = df_all_terms.apply(lambda x: -np.log10(x))
                df_all_terms[df_all_terms.isna()]=0
                df_all_terms_th=df_all_terms.apply(lambda r: np.percentile(r.values, 100, interpolation='higher'), axis=1)
                df_temp=df_all_terms.apply(lambda r: np.percentile(r.values, 100, interpolation='higher'), axis=1)
                df_temp=df_temp.to_frame()
                df_temp.columns=["pval"]
                df_temp["GO name"] = pd.Series(goids2gonames.get_go_names(list(df_temp.index)),
                                                   index=df_temp.index)

                df_temp.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "emp_go_terms_th_values.tsv"), sep='\t', index_label="GO id")

                df_agg_statistic = pd.concat((df_all_terms.mean(axis=1), df_all_terms.std(axis=1, ddof=0)), axis=1)
                df_agg_statistic.columns = ['mean', 'std']
                df_agg_statistic.loc[df_agg_statistic['std'].isna().values, 'std'] = 1


                # if os.path.exists(os.path.join(constants.OUTPUT_GLOBAL_DIR,  "{}_{}".format(prefix, dataset), algo)):
                #     shutil.rmtree(os.path.join(constants.OUTPUT_GLOBAL_DIR, "{}_{}".format(prefix, dataset), algo))
                # run_dataset("{}_{}".format(prefix, dataset), score_method=score_method,
                #             algos=[algo], network_file_name=network_file_name)

                # plot_dist(pval, [algo, random_ds])

                pval, df_go, df_agg_pval = calc_dist([algo], [cur_real_ds], is_plot=False, empirical_th=None)
                df_agg_pval=df_agg_pval.apply(lambda x: -np.log10(x))
                df_max_pval=df_agg_pval.max(axis=1).to_frame()
                df_pvals = pd.DataFrame()
                beta_iteration=0
                total_beta_iterations=len(df_all_terms.index)
                for index, row in df_all_terms.iterrows():
                    print "current beta iteration: {}/{}".format(beta_iteration,total_beta_iterations)
                    beta_iteration+=1
                    bata_params = scipy.stats.beta.fit(row.values)
                    real_value=df_max_pval.loc[index, :] if index in df_max_pval.index else 0
                    df_pvals.loc[index,"emp_pval"]=1-scipy.stats.beta.cdf(real_value, *bata_params)
                    df_pvals.loc[index, "dist_n_samples"] = np.size(row.values)

                for index, row in df_max_pval.iterrows():
                    if index not in df_pvals.index:
                        df_pvals.loc[index, "emp_pval"] =0
                        df_pvals.loc[index, "dist_n_samples"] = 0

                df_max_pval.columns=['pval']
                df_max_pval["GO name"] = pd.Series(goids2gonames.get_go_names(list(df_max_pval.index)),
                                                   index=df_max_pval.index)


                x = df_all_terms.loc[df_agg_pval.index.values]
                x[x.isna()]=0
                x=x.values.flatten()
                x = x[x != 0]
                # fig, ax = plt.subplots(figsize=(12, 10))
                # ax.set_yscale('log')
                # sns.distplot(x, kde=False)
                # plt.title("joint dist real data terms")
                # plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR,
                #                          "joint_dist_real_terms_{}.png".format(algo)))
                # plt.clf()

                joint_dist_mean=x.mean()
                joint_dist_std = x.std()
                # x= (x - x.mean()) / x.std()
                # fig, ax = plt.subplots(figsize=(12, 10))
                # ax.set_yscale('log')
                # sns.distplot(x, kde=False)
                # plt.title("joint dist real data terms")
                # plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR,
                #                          "joint_dist_real_terms_normalized_{}.png".format(algo)))
                # plt.clf()




                print "total # of rows: {}".format(df_max_pval.shape)
                print "total # rows after 1 filter: {}".format(df_max_pval.loc[~df_max_pval.index.isin(df_all_terms_th.index.values)].shape)
                df_agg_pval_filtered=df_max_pval.loc[~df_max_pval.index.isin(df_all_terms_th.index.values) | (df_max_pval.iloc[:, 0] > df_all_terms_th.loc[df_max_pval.index]).values]
                print "total # row after both: {}".format(df_agg_pval_filtered.shape[0])

                df_max_pval["filtered_in"]=pd.Series(df_max_pval.index.isin(df_agg_pval_filtered.index.values), df_max_pval.index)
                df_max_pval.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "emp_go_terms.tsv"), sep='\t', index_label="GO id")
                df_agg_pval_filtered.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "emp_filtered_go_terms.tsv"),
                                            sep='\t', index_label="GO id")
                df_diff=pd.concat((df_temp.rename(columns={"pval": "th_pval"}), df_max_pval.rename(columns={"pval": "filtered_pval"})) ,axis=1)
                df_diff.loc[df_diff["th_pval"].isna().values, "th_pval"] = 0
                df_diff.loc[df_diff["filtered_pval"].isna().values, "filtered_pval"] = 0
                df_diff["diff"]=df_diff["filtered_pval"]-df_diff["th_pval"]
                df_diff["emp_pval"]=df_pvals["emp_pval"]
                df_diff["dist_n_samples"] = df_pvals["dist_n_samples"]
                # df_diff["zscore"] = (df_diff["filtered_pval"] - df_agg_statistic['mean'])/df_agg_statistic['std']
                df_diff["zscore"] = (df_diff["filtered_pval"] - joint_dist_mean) / joint_dist_std
                df_diff.loc[df_diff["zscore"].isna().values, "zscore"] = np.inf
                df_diff=df_diff.sort_values(ascending=False, by=['zscore', 'diff'])
                df_diff.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "emp_diff.tsv"),
                                            sep='\t', index_label="GO id")

                # sns.distplot(df_diff["diff"].values, kde=False)
                # fig, ax = plt.subplots(figsize=(12, 10))
                # ax.set_yscale('log')
                # ax.set_xlabel('-log10(pval)')
                # ax.set_ylabel('log10(count)')
                # sns.distplot(df_diff["diff"].values, kde=False)
                # plt.title("diff dist")
                # txt = "total # terms: {}. # terms passed filter={}. # terms passed 2 OOM: {}".format(len(df_diff.index), np.sum(df_diff["diff"].values >0), np.sum(df_diff["diff"].values>2))
                # fig.text(.5, .05, txt, ha='center')
                # plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR,
                #                          "diff_dist_{}.png".format(algo)))
                # plt.clf()

                total_n_terms.append(len(df_diff.index))
                n_terms_filtered_in.append(np.sum(df_diff["diff"].values > 0))
                diff_values=np.append(diff_values, df_diff["diff"][df_diff["diff"]>0].values)
                n_terms_2_oom.append(np.sum(df_diff["diff"].values > 2))
                # summary.append([len(df_go.index), empirical_th])
                # print df_go[["GO name","pval"]]
                # print "percentile: {}".format(empirical_th)

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

        print summary



