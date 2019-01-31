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
        emp_pval=np.array([1])
        for algo in algos:
            corrupted_datsets=[]
            for cur_real_from_rand in range(0, 100):
                df_all_terms = pd.DataFrame()
                cur_real_ds= "{}_random_{}_{}".format(prefix, dataset, cur_real_from_rand)
                print "\ncur outer iteration: {}".format(cur_real_from_rand)

                # if not os.path.exists(os.path.join(constants.OUTPUT_GLOBAL_DIR, cur_real_ds, algo, "modules_summary.tsv")):
                #     print "outer {} is corrupted. continue...".format(cur_real_ds)
                #     corrupted_datsets.append(cur_real_ds)
                #     continue
                pval = np.array([])
                random_ds = "{}_random_{}".format(prefix, dataset)
                for cur in range(0,1000):

            
                    print "\ncur inner iteration: {}".format(cur)
                    random_ds = "{}_random_{}_{}".format(prefix, dataset, cur)
                    # random_ds=create_random_ds(prefix, "{}_{}".format(prefix, dataset),cur)
                    # permuted_network_file_name="dip" # _perm

                    if cur==cur_real_from_rand:
                        print "same datasets: {}=={} continue...".format(cur,cur_real_from_rand)
                        continue
                    
                    # if not os.path.exists(
                    #         os.path.join(constants.OUTPUT_GLOBAL_DIR, random_ds, algo, "modules_summary.tsv")):
                    #     print "inner {} is corrupted. continue...".format(cur)
                    #     corrupted_datsets.append(cur_real_ds)
                    #     continue




                    # if cur==0:
                    #     permuted_network_file_name=create_permuted_network(network_file_name=network_file_name)
                    # run_dataset(random_ds, score_method=score_method,
                    #             algos=[algo], network_file_name=permuted_network_file_name)
                    cur_pval, df_terms, df_pval_terms = calc_dist([algo], [random_ds.format(prefix, dataset)])
                    pval = np.append(pval, cur_pval)
                    df_all_terms=pd.concat((df_all_terms, df_pval_terms.max(axis=1).to_frame()), axis=1)
        #
                df_all_terms = df_all_terms.apply(lambda x: -np.log10(x))
                df_all_terms[df_all_terms.isna()]=0
#               df_all_terms_th=df_all_terms.apply(lambda r: np.percentile(r.values, 100, interpolation='higher'), axis=1)
#               df_temp=df_all_terms.apply(lambda r: np.percentile(r.values, 100, interpolation='higher'), axis=1)
#               df_temp=df_temp.to_frame()
#               df_temp.columns=["pval"]
#               df_temp["GO name"] = pd.Series(goids2gonames.get_go_names(list(df_temp.index)),
#                                                  index=df_temp.index)

#               df_temp.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "emp_go_terms_th_values.tsv"), sep='\t', index_label="GO id")


                pval, df_go, df_agg_pval = calc_dist([algo], [cur_real_ds], is_plot=False, empirical_th=None)
                df_agg_pval=df_agg_pval.apply(lambda x: -np.log10(x))
                df_max_pval=df_agg_pval.max(axis=1).to_frame()
                df_pvals = pd.DataFrame()
                df_max_pval.columns = ['pval']
                # print "testme: {}".format(df_all_terms.shape)
                for index, row in df_all_terms.iterrows():
                    pos = np.size(row) - np.searchsorted(np.sort(row), df_max_pval.loc[index, :].iloc[0] if index in df_max_pval.index else 0, side='left')
                    df_pvals.loc[index, "emp_pval"] = pos / float(np.size(row))
                    # df_pvals.loc[index, "dist_n_samples"] = str(list(row.values))
                    # df_pvals.loc[index, "sample_pos"] = pos
                    # print "current pval: {}".format(df_pvals.loc[index, "emp_pval"])

                for index, row in df_max_pval.iterrows():
                    if index not in df_pvals.index:
                        df_pvals.loc[index, "emp_pval"] = 0
                        # df_pvals.loc[index, "dist_n_samples"] = 0
                        # df_pvals.loc[index, "sample_pos"] = 0

                emp_pval=np.append(emp_pval, df_pvals["emp_pval"].values)
                print "current pvals: {}".format(emp_pval)
                txt = "total # diff: {}. 90'th percentile ={}. 95'th percentile : {}, min: {}".format(np.size(emp_pval), np.percentile(emp_pval, 10, interpolation='lower'), np.percentile(emp_pval, 5, interpolation='lower'), np.min(emp_pval))
                print txt
