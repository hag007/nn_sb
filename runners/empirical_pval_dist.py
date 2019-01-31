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
                except EmptyDataError:
                    pass

        df_go = df_go[df_go['qval'] < 0.05]
        if empirical_th:
            df_go = df_go[df_go['pval'].apply(lambda x:-np.log10(x)) > empirical_th]

        # df_go = df_go[~df_go.index.duplicated(keep='first')]
        pval = -np.log10(df_go["pval"].values)
        if np.size(pval) == 0:
            pval = np.array([0])

        if is_plot:
            plot_dist(pval, algos_filter=algos+datasets)

        return pval, df_go

if __name__ == "__main__":
    # prefix="GE"
    # datasets = ["{}_TNFa_2".format(prefix)]
    # algos = np.array(["jactivemodules_greedy"])
    #
    # calc_dist(algos, datasets)
    network_file_name="dip"
    prefix="GE"
    datasets = ["SOC"] # alzheimers, schizophrenia
    algos = ["jactivemodules_greedy"] # "bionet"

    score_method = constants.PREDEFINED_SCORE
    if prefix=="GE":
        score_method=constants.DEG_EDGER

    summary = []
    for dataset in datasets:
        for algo in algos:
            pval = np.array([])
            random_ds = "{}_random_{}".format(prefix, dataset)
            for cur in range(1):
                random_ds=create_random_ds(prefix, "{}_{}".format(prefix, dataset))
                permuted_network_file_name="dip" # _perm
                # if cur==0:
                #     permuted_network_file_name=create_permuted_network(network_file_name=network_file_name)
                # run_dataset(random_ds, score_method=score_method,
                #             algos=[algo], network_file_name=permuted_network_file_name)
                pval = np.append(pval, calc_dist([algo], [random_ds])[0])
            empirical_th=np.percentile(np.sort(pval), 40)

            # run_dataset("{}_{}".format(prefix, dataset), score_method=score_method,
            #             algos=[algo], network_file_name=network_file_name)
            plot_dist(pval, [algo, random_ds])

            pval, df_go = calc_dist([algo], ["{}_{}".format(prefix, dataset)], is_plot=True, empirical_th=empirical_th)
            summary.append([len(df_go.index), empirical_th])
            print df_go[["GO name","pval"]]
            print "percentile: {}".format(empirical_th)

    print summary


