import sys
sys.path.insert(0, '../')

import argparse

import pandas as pd
import numpy as np
import shutil

from pandas.errors import EmptyDataError

import constants
import os

import utils.goids2gonames as goids2gonames

import multiprocessing
from utils.daemon_multiprocessing import func_star
from utils.randomize_data import get_permutation_name

from runners.datasets_multithread_runner import run_dataset

def calc_dist(algos, datasets, shared_list=None, is_max=True):
    try:

        for cur_algo in algos:
            algos_filter = cur_algo

            df_go_pvals = pd.DataFrame()
            df_go_pvals.index.name="GO id"
            for cur_ds in datasets:
                print "fetch permutation {}".format(cur_ds)
                n_modules=len(pd.read_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, cur_ds, cur_algo, "modules_summary.tsv"), sep='\t').index)
                go_results = [os.path.join(constants.OUTPUT_GLOBAL_DIR, cur_ds, cur_algo, cur_module) for cur_algo in
                              os.listdir(os.path.join(constants.OUTPUT_GLOBAL_DIR, cur_ds))
                              if os.path.isdir(
                        os.path.join(constants.OUTPUT_GLOBAL_DIR, cur_ds, cur_algo)) and cur_algo in algos_filter for
                              cur_module in os.listdir(os.path.join(constants.OUTPUT_GLOBAL_DIR, cur_ds, cur_algo)) if
                              "separated_modules" in cur_module and int(cur_module.split("_")[1]) < n_modules]

                for cur in go_results:
                    try:
                        df_go_pvals = pd.concat((df_go_pvals, pd.read_csv(cur, sep='\t').set_index("GO id")['pval']), axis=1)

                        if is_max:
                            df_go_pvals[df_go_pvals.isna()] = 1
                            df_go_pvals=df_go_pvals.min(axis=1).to_frame()

                    except EmptyDataError,e:
                        print e
                        pass
                if len(go_results)==0:
                    df_go_pvals=pd.DataFrame(data=np.array([[1]]),index=["GO:0008150"])
 
            if not is_max:
                df_go_pvals[df_go_pvals.isna()] = 1

            if shared_list is not None:
                shared_list.append(df_go_pvals)
                print "done aggregate {} permutations".format(len(shared_list)) 
            return df_go_pvals
    except Exception, e:
       print Exception, e 
       pass


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
    parser.add_argument('--pf', dest='pf', help="parallelization factor", default=20)
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
    recalc_true_modules=args.recalc_true_modules.lower()=="true"
    pval_measurement=args.pval_measurement
    pf=args.pf


    for dataset in datasets:


        score_method = constants.PREDEFINED_SCORE
        if prefix == "GE":
            score_method = constants.DEG_EDGER
            if dataset.startswith("IE"):
                score_method = constants.DEG_T


        for algo in algos:

            manager = multiprocessing.Manager()
            pvals = manager.list()
            params = []
            p = multiprocessing.Pool(int(pf))

            total_n_terms = []
            n_terms_filtered_in = []
            n_terms_2_oom = []
            diff_values = np.array([0])
            df_all_terms = pd.DataFrame()

            
            params=[[calc_dist, [[algo], [get_permutation_name(prefix, dataset, algo,  cur)], pvals]] for cur in range(n_start,n_end)]
            print "test"
            p.map(func_star, params)
            pvals=list(pvals)
            df_all_terms = pd.concat(pvals, axis=1)
            df_all_terms=df_all_terms.fillna(1)
            print "total # permutations: {}/{}".format(len(pvals), n_end-n_start)
            
            if recalc_true_modules:
                if os.path.exists(os.path.join(constants.OUTPUT_GLOBAL_DIR,  "{}_{}".format(prefix, dataset), algo)):
                    shutil.rmtree(os.path.join(constants.OUTPUT_GLOBAL_DIR, "{}_{}".format(prefix, dataset), algo))
                run_dataset("{}_{}".format(prefix, dataset), score_method=score_method,
                            algos=[algo], network_file_name=network_file_name)


            ## empirical pvalues
            df_real_agg_pval = calc_dist([algo], ["{}_{}".format(prefix, dataset)])
            df_real_agg_pval=df_real_agg_pval.apply(lambda x: -np.log10(x))
            df_real_max_pval=df_real_agg_pval.max(axis=1).to_frame()
            print "total # real terms: {}".format(len(df_real_max_pval.index))
            df_counter=0
            df_results=df_all_terms.apply(lambda row : str(list(-np.log10(row.values.astype(np.float)))), axis=1).to_frame()
            df_results.columns = ['dist_n_samples']
            missing_indices=set(df_real_max_pval.index).difference(df_results.index)
            df_results.loc[missing_indices, "dist_n_samples"]=str([0])
            df_results['hg_pval']= df_real_agg_pval.max(axis=1)
            df_results["GO name"] = pd.Series(goids2gonames.get_go_names(list(df_real_max_pval.index)),
                                           index=df_real_max_pval.index)
            df_results.loc[df_results["hg_pval"].isna().values, "hg_pval"] = 0

            df_results.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "emp_diff_{}_{}.tsv".format(dataset, algo)),
                                        sep='\t', index_label="GO id")

            if args.max_dist:
                 df_results.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "emp_fdr","MAX","emp_diff_{}_{}.tsv".format(dataset, algo)),  sep='\t', index_label="GO id")
     
            print "permutation shape: {}".format(df_all_terms)
   
