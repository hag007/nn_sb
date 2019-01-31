import sys
sys.path.insert(0, '../')

import os
import subprocess
import constants

import pandas as pd
import argparse


current_path=os.path.dirname(os.path.realpath(__file__))


def execute_stage(py_script, params):
    try:
        df_status_report = pd.read_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "status_report.tsv"), sep="\t", index_col=0)
    except:
        print "no status report found. creating new one..."
        df_status_report=pd.DataFrame()

    df_status_report.loc[cur_alg, cur_ds] = py_script
    df_status_report.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "status_report.tsv"), sep="\t")

    params = " ".join(params)
    print "about to start script {} with params:\n{}".format(py_script, params)
    out=subprocess.Popen("{}/../python27/bin/python {} {}".format(constants.REPO_DIR, py_script, params), shell=True,
                           stdout=subprocess.PIPE, cwd=current_path)
    print out.stdout.read()
    out.wait()
    if out.returncode!=0:
        raise Exception, "Error in {}: expected return code 0 but got {}".format(py_script, out.returncode)

    return out.returncode



if __name__=="__main__":

    parser = argparse.ArgumentParser(description='args')
    parser.add_argument('--datasets', dest='datasets', default="SOC")
    parser.add_argument('--prefix', dest='prefix', default="GE")
    parser.add_argument('--algos', dest='algos', default="jactivemodules_greedy")
    parser.add_argument('--network', dest='network', default="dip.sif")
    parser.add_argument('--n_start', help="number of iterations (total n permutation is pf*(n_end-n_start))", dest='n_start', default=0)
    parser.add_argument('--n_end', help="number of iterations (total n permutation is pf*(n_end-n_start))", dest='n_end', default=1000)
    parser.add_argument('--pf', dest='pf', help="parallelization factor", default=20)
    parser.add_argument('--recalc_true_modules', dest='recalc_true_modules', default="false")
    parser.add_argument('--n_iteration', dest='n_iteration', default=100)
    parser.add_argument('--n_total_samples', help="n_total_samples", dest='n_total_samples', default=1000)
    parser.add_argument('--n_dist_samples', help="n_dist_samples", dest='n_dist_samples', default=200)
    parser.add_argument('--override_permutations', help="override_permutations", dest='override_permutations', default='false') 
    args = parser.parse_args()

    datasets=args.datasets.split(",")
    algos=args.algos.split(",")
    prefix = args.prefix
    network_file_name = args.network
    n_start=int(args.n_start)
    n_end=int(args.n_end)
    recalc_true_modules=args.recalc_true_modules.lower()
    override_permutations=args.override_permutations
    n_iteration = int(args.n_iteration)
    n_total_samples = int(args.n_total_samples)
    n_dist_samples = int(args.n_dist_samples)
    n_permutations = int(args.n_total_samples)
    pf=args.pf

    prefix_param = "--prefix {}".format(prefix)
    n_start_param = "--n_start {}".format(n_start)
    n_end_param = "--n_end {}".format(n_end)
    pf_param = "--pf {}".format(pf)
    recalc_true_modules_param = "--recalc_true_modules {}".format(recalc_true_modules)
    n_iteration_param = "--n_iteration {}".format(n_iteration)
    n_total_samples_param = "--n_total_samples {}".format(n_total_samples)
    n_dist_samples_param = "--n_dist_samples {}".format(n_dist_samples)
    override_permutations_param = "--override_permutations {}".format(override_permutations)
    n_permutations_param = "--n_permutations {}".format(n_permutations)
    
    for cur_ds in datasets:
        for cur_alg in algos:
            try:
               datasets_param="--datasets {}".format(cur_ds)
               algos_param="--algos {}".format(cur_alg)


               # py_script = "generate_bg_dist.py"
               # execute_stage(py_script, [datasets_param, algos_param, prefix_param, n_start_param, n_end_param, pf_param, override_permutations_param])

               # py_script = "aggregate_bg_dist.py"
               # execute_stage(py_script, [datasets_param, algos_param, prefix_param, n_start_param, n_end_param, pf_param, recalc_true_modules_param])

               # py_script = "add_terms_md.py"
               # execute_stage(py_script, [datasets_param, algos_param, prefix_param, n_permutations_param])

               py_script = "fdr_consistent_terms.py"
               execute_stage(py_script, [datasets_param, algos_param, prefix_param, pf_param, n_iteration_param, n_total_samples_param, n_dist_samples_param, n_iteration_param])

               py_script = "report_result.py"
               execute_stage(py_script, [datasets_param, algos_param, prefix_param, pf_param])
    
               try:
                      df_status_report = pd.read_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "status_report.tsv"), sep="\t", index_col=0)
               except:
                      print "no status report found. creating new one..."
                      df_status_report=pd.DataFrame()
               df_status_report.loc[cur_alg, cur_ds] = "done!"
               df_status_report.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "status_report.tsv"), sep="\t") 
            except Exception, e:
               print "error in {}, {}: {}".format(cur_ds, cur_alg, e)
               try:
                  df_status_report = pd.read_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "status_report.tsv"), sep="\t", index_col=0)
               except:
                  print "no status report found. creating new one..."
                  df_status_report=pd.DataFrame()
               df_status_report.loc[cur_alg, cur_ds] = df_status_report.loc[cur_alg, cur_ds]+"\nerror! {}".format(e)
               df_status_report.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "status_report.tsv"), sep="\t")
