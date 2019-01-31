import sys
sys.path.insert(0, '../')
import utils.add_GO_terms_metadata_agg
import pandas as pd
import numpy as np
import os
import constants
import argparse




if __name__=="__main__":


    parser = argparse.ArgumentParser(description='args')
    parser.add_argument('--datasets', dest='datasets', default="SOC")
    parser.add_argument('--prefix', dest='prefix', default="GE")
    parser.add_argument('--algos', dest='algos', default="jactivemodules_greedy")
    parser.add_argument('--n_permutations', dest='n_permutations', default=1000)

    args = parser.parse_args()

    datasets=args.datasets.split(",")
    algos=args.algos.split(",")
    prefix = args.prefix
    n_permutations=int(args.n_permutations)

    n_terms = pd.DataFrame(index=algos, columns=datasets)
    hg_cutoffs = pd.DataFrame(index=algos, columns=datasets)
    emp_cutoffs = pd.DataFrame(index=algos, columns=datasets)
    for cur_ds in datasets:
        for cur_alg in algos:
            print "{}_{}".format(cur_ds, cur_alg)
            results=utils.add_GO_terms_metadata_agg.main(cur_ds, cur_alg, n_permutations)
            if results is None:
                continue
            n_filtered_terms, n_corrected_terms, hg_cutoff, emp_cutoff=results
            hg_cutoffs.loc[cur_alg ,cur_ds] = hg_cutoff
            emp_cutoffs.loc[cur_alg, cur_ds] = emp_cutoff
            n_terms.loc[cur_alg, cur_ds] = "{}/{} ({})".format(n_corrected_terms, n_filtered_terms, n_corrected_terms/float(n_filtered_terms))


        hg_cutoffs.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "emp_fdr", "hg_cutoffs_summary.tsv"), sep='\t')
        emp_cutoffs.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "emp_fdr","emp_cutoffs_summary.tsv"), sep='\t')
        n_terms.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "emp_fdr","n_terms_summary.tsv"), sep='\t')






