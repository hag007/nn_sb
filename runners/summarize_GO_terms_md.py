import utils.add_GO_terms_metadata
import pandas as pd
import numpy as np
import os
import constants

datasets=["SHEZH_1"] # , "IEM" , "IES", "ROR_2", "SHEZH_1", "SHEZH_2", "ERS_1", "ERS_2"] # "SOC"
algos=["bionet", "jactivemodules_greedy", "jactivemodules_sa"] # , "bionet" # "hotnet2"


n_terms = pd.DataFrame(index=algos, columns=datasets)
hg_cutoffs = pd.DataFrame(index=algos, columns=datasets)
emp_cutoffs = pd.DataFrame(index=algos, columns=datasets)
for cur_ds in datasets:
    for cur_alg in algos:
        results=utils.add_GO_terms_metadata.main(cur_ds, cur_alg)
        if results is None:
            continue
        n_filtered_terms, n_corrected_terms, hg_cutoff, emp_cutoff=results
        hg_cutoffs.loc[cur_alg ,cur_ds] = hg_cutoff
        emp_cutoffs.loc[cur_alg, cur_ds] = emp_cutoff
        n_terms.loc[cur_alg, cur_ds] = "{}/{} ({})".format(n_corrected_terms, n_filtered_terms, n_corrected_terms/float(n_filtered_terms))


    hg_cutoffs.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "emp_fdr", "hg_cutoffs_summary.tsv"), sep='\t')
    emp_cutoffs.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "emp_fdr","emp_cutoffs_summary.tsv"), sep='\t')
    n_terms.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "emp_fdr","n_terms_summary.tsv"), sep='\t')





