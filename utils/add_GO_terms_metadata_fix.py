import os

import numpy as np
import pandas as pd
from statsmodels.sandbox.stats.multicomp import fdrcorrection0

import constants
import go_hierarcies
from utils.ensembl2entrez import entrez2ensembl_convertor

# dict_result, go2geneids, geneids2go, entrez2ensembl = go_hierarcies.build_hierarcy(
#             roots=['GO:0008150'])
# vertices=dict_result.values()[0]['vertices']
#
# N_PERMUTATIONS=300

# def mean_difference(row, dataset_data, classes_data):
#     try:
#         return dataset_data.loc[entrez2ensembl_convertor(get_all_genes_for_term(vertices, row["index"], row["index"],
#                                                                              True)), classes_data == 2].dropna().values.mean() - \
#         dataset_data.loc[entrez2ensembl_convertor(get_all_genes_for_term(vertices, row["index"], row["index"],
#                                                                              True)), classes_data == 1].dropna().values.mean()
#     except:
#         print "no gene were found for {}, {} (pval={})".format(row["index"], row["GO name"],
#                                                                         row["hg_pval"])
#
# def calc_empirical_pval(row):
#
#     pval = np.array([float(x) for x in row["dist_n_samples"][1:-1].split(", ")])
#
#     if len(pval) < N_PERMUTATIONS:
#         raise ValueError
#
#     else:
#         pval = pval[:N_PERMUTATIONS]
#         pos = np.size(pval) - np.searchsorted(np.sort(pval), row["hg_pval"], side='left')
#         emp_pval = pos / float(np.size(pval))
#
#     return emp_pval
#
#
# def get_all_genes_for_term(vertices, cur_root, term, in_subtree):
#
#     in_subtree= in_subtree or term==cur_root
#     all_genes = set()
#     if in_subtree:
#         all_genes.update(go2geneids[cur_root])
#     for cur_child in vertices[cur_root]["obj"].children:
#         all_genes.update(get_all_genes_for_term(vertices, cur_child.id, term, in_subtree))
#     return all_genes


def main(dataset="SOC", algo="jactivemodules_sa", csv_file_name=os.path.join(constants.OUTPUT_GLOBAL_DIR, "emp_fdr","MAX/emp_diff_{dataset}_{algo}.tsv" )):

    # dataset_data=pd.read_csv(os.path.join(constants.DATASETS_DIR, "GE_{}".format(dataset),"data", "ge.tsv"), sep='\t', index_col=0)
    # classes_data=np.array(file(os.path.join(constants.DATASETS_DIR, "GE_{}".format(dataset), "data", "classes.tsv")).readlines()[0].strip().split("\t")).astype(np.int)
    csv_file_name=csv_file_name.format(dataset=dataset, algo=algo)
    print csv_file_name
    df=None
    try:
        df=pd.read_csv(csv_file_name, sep='\t',index_col=0)
    except:
        return None
    if np.max(np.array(df["dist_n_samples"].iloc[0][1:-1].split(", ")).astype(np.float)) <=1: 
       df["dist_n_samples"]=df["dist_n_samples"].apply(lambda x: str(list(-np.log10(np.array(x[1:-1].split(", ")).astype(np.float)))))
       df.to_csv(csv_file_name, sep='\t',index_label="GO id")
       print "transformation completed"
    else:
       print "already in -log10 transform"
    


if __name__ == "__main__":
    csv_file_name=os.path.join(constants.OUTPUT_GLOBAL_DIR,"emp_fdr", "MAX", "emp_diff_{dataset}_{algo}.tsv")
    main(dataset="ERS_1", algo="hotnet2", csv_file_name=csv_file_name)
