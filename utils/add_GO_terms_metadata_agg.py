import os

import numpy as np
import pandas as pd
from statsmodels.sandbox.stats.multicomp import fdrcorrection0

import constants
import go_hierarcies
from utils.ensembl2entrez import entrez2ensembl_convertor

dict_result, go2geneids, geneids2go, entrez2ensembl = go_hierarcies.build_hierarcy(
            roots=['GO:0008150'])
vertices=dict_result.values()[0]['vertices']


def mean_difference(row, dataset_data, classes_data):
    try:
        return dataset_data.loc[entrez2ensembl_convertor(get_all_genes_for_term(vertices, row["index"], row["index"],
                                                                             True)), classes_data == 2].dropna().values.mean() - \
        dataset_data.loc[entrez2ensembl_convertor(get_all_genes_for_term(vertices, row["index"], row["index"],
                                                                             True)), classes_data == 1].dropna().values.mean()
    except:
        print "no gene were found for {}, {} (pval={})".format(row["index"], row["GO name"],
                                                                        row["hg_pval"])

def calc_empirical_pval(row, n_permutation):

    pval = np.array([float(x) for x in row["dist_n_samples"][1:-1].split(", ")])

    if len(pval) < n_permutation:
        raise ValueError, "too few samples: {} (expected at least {})".format(len(pval), n_permutation)

    else:
        pval = pval[:n_permutation]
        pos = np.size(pval) - np.searchsorted(np.sort(pval), row["hg_pval"], side='left')
        emp_pval = pos / float(np.size(pval))

    return emp_pval


def get_all_genes_for_term(vertices, cur_root, term, in_subtree):

    in_subtree= in_subtree or term==cur_root
    all_genes = set()
    if in_subtree:
        all_genes.update(go2geneids[cur_root])
    for cur_child in vertices[cur_root]["obj"].children:
        all_genes.update(get_all_genes_for_term(vertices, cur_child.id, term, in_subtree))
    return all_genes


def main(dataset="SOC", algo="jactivemodules_sa", n_permutations=300, csv_file_name=os.path.join(constants.OUTPUT_GLOBAL_DIR, "emp_fdr", "MAX/emp_diff_{dataset}_{algo}.tsv")):

    dataset_data=pd.read_csv(os.path.join(constants.DATASETS_DIR, "GE_{}".format(dataset),"data", "ge.tsv"), sep='\t', index_col=0)
    classes_data=np.array(file(os.path.join(constants.DATASETS_DIR, "GE_{}".format(dataset), "data", "classes.tsv")).readlines()[0].strip().split("\t")).astype(np.int)

    csv_file_name=csv_file_name.format(dataset=dataset, algo=algo)
    df=None
    try:
        df=pd.read_csv(csv_file_name,sep='\t',index_col=0)
    except:
        return None
    df=df.dropna()

    n_genes=[len(get_all_genes_for_term(vertices, cur_go_id, cur_go_id, cur_go_id==cur_go_id)) for i, cur_go_id in enumerate(df.index.values)]
    depth=[dict_result.values()[0]['vertices'][cur_go_id]['D'] for i, cur_go_id in enumerate(df.index.values)]
    df["n_genes"]=pd.Series(n_genes, index=df.index)
    df["depth"]=pd.Series(depth, index=df.index)
    df=df.rename(columns={"filtered_pval" : "hg_pval"})

    n_genes_pvals=df.loc[np.logical_and.reduce([df["n_genes"].values > 5, df["n_genes"].values < 500]), "hg_pval"].values

    print "total n_genes with pval:{}/{}".format(np.size(n_genes_pvals), 7435)
    n_genes_pvals=np.append(n_genes_pvals,np.zeros(7435-np.size(n_genes_pvals)))
    n_genes_pvals = [10**(-x) for x in n_genes_pvals]
    fdr_results = fdrcorrection0(n_genes_pvals, alpha=0.05, method='indep', is_sorted=False)
    true_counter = len([cur for cur in fdr_results[0] if cur == True])
    HG_CUTOFF=-np.log10(np.sort(n_genes_pvals))[true_counter - 1] if true_counter > 0 else 0
    print "cutoff: {}".format(HG_CUTOFF)

    df_filtered_in=df.loc[np.logical_and.reduce([df["n_genes"].values > 5, df["n_genes"].values < 500, df["hg_pval"].values >= HG_CUTOFF]), :]
    df_filtered_out = df.loc[~np.logical_and.reduce([df["n_genes"].values > 5, df["n_genes"].values < 500, df["hg_pval"].values >= HG_CUTOFF]), :]

    df_filtered_in["index"] = df_filtered_in.index.values
    df_filtered_in["emp_pval"] = df_filtered_in.apply(lambda row: calc_empirical_pval(row, n_permutations), axis=1)
    df_filtered_in["mean_difference"] = df_filtered_in.apply(lambda x: mean_difference(x, dataset_data, classes_data), axis=1)

    pvals_corrected = df_filtered_in["emp_pval"].values
    fdr_results = fdrcorrection0(pvals_corrected, alpha=0.05, method='indep', is_sorted=False)
    true_counter = len([cur for cur in fdr_results[0] if cur == True])
    emp_cutoff = np.sort(np.sort(pvals_corrected))[true_counter - 1] if true_counter > 0 else 0
    print "emp true hypothesis: {} (emp cutoff: {}, n={})".format(true_counter, emp_cutoff, len(fdr_results[0]))

    df_filtered_in["passed_fdr"]=df_filtered_in["emp_pval"].apply(lambda x: x<=emp_cutoff)

    df_filtered_in["emp_rank"]=df_filtered_in["emp_pval"].rank(ascending=1)
    df_filtered_in["hg_rank"]=df_filtered_in["hg_pval"].rank(ascending=0)

    df_filtered_in=df_filtered_in.sort_values(by=["emp_rank", "hg_rank"])

    df_all=pd.concat((df_filtered_in, df_filtered_out), axis=0)
    df_all.loc[df_all["hg_pval"].values > 0, :][["GO name", "hg_pval", "emp_pval", "hg_rank", "emp_rank",  "n_genes", "depth", "mean_difference", "passed_fdr"]].to_csv(csv_file_name[:-4]+"_md.tsv",sep='\t')
    return len(df_filtered_in.index), true_counter, HG_CUTOFF, emp_cutoff

if __name__ == "__main__":
    csv_file_name=os.path.join(constants.OUTPUT_GLOBAL_DIR,"emp_fdr", "MAX", "emp_diff_{dataset}_{algo}.tsv")
    main(dataset="ERS_1", algo="hotnet2", csv_file_name=csv_file_name)
    # csv_file_name=os.path.join(constants.OUTPUT_GLOBAL_DIR,"emp_fdr","MAX", "emp_diff_{dataset}_{algo}.tsv")
    # main(dataset="TNFa_2", algo="jactivemodules_sa", csv_file_name=csv_file_name)


