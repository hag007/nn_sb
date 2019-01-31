#!/usr/bin/env python
"""Calculate differentially expressed genes using EdgeR from bioconductor.
http://bioconductor.org/packages/2.5/bioc/html/edgeR.html
Usage:
    count_diffexp.py <count_file>
"""
import os
import pandas as pd

import numpy as np
# import rpy2.robjects.numpy2ri  as numpy2ri
# numpy2ri.activate()

from rpy2.robjects import pandas2ri
pandas2ri.activate()


import constants
import infra
import shutil
import time
import random
import utils.add_t_test_to_ge as t_test
from utils.r_runner import run_rscript


def run_DEG(method, conditions, data, genes, group, dataset=constants.DATASET_NAME):
    if method == constants.DEG_T:
        return t_test.calc_ttest(dataset)
    else:
        script = file(os.path.join(constants.REPO_DIR,"r","scripts","{}.r".format(method))).read()
        print group.shape
        print data.shape
        return run_rscript(script=script, data=data, genes=genes, conditions=conditions, group=group)


def prepare_input(gene_expression_file_name="ge.tsv", groups = None):
    if groups is not None:
        groups = infra.load_classes()
    elif os.path.exists(os.path.join(constants.DATA_DIR, "classes.tsv")):
        groups = infra.load_classes()
    else:
        groups = [1,1,1,2,2,2]


    ge_raw = infra.load_gene_expression_profile_by_genes(gene_expression_file_name=gene_expression_file_name)
    genes, conditions, data = infra.separate_headers(ge_raw)
    conditions = np.array(conditions)
    groups = np.array(groups,dtype=np.int)
    data = pd.DataFrame(data, index=genes, columns=conditions, dtype=np.float)
    return conditions, data, genes, groups



def main(method=constants.DEG_EDGER):
    conditions, data, genes, group = prepare_input(gene_expression_file_name="ge.tsv")

    results = run_DEG(method, conditions, data, genes, group, constants.DATASET_NAME)
    if method==constants.DEG_EDGER:
        res = results["result"][["PValue", "FDR"]]
        res = res.rename(columns={"PValue": "pval", "FDR": "qval"})
    elif method==constants.DEG_DESEQ:
        res = results["result"][["pvalue", "padj"]]
        res = res.rename(columns={"pvalue": "pval", "padj": "qval"})
    elif method==constants.DEG_T:
        res = results["result"][["pval", "qval"]]

    combined_DEG = pd.concat([data, res], axis=1, sort=False)
    combined_DEG = combined_DEG.sort_values("pval".format(method))
    rand = str(random.random())
    combined_DEG.to_csv(os.path.join(constants.CACHE_DIR, "deg_{}.tsv.tmp_{}".format(method, rand)), sep="\t", index_label="id")
    print "prepare {}".format(os.path.join(constants.CACHE_DIR, "deg_{}.tsv.tmp_{}".format(method, rand)))
    while not os.path.exists(os.path.join(constants.CACHE_DIR, "deg_{}.tsv.tmp_{}".format(method, rand))):
        print "waiting for {}".format(os.path.join(constants.CACHE_DIR, "deg_{}.tsv.tmp_{}".format(method,rand)))
    shutil.move(os.path.join(constants.CACHE_DIR, "deg_{}.tsv.tmp_{}".format(method, rand)), os.path.join(constants.CACHE_DIR, "deg_{}.tsv".format(method)))

if __name__ == "__main__":
    main()

