#!/usr/bin/env python
"""Calculate differentially expressed genes using EdgeR from bioconductor.
http://bioconductor.org/packages/2.5/bioc/html/edgeR.html
Usage:
    count_diffexp.py <count_file>
"""
import os
import sys
import csv
import collections

import pandas as pd

import numpy as np
import rpy2.robjects as robjects
import rpy2.robjects.numpy2ri  as numpy2ri
numpy2ri.activate()
from rpy2.robjects import r, pandas2ri
pandas2ri.activate()
# import pandas.rpy.common as com


import constants
import infra


def main(count_file):
    base, ext = os.path.splitext(count_file)
    outfile = "%s-diffs.csv" % (base)
    ge_raw = infra.load_gene_expression_profile_by_genes()

    genes, conditions, data = infra.separate_headers(ge_raw)

    group = [1, 1, 1, 2, 2, 2]
    conditions=np.array(conditions)
    group=np.array(group)
    data = pd.DataFrame(data, index=genes, columns=conditions, dtype=np.int, )
    probs = run_rscript(data=data, genes=genes, conditions=conditions, group=group)


def write_outfile(outfile, genes, conditions, work_counts, probs):
    with open(outfile, "w") as out_handle:
        writer = csv.writer(out_handle)
        writer.writerow(["Region"] +
                        ["%s count" % c for c in conditions] + ["edgeR p-value"])
        out_info = []
        for i, gene in enumerate(genes):
            counts = [int(work_counts[c][gene]) for c in conditions]
            out_info.append((probs[i], [gene] + counts))
        out_info.sort()
        [writer.writerow(start + [prob]) for prob, start in out_info]


def run_rscript(script, **kwargs):
    """Call edgeR in R and organize the resulting differential expressed genes."""

    for k,v in kwargs.iteritems():
        robjects.globalenv[k] = v



    robjects.r(script)

    edgeR_results = robjects.r['results']

  #  edgeR_results = dict(zip(edgeR_results.names, list(edgeR_results)))
  #  edgeR_data = pandas2ri.ri2py(edgeR_results['table'])

    robjects.r('''
        
        
    ''')

    DESeq_results = robjects.r['DESeq_results']

    x=1






def get_conditions_and_genes(work_counts):
    conditions = work_counts.keys()
    conditions.sort()
    all_genes = []
    for c in conditions:
        all_genes.extend(work_counts[c].keys())
    all_genes = list(set(all_genes))
    all_genes.sort()
    sizes = [sum(work_counts[c].values()) for c in conditions]
    # all_genes.remove("Total")
    return conditions, all_genes, sizes



if __name__ == "__main__":
    constants.DATASET_DIR = os.path.join(constants.DATASETS_DIR, "TNFa_2")
    main(os.path.join(constants.DATASET_DIR, "ge.tsv"))
    print