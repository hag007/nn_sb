#!/usr/bin/env python
"""Calculate differentially expressed genes using EdgeR from bioconductor.
http://bioconductor.org/packages/2.5/bioc/html/edgeR.html
Usage:
    count_diffexp.py <count_file>
"""

import sys
sys.path.insert(0, '../')

import os
import numpy as np
import pandas as pd

# import rpy2.robjects.numpy2ri  as numpy2ri
# numpy2ri.activate()
from rpy2.robjects import pandas2ri
pandas2ri.activate()

import constants
from utils.r_runner import run_rscript
from utils.server import get_parameters

import infra

import utils.go


def prepare_input(method=constants.DEG_EDGER, network_name="dip"):
    ge_file_name = os.path.join(constants.CACHE_DIR, "deg_{}.tsv".format(method).format(method).format(method))
    network_file_name = os.path.join(constants.NETWORKS_DIR, "{}.sif".format(network_name))

    network_df = pd.read_csv(network_file_name, sep="\t")
    src = np.array(network_df["ID_interactor_A"])
    dst = np.array(network_df["ID_interactor_B"])

    vertices = list(set(np.append(src,dst)))
    ppi_i = []

    deg = infra.load_gene_expression_profile_by_genes(gene_expression_path=ge_file_name)
    h_rows, h_cols, deg_data = infra.separate_headers(deg)
    bg_genes = list(set(vertices).intersection(set(h_rows)))

    for i, cur_r in network_df.iterrows():
        if (cur_r["ID_interactor_A"] in bg_genes and cur_r["ID_interactor_B"] in bg_genes):
            ppi_i.append("\t".join([str(bg_genes.index(cur_r["ID_interactor_A"])+1),str(bg_genes.index(cur_r["ID_interactor_B"])+1)]))



    avg_rd = np.average(deg_data[:,0:5])
    avg_p_deg = np.average(deg_data[:, 6])
    avg_q_deg = np.average(deg_data[:, 7])

    normalized_ge = []

    for cur_v in bg_genes:
        normalized_ge.append(deg_data[np.where(h_rows==cur_v)[0][0]])

    pd.DataFrame(normalized_ge, index=bg_genes, columns=h_cols).to_csv(sep="\t",path_or_buf=os.path.join(constants.OUTPUT_DIR, "cosine_ge.tsv"))

    bg_genes_file_name=os.path.join(constants.OUTPUT_DIR, "cosine_bg_genes.txt")
    file(os.path.join(constants.OUTPUT_DIR, bg_genes_file_name), "w+").write("\n".join(bg_genes))


    file(os.path.join(constants.OUTPUT_GLOBAL_DIR,"ppi_i.txt"),"w+").write("\n".join(ppi_i))

    return network_file_name, bg_genes, vertices, \
           os.path.join(constants.OUTPUT_GLOBAL_DIR, "ppi_i.txt"), os.path.join(constants.OUTPUT_DIR, "cosine_ge.tsv")

def run_cosine(ppi_i_file, ge_file):
    script = file("../r/scripts/cosine.r").read()
    return run_rscript(script=script, output_vars = ["max_subnet", "subnets", "adjs"], ppi_i_file=ppi_i_file, ge_file=ge_file)

if __name__ == "__main__":
    params = get_parameters()
    if params != None:
        args, NETWORK_NAME, dataset_name = params

    score_method=constants.DEG_EDGER
    network_file_name, bg_genes, vertices_ids, ppi_i_file, ge_file = prepare_input(method=score_method)
    results = run_cosine(ppi_i_file=ppi_i_file, ge_file=ge_file)

    for i, subnet in enumerate(results['subnets']):
        result = [vertices_ids[x-1] for x in subnet]
        module_genes = list(set(result))
        file(os.path.join(constants.OUTPUT_DIR,"cosine_module_genes_{}_r.txt".format(i)), "w+").write("\n".join(module_genes))
        utils.go.check_group_enrichment(module_genes, bg_genes)

    result = [vertices_ids[x-1] for x in results['max_subnet']]
    module_genes = list(set(result))
    file(os.path.join(constants.OUTPUT_DIR,"cosine_module_genes_max_r.txt"), "w+").write("\n".join(module_genes))
    utils.go.check_group_enrichment(module_genes, bg_genes)

    print results['adjs']






