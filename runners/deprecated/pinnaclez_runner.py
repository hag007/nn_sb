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
import subprocess

# import rpy2.robjects.numpy2ri  as numpy2ri
# numpy2ri.activate()

from rpy2.robjects import pandas2ri
pandas2ri.activate()


import constants

from utils.r_runner import run_rscript
from utils.server import get_parameters

import infra

import utils.go

def sif2hotnet2(network_name):
    network_file_name = os.path.join(constants.NETWORKS_DIR, "{}.sif".format(network_name))

    network_df = pd.read_csv(network_file_name, sep="\t")
    src = np.array(network_df["ID_interactor_A"])
    dst = np.array(network_df["ID_interactor_B"])

    vertices = list(set(np.append(src,dst)))

    lns = ["{} {}".format(i+1, cur) for i, cur in enumerate(vertices)]
    file(os.path.join(constants.CACHE_DIR, "hotnet2_vertices.txt"), "w+").write("\n".join(lns))

    inds = map(lambda x: (vertices.index(x[0])+1, vertices.index(x[1])+1), zip(src,dst))


    lns = ["{} {} {}".format(cur[0], cur[1], 1) for cur in inds]
    file(os.path.join(constants.CACHE_DIR, "hotnet2_edges.txt"), "w+").write("\n".join(lns))

    print subprocess.Popen("bash ../sh/scripts/prepare_hotnet2.sh.format", shell=True, stdout=subprocess.PIPE).stdout.read() # cwd=dir_path

def run_hotnet2(deg_file_name, network_file_name):
    script = file("scripts/bionet.r").read()
    return run_rscript(script=script, output_vars = ["module_genes", "bg_genes"], network_file_name=network_file_name, deg_file_name=deg_file_name)


if __name__ == "__main__":
    constants.update_dirs(DATASET_NAME_u="TNFa")
    params = get_parameters()
    if params != None:
        args, NETWORK_NAME, dataset_name = params

    print subprocess.Popen("bash ../sh/scripts/run_pinnaclez.sh.format", shell=True,
                           stdout=subprocess.PIPE).stdout.read()  # cwd=dir_path

    results = file(os.path.join(constants.OUTPUT_DIR,"pinnaclez_results.txt")).read().split()
    module_genes = list(set([x for x in results if x.startswith("ENSG")]))
    dip_network = pd.read_csv(os.path.join(constants.NETWORKS_DIR, "dip_out.sif"), sep="\t", index_col=False, header=None)
    bg_genes = set(dip_network.ix[:,0]).union(set(dip_network.ix[:,2]))
    exp_genes, _1, _2 = infra.separate_headers(infra.load_gene_expression_profile())
    # bg_genes = list(bg_genes.union(set(exp_genes)))
    bg_genes = list(bg_genes)
    file(os.path.join(constants.OUTPUT_DIR, "pinnaclez_bg_genes.txt"), "w+").write("\n".join(bg_genes))
    file(os.path.join(constants.OUTPUT_DIR,"pinnaclez_module_genes.txt"), "w+").write("\n".join(module_genes))

    utils.go.check_group_enrichment(module_genes, bg_genes)






