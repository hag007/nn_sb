#!/usr/bin/env python
"""Calculate differentially expressed genes using EdgeR from bioconductor.
http://bioconductor.org/packages/2.5/bioc/html/edgeR.html
Usage:
    count_diffexp.py <count_file>
"""
import sys

from scipy.optimize._remove_redundancy import bg_update_dense

sys.path.insert(0, '../')

import os
import numpy as np
from numpy import log10
import pandas as pd
import subprocess
import simplejson as json

# import rpy2.robjects.numpy2ri  as numpy2ri
# numpy2ri.activate()

from rpy2.robjects import pandas2ri
pandas2ri.activate()


import constants

from utils.r_runner import run_rscript
from utils.server import get_parameters
from utils.scripts import format_script
from utils.network import build_all_reports
import utils.server as server
import DEG_runner

import infra

import utils.go

from utils.network import output_modules

ALGO_NAME = "hotnet2"
ALGO_DIR = os.path.join(constants.ALGO_BASE_DIR, ALGO_NAME)

def sif2hotnet2(network_file_name, script_file_name):

    network_df = pd.read_csv(network_file_name, sep="\t")
    src = np.array(network_df["ID_interactor_A"])
    dst = np.array(network_df["ID_interactor_B"])

    vertices = list(set(np.append(src,dst)))

    lns = ["{} {}".format(i+1, cur) for i, cur in enumerate(vertices)]
    file(os.path.join(constants.CACHE_DIR, "hotnet2_vertices.txt"), "w+").write("\n".join(lns))

    inds = map(lambda x: (vertices.index(x[0])+1, vertices.index(x[1])+1), zip(src,dst))


    lns = ["{} {} {}".format(cur[0], cur[1], 1) for cur in inds]
    file(os.path.join(constants.CACHE_DIR, "hotnet2_edges.txt"), "w+").write("\n".join(lns))

    print subprocess.Popen("bash {}".format(script_file_name) , shell=True, stdout=subprocess.PIPE).stdout.read() # cwd=dir_path


def init_specific_params(score_file_name, method=constants.DEG_EDGER, network_file_name=os.path.join(constants.NETWORKS_DIR, "dip.sif")):

    script_file_name=format_script(os.path.join(constants.SH_DIR, "prepare_hotnet2.sh"), ALGO_DIR=ALGO_DIR,
                  CACHE_DIR=constants.CACHE_DIR, cwd=ALGO_DIR)

    heat_file_name = os.path.join(constants.CACHE_DIR, "heatfile.txt")


    deg = infra.load_gene_expression_profile_by_genes(gene_expression_path=score_file_name)
    h_rows, h_cols, deg_data = infra.separate_headers(deg)
    ind = np.where(h_cols=="qval")[0][0]

    lns = []
    if method==constants.PREDEFINED_SCORE and constants.IS_PVAL_SCORES:
        for i, cur in enumerate(deg_data):
            lns.append(" ".join([str(h_rows[i]), str(-log10(cur[ind]))]))
    else:
        for i, cur in enumerate(deg_data):
            lns.append(" ".join([str(h_rows[i]), str(cur[ind])]))

    file(heat_file_name,"w+").write("\n".join(lns))

    sif2hotnet2(network_file_name, script_file_name)
    os.remove(script_file_name)
    # file(os.path.join(constants.OUTPUT_DIR, "hotnet2_bg_genes.txt"), "w+").write("\n".join(bg_genes))
    return heat_file_name, network_file_name


def extract_modules_and_bg(bg_genes):
    results = json.load(file(os.path.join(constants.OUTPUT_DIR, "results", "consensus", "subnetworks.json")))
    modules = [x["core"] for x in results["consensus"] if len(x["core"]) > 3 ]
    all_bg_genes = [bg_genes for x in modules]
    print "extracted {} modules".format(len(modules))
    return modules, all_bg_genes


def main(dataset_name=constants.DATASET_NAME, disease_name=None, expected_genes = None, score_method=constants.DEG_EDGER, network_file_name="dip.sif"):

    constants.update_dirs(DATASET_NAME_u=dataset_name)
    network_file_name, score_file_name, score_method, bg_genes = server.init_common_params(network_file_name,score_method)

    heat_file_name, network_file_name = init_specific_params(score_file_name, score_method, network_file_name)

    script_file_name=format_script(os.path.join(constants.SH_DIR, "run_{}.sh".format(ALGO_NAME)), ALGO_DIR=ALGO_DIR,
                  CACHE_DIR=constants.CACHE_DIR, OUTPUT_DIR=constants.OUTPUT_DIR, NETWORK_NAME=os.path.splitext(os.path.basename(network_file_name))[0])
    print subprocess.Popen("bash {}".format(script_file_name), shell=True,
                           stdout=subprocess.PIPE).stdout.read()  # cwd=dir_path
    os.remove(script_file_name)
    modules, all_bg_genes = extract_modules_and_bg(bg_genes)
    output_base_dir = ""
    if constants.REPORTS:
        output_base_dir = build_all_reports(ALGO_NAME, dataset_name, modules, all_bg_genes, score_file_name, network_file_name, disease_name, expected_genes)

    output_file_name = os.path.join(constants.OUTPUT_DIR,
                                    "{}_client_output.txt".format(ALGO_NAME))
    output_modules(output_file_name, modules, score_file_name,output_base_dir)


if __name__ == "__main__":
    constants.update_dirs(DATASET_NAME_u="TNFa_2")
    main()







