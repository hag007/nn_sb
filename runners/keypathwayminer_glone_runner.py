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

import shutil
import constants

from utils.scripts import format_script
from utils.network import get_network_genes
from utils.network import build_all_reports

from utils.server import get_parameters

import infra
import DEG_runner
import utils.server as server

import utils.go

from utils.network import output_modules

ALGO_NAME = "keypathwayminer"
ALGO_DIR = os.path.join(constants.ALGO_BASE_DIR, ALGO_NAME)

NETWORK_NAME = "dip"

def init_common_params(score_file_name, score_method):
    if os.path.exists(os.path.join(ALGO_DIR, "results")):
        shutil.rmtree(os.path.join(ALGO_DIR, "results"))

    if score_method != constants.PREDEFINED_SCORE:
        deg = infra.load_gene_expression_profile_by_genes(gene_expression_path=score_file_name)
        h_rows, h_cols, deg_data = infra.separate_headers(deg)
        ind = np.where(h_cols=="qval")[0][0]
        ordered_ind = np.argsort(deg_data[:,ind])
        deg_data=deg_data[ordered_ind,:]
        h_rows=h_rows[ordered_ind]
        sig_binary_col = deg_data[:,np.where(h_cols=="qval")[0][0]]<0.05
        sig_binary_output = np.c_[h_rows,  np.array(sig_binary_col, dtype=np.int)]
        score_file_name = os.path.join(constants.CACHE_DIR, "binary_score_{}.txt".format(score_method))
        file(score_file_name, "w+").write("\n".join(["\t".join(["id", "pval", "qval"])] + ["\t".join(list(x) + list([x[-1]])) for x in sig_binary_output]))

    return score_file_name




def format_scripts(algo_name, score_file_name, network_name="dip", STRATEGY="INES"):
    format_script(os.path.join(constants.SH_DIR, "run_{}.sh".format(ALGO_NAME)), BASE_FOLDER=constants.BASE_PROFILE, DATASET_DIR=constants.DATASET_DIR, STRATEGY=STRATEGY)
    format_script(os.path.join(constants.ALGO_BASE_DIR,algo_name, "kpm.properties"), base_folder=constants.BASE_PROFILE, network_name=network_name, algo_dir=ALGO_DIR)
    format_script(os.path.join(constants.ALGO_BASE_DIR, algo_name, "datasets_file.txt"), base_folder=constants.BASE_PROFILE, dataset=constants.DATASET_NAME, score_file_name=score_file_name)


def extract_module_genes(bg_genes, STRATEGY):
    i = 1
    modules = []
    while os.path.exists(os.path.join(ALGO_DIR, "results", "Pathway-{}-NODES-.txt".format("%02d" % (i,)))):
        results = file(
            os.path.join(ALGO_DIR, "results", "Pathway-{}-NODES-.txt".format("%02d" % (i,)))).readlines()
        results = map(lambda x: x.strip(), results)
        modules.append(results)
        i += 1

    module_genes = list(set([y for x in modules for y in x]))
    module_genes = list(set(module_genes))
    module_genes_file_name = os.path.join(constants.OUTPUT_DIR, "{}_{}_module_genes.txt".format(ALGO_NAME, STRATEGY))
    file(module_genes_file_name, "w+").write("\n".join(module_genes))

    return modules, [bg_genes for x in modules]


def main(dataset_name=constants.DATASET_NAME, disease_name=None, expected_genes = None):
    global NETWORK_NAME
    constants.update_dirs(DATASET_NAME_u=dataset_name)
    network_file_name, score_file_name, score_method, bg_genes = server.init_common_params(NETWORK_NAME)
    STRATEGY = "GLONE"
    binary_score_file_name = init_common_params(score_file_name, score_method)
    format_scripts(algo_name=ALGO_NAME, score_file_name=binary_score_file_name, network_name=NETWORK_NAME, STRATEGY=STRATEGY)
    print subprocess.Popen("bash {}/run_{}.sh".format(constants.SH_DIR, ALGO_NAME), shell=True,
                           stdout=subprocess.PIPE, cwd=ALGO_DIR).stdout.read()
    modules, all_bg_genes = extract_module_genes(bg_genes, STRATEGY)

    output_base_dir = ""
    if constants.REPORTS:
        output_base_dir = build_all_reports(ALGO_NAME + "_" + STRATEGY, modules, all_bg_genes, score_file_name, network_file_name, disease_name, expected_genes)

    output_file_name = os.path.join(constants.OUTPUT_DIR,
                                    "{}_client_output.txt".format(ALGO_NAME))
    output_modules(output_file_name, modules, score_file_name, output_base_dir )


if __name__ == "__main__":
    main()
