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
from utils.network import get_network_genes
from utils.network import build_all_reports
from utils.server import get_score
import infra

import utils.go

import DEG_runner

from utils.ensembl2entrez import ensembl2entrez_convertor
from utils.scripts import format_script
import utils.server as server
from utils.network import output_modules

ALGO_NAME = "netbox"
ALGO_DIR = os.path.join(constants.ALGO_BASE_DIR, ALGO_NAME)

import shutil
import random

def init_specific_params(score_file_name, dest_algo_dir):

    deg = infra.load_gene_expression_profile_by_genes(gene_expression_path=score_file_name)
    h_rows, h_cols, deg_data = infra.separate_headers(deg)

    ind = np.where(h_cols=="qval")[0][0]
    ordered_ind = np.argsort(deg_data[:,ind])
    deg_data=deg_data[ordered_ind,:]
    h_rows=h_rows[ordered_ind]
    last_q_index = np.where(deg_data[:,np.where(h_cols=="qval")[0][0]]>0.05)[0][0]
    ge_list_file_name = os.path.join(constants.OUTPUT_DIR, "ge_list.txt")
    file(ge_list_file_name, "w+").write("\n".join([x for x in h_rows[:last_q_index] if len(ensembl2entrez_convertor([x]))>0 ])) # ensembl2entrez_convertor([x])[0]

    conf_file = "conf.props"
    conf_file_name=format_script(os.path.join(dest_algo_dir, conf_file), pval_threshold=0.05, sp_threshold=2, gene_file=ge_list_file_name)
    return conf_file_name

def extract_modules_and_bg(bg_genes, dest_algo_dir):
    results = file(os.path.join(dest_algo_dir, "modules.txt")).readlines()
    modules = [[] for x in range(max([int(x.strip().split(" =")[1]) for x in results[1:]]) + 1)]
    for x in results[1:]:
        if int(x.strip().split(" =")[1]) != -1:
            modules[int(x.strip().split(" =")[1])].append(x.strip().split(" =")[0])
        else:
            modules.append([x.strip().split(" =")[0]])
    modules = filter(lambda x: len(x) > 3, modules)
    all_bg_genes = [bg_genes for x in modules]
    print "extracted {} modules".format(len(modules))
    return modules, all_bg_genes


def main(dataset_name=constants.DATASET_NAME, disease_name=None, expected_genes = None, score_method=constants.DEG_EDGER, network_file_name="dip.sif"):
    constants.update_dirs(DATASET_NAME_u=dataset_name)
    network_file_name, score_file_name, score_method, bg_genes = server.init_common_params(network_file_name, score_method)

    
    script_name = "run_{}.sh".format(ALGO_NAME)
    dest_algo_dir="{}_{}".format(ALGO_DIR,random.random())
    shutil.copytree(ALGO_DIR, dest_algo_dir)
    conf_file_name=init_specific_params(score_file_name, dest_algo_dir)
    script_file_name=format_script(os.path.join(constants.SH_DIR, script_name), BASE_FOLDER=constants.BASE_PROFILE,
                  DATASET_DIR=constants.DATASET_DIR, CONFIG_FILE_NAME=conf_file_name, NETBOX_DIR=dest_algo_dir)
    print subprocess.Popen("bash {}".format(script_file_name), shell=True,
                           stdout=subprocess.PIPE, cwd=dest_algo_dir).stdout.read()

    modules, all_bg_genes = extract_modules_and_bg(bg_genes, dest_algo_dir)
    os.remove(script_file_name)
    os.remove(conf_file_name)
    shutil.rmtree(dest_algo_dir)
    output_base_dir = ""
    if constants.REPORTS:
        output_base_dir = build_all_reports(ALGO_NAME, dataset_name, modules, all_bg_genes, score_file_name, network_file_name, disease_name, expected_genes)
    output_file_name = os.path.join(constants.OUTPUT_DIR,
                                    "{}_client_output.txt".format(ALGO_NAME))
    # output_modules(output_file_name, modules, score_file_name, output_base_dir )


if __name__ == "__main__":
    constants.update_dirs(DATASET_NAME_u="MCF7_2")
    main()







