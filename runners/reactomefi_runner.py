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
import utils.server as server
import infra

import utils.go

import DEG_runner

from utils.ensembl2entrez import ensembl2entrez_convertor
from utils.scripts import format_script
from utils.server import get_parameters
from utils.server import get_score
from utils.network import output_modules

ALGO_NAME = "reactomefi"
ALGO_DIR = os.path.join(constants.ALGO_BASE_DIR, ALGO_NAME)
NETWORK_NAME = "dip"


def extract_modules_and_bg(bg_genes):
    results = file(os.path.join(constants.OUTPUT_DIR, "reactomefi_modules.txt")).readlines()
    modules = [x.split("\t")[0].split() for x in results]
    modules = [cur for cur in modules if len(modules) > 2]
    all_bg_genes = [bg_genes for x in modules]
    print "extracted {} modules".format(len(modules))
    return modules, all_bg_genes


def init_specific_params(NETWORK_NAME):
    network_file_name = os.path.join(constants.NETWORKS_DIR, "{}.sif".format(NETWORK_NAME))
    bg_genes = get_network_genes()
    return bg_genes, network_file_name


def main(dataset_name=constants.DATASET_NAME, disease_name=None, expected_genes = None, score_method=constants.DEG_EDGER):
    global NETWORK_NAME
    constants.update_dirs(DATASET_NAME_u=dataset_name)
    network_file_name, score_file_name, score_method, bg_genes = server.init_common_params(NETWORK_NAME, score_method)
    if score_method == constants.PREDEFINED_SCORE:
        raise Exception("Cannot run this algo on scor-based metrics. please provide gene expression file")

    bg_genes, network_file_name = init_specific_params(NETWORK_NAME)

    format_script(os.path.join(constants.SH_DIR, "run_{}.sh".format(ALGO_NAME)), BASE_FOLDER=constants.BASE_PROFILE,
                  DATASET_DIR=constants.DATASET_DIR, ALGO_DIR=ALGO_DIR, NETWORK_NAME=NETWORK_NAME)

    subprocess.Popen("bash {}/run_{}.sh".format(constants.SH_DIR, ALGO_NAME), shell=True,
                     stdout=subprocess.PIPE, cwd=ALGO_DIR).stdout.read()

    modules, all_bg_genes = extract_modules_and_bg(bg_genes)
    output_base_dir = ""
    if constants.REPORTS:
        output_base_dir = build_all_reports(ALGO_NAME, modules, all_bg_genes, score_file_name, network_file_name, disease_name, expected_genes)

    output_file_name = os.path.join(constants.OUTPUT_DIR,
                                    "{}_client_output.txt".format(ALGO_NAME))
    output_modules(output_file_name, modules, score_file_name, output_base_dir )

if __name__ == "__main__":
    constants.update_dirs(DATASET_NAME_u="SOC")
    main()
