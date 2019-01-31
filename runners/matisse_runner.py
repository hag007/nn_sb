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
import json
from utils.ensembl2entrez import ensembl2entrez_convertor
from utils.ensembl2gene_symbol import e2g_convertor

# import rpy2.robjects.numpy2ri  as numpy2ri
# numpy2ri.activate()

from rpy2.robjects import pandas2ri
pandas2ri.activate()

import constants

import infra

import utils.go

import DEG_runner
from utils.scripts import format_script
from utils.network import get_network_genes
from utils.network import build_all_reports
import utils.server as server
from utils.network import output_modules

ALGO_NAME = "matisse"
ALGO_DIR = os.path.join(constants.ALGO_BASE_DIR, "expander")
NETWORK_NAME = "dip"


def init_specific_params(ge_file_name=os.path.join(constants.DATA_DIR, "ge.tsv"), network_file_name=os.path.join(constants.NETWORKS_DIR, NETWORK_NAME + ".sif")):
    h_rows, h_columns, values = infra.separate_headers(infra.load_gene_expression_profile_by_genes(gene_expression_file_name=ge_file_name))
    df_ge = pd.DataFrame(columns=h_columns, index=h_rows, data=values)
    df_ge_cond_col = df_ge.columns
    df_ge["gene ID"] = df_ge.index
    df_ge["GeneName"] = [e2g_convertor([cur])[0] if len(e2g_convertor([cur])) > 0 else np.NAN for cur in
                   df_ge.index]
    df_ge = df_ge[["gene ID", "GeneName"] + list(df_ge_cond_col)]
    df_ge= df_ge[~df_ge['gene ID'].duplicated(keep='first')]
    ge_file_name_mts =os.path.splitext(ge_file_name)[0]+"_mts.tsv"
    df_ge.to_csv(ge_file_name_mts, index=False, sep="\t")

    output_file_name = os.path.join(constants.OUTPUT_DIR, "matisse_output.txt")
    return ge_file_name_mts, network_file_name, output_file_name


def extract_modules_and_bg(bg_genes, results_file_name):
    results = file(results_file_name).readlines()
    modules = [[] for x in range(max([int(x.strip().split()[2]) for x in results[1:]])+1    )]
    for x in results:
        modules[int(x.strip().split()[2])].append(x.strip().split()[0])
    all_bg_genes = [bg_genes for x in modules]
    print "extracted {} modules".format(len(modules))
    return modules, all_bg_genes


def main(dataset_name=constants.DATASET_NAME, disease_name=None, expected_genes = None, score_method=constants.DEG_EDGER):
    global NETWORK_NAME
    constants.update_dirs(DATASET_NAME_u=dataset_name)
    network_file_name, score_file_name, score_method, bg_genes= server.init_common_params(NETWORK_NAME, score_method)

    ge_file_name, network_file_name, output_file_name = init_specific_params(ge_file_name=os.path.join(constants.DATA_DIR, "ge.tsv"), network_file_name=os.path.join(constants.NETWORKS_DIR, NETWORK_NAME + ".sif"))

    format_script(os.path.join(constants.SH_DIR, "run_{}.sh".format(ALGO_NAME)),ALGO_BASE_DIR=constants.ALGO_BASE_DIR, GE_FILE_NAME=ge_file_name, NETWORK_FILE_NAME=network_file_name, BETA=0.95, MINIMAL_MODULE_SIZE=4, MAXIMAL_MODULE_SIZE=1000, OUTPUT_FILE_NAME=output_file_name)

    subprocess.Popen("bash {}/run_{}.sh".format(constants.SH_DIR, ALGO_NAME), shell=True,
                     stdout=subprocess.PIPE, cwd=ALGO_DIR).stdout.read()

    modules, all_bg_genes= extract_modules_and_bg(bg_genes, output_file_name)

    output_base_dir = ""
    if constants.REPORTS:
        output_base_dir = build_all_reports(ALGO_NAME, modules, all_bg_genes, score_file_name, network_file_name, disease_name, expected_genes)

    output_file_name=os.path.join(constants.OUTPUT_DIR,
                 "{}_client_output.txt".format(ALGO_NAME))
    output_modules(output_file_name, modules, score_file_name, output_base_dir)


# if __name__ == "__main__":
    # main(dataset_name="TNFa_2")
    # main(dataset_name="MCF7_2")
    # main(dataset_name="SOC")
    # main(dataset_name="HC12")
    # main(dataset_name="IEM", score_method=constants.DEG_T)
    # main(dataset_name="IEN", score_method=constants.DEG_T)
    # main(dataset_name="IES", score_method=constants.DEG_T)









