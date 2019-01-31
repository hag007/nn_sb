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
from utils.network import remove_subgraph_by_nodes

from utils.server import get_parameters

import infra
import DEG_runner
import utils.server as server

import utils.go

from utils.network import output_modules

ALGO_NAME = "keypathwayminer"
ALGO_DIR = os.path.join(constants.ALGO_BASE_DIR, ALGO_NAME)

NETWORK_NAME = "dip"

import random

def init_specific_params(score_file_name, score_method, omitted_genes, network_file_name, ts, dest_algo_dir):
    if os.path.exists(os.path.join(dest_algo_dir, "results")):
        shutil.rmtree(os.path.join(dest_algo_dir, "results"))

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

    new_network_file_name = remove_subgraph_by_nodes(omitted_genes, network_file_name, ts=ts)
    return score_file_name, new_network_file_name


def format_scripts(score_file_name, network_name="dip", STRATEGY="INES", algorithm="GREEDY", algo_dir=None, dataset_name=None):
    script_file_name=format_script(os.path.join(constants.SH_DIR, "run_{}.sh".format(ALGO_NAME)), BASE_FOLDER=constants.BASE_PROFILE, DATASET_DIR=os.path.join(constants.DATASETS_DIR,dataset_name), STRATEGY=STRATEGY, ALGORITHM=algorithm)
    format_script(os.path.join(algo_dir, "kpm.properties"), uid=False, base_folder=constants.BASE_PROFILE, network_name=network_name, algo_dir=algo_dir, algorithm=algorithm)
    format_script(os.path.join(algo_dir, "datasets_file.txt"), uid=False, base_folder=constants.BASE_PROFILE, dataset=dataset_name, score_file_name=score_file_name)

    return script_file_name


def extract_module_genes(bg_genes, STRATEGY, algorithm, dest_algo_dir):
    i = 1
    modules = []
    while os.path.exists(os.path.join(dest_algo_dir, "results", "Pathway-{}-NODES-.txt".format("%02d" % (i,)))) and i<2:
        results = file(
            os.path.join(dest_algo_dir, "results", "Pathway-{}-NODES-.txt".format("%02d" % (i,)))).readlines()
        results = map(lambda x: x.strip(), results)
        modules.append(results)
        i += 1

    module_genes = list(set([y for x in modules for y in x]))
    module_genes = list(set(module_genes))
    module_genes_file_name = os.path.join(constants.OUTPUT_DIR, "{}_{}_{}_module_genes.txt".format(ALGO_NAME, STRATEGY, algorithm))
    file(module_genes_file_name, "w+").write("\n".join(module_genes))

    return modules, [bg_genes for x in modules]


def main(dataset_name=constants.DATASET_NAME, disease_name=None, expected_genes = None, score_method=constants.DEG_EDGER, network_file_name="dip.sif"):
    constants.update_dirs(DATASET_NAME_u=dataset_name)
    network_file_name, score_file_name, score_method, bg_genes = server.init_common_params(network_file_name, score_method)
    strategy = "INES"
    algorithm = "GREEDY"
    omitted_genes = []
    modules = []
    all_bg_genes = []
    dest_algo_dir = "{}_{}".format(ALGO_DIR, random.random())
    shutil.copytree(ALGO_DIR, dest_algo_dir)
    empty_counter = 0
    for cur_i_module in range(40):
        binary_score_file_name, cur_network_file_name = init_specific_params(score_file_name, score_method, omitted_genes,
                                                                         network_file_name, str(random.random()), dest_algo_dir)

        script_file_name=format_scripts(score_file_name=binary_score_file_name, network_name=cur_network_file_name,
                       STRATEGY=strategy, algorithm=algorithm, algo_dir=dest_algo_dir, dataset_name=dataset_name)
        print subprocess.Popen("bash {}".format(script_file_name), shell=True,
                               stdout=subprocess.PIPE, cwd=dest_algo_dir).stdout.read()
        module, all_bg_gene = extract_module_genes(bg_genes, strategy, algorithm, dest_algo_dir)

        if len(module[0]) > 3:
            empty_counter=0
            modules.append(module[0])
            all_bg_genes.append(all_bg_gene[0])
        else:
            empty_counter+=1
        omitted_genes += list(module[0])
        os.remove(script_file_name)

        if empty_counter>3:
            print "got more that 3 smalle modules in row. continue..."
            break

    shutil.rmtree(dest_algo_dir)

    output_base_dir = ""
    if constants.REPORTS:
        output_base_dir = build_all_reports("{}_{}_{}".format(ALGO_NAME,strategy, algorithm), dataset_name, modules, all_bg_genes, score_file_name, network_file_name, disease_name, expected_genes)
    output_file_name = os.path.join(constants.OUTPUT_DIR,
                                    "{}_client_output.txt".format("{}_{}_{}".format(ALGO_NAME,strategy, algorithm)))
    output_modules(output_file_name, modules, score_file_name, output_base_dir )



if __name__ == "__main__":
    main(dataset_name="HC12")

