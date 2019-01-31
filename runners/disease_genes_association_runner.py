#!/usr/bin/env python
"""Calculate differentially expressed genes using EdgeR from bioconductor.
http://bioconductor.org/packages/2.5/bioc/html/edgeR.html
Usage:
    count_diffexp.py <count_file>
"""

import sys

sys.path.insert(0, '../')

import os
import time
# import rpy2.robjects.numpy2ri  as numpy2ri
# numpy2ri.activate()

from rpy2.robjects import pandas2ri

pandas2ri.activate()

import constants

from runners.datasets_multithread_runner import run_dataset

from utils.server import get_parameters
from utils.network import get_network_genes
from utils.disease_genes import get_genes_by_disease
from utils.disease_genes import get_all_disease_names

import pandas as pd
import numpy as np
import random

import utils.aggregate_reports as aggregate_reports

NETWORK_NAME = "dip"
MAX_DISEASES = 20
TOTAL_TO_NW_RATIO = 0.6
ABS_NW_GENES = 20


def create_ge(disease_nw_genes, bg_genes, dataset_name):
    df = pd.DataFrame(columns=["id", "CONTROL", "CASE"])
    for cur in disease_nw_genes:
        std = 1
        mean_case_up = np.floor(random.random() * 3000)
        mean_case_down = np.floor(random.random() * 5)
        cur_gene = {'id': cur}
        for cur_rep in range(3):
            cur_gene["CASE{}".format(cur_rep)] = max(0, int(np.random.randn() * std + mean_case_up))
            cur_gene["CONTROL{}".format(cur_rep)] = max(0, int(np.random.randn() * std + mean_case_down))

        df = df.append(cur_gene, ignore_index=True)
    for cur in (set(bg_genes) - set(disease_nw_genes)):
        mean_control = int(round(200 * random.random()))
        std = 1
        cur_gene = {'id': cur}
        for cur_rep in range(3):
            cur_gene["CASE{}".format(cur_rep)] = max(0, int(np.random.randn() * std + mean_control))
            cur_gene["CONTROL{}".format(cur_rep)] = max(0, int(np.random.randn() * std + mean_control))
        df = df.append(cur_gene, ignore_index=True)

    df = df[["id", "CASE0", "CASE1", "CASE2", "CONTROL0", "CONTROL1", "CONTROL2"]]

    df.to_csv(os.path.join(os.path.join(constants.DATASETS_DIR, dataset_name, "data"), "ge.tsv"), sep="\t", index=False)


def create_ds_folders(dataset_name):
    os.makedirs(os.path.join(os.path.join(constants.DATASETS_DIR, dataset_name, "data")))
    os.makedirs(os.path.join(os.path.join(constants.DATASETS_DIR, dataset_name, "cache")))
    os.makedirs(os.path.join(os.path.join(constants.DATASETS_DIR, dataset_name, "output")))


def main():
    global NETWORK_NAME
    params = get_parameters()
    if params != None:
        args, NETWORK_NAME, dataset_name = params
    ts = str(time.time())
    bg_genes = get_network_genes()
    disease_counter = 0
    for cur_disease_name in get_all_disease_names():
        disease_genes = get_genes_by_disease(disease_name=cur_disease_name)
        sys.stdout.write("total disease genes: {}\n".format(len(disease_genes)))
        disease_nw_genes = set(disease_genes).intersection(bg_genes)
        sys.stdout.write("total disease genes in network: {}\n".format(len(disease_nw_genes)))
        if len(disease_nw_genes) == 0:
            print "no genes for disease {}... skipping".format(cur_disease_name)
            continue

        total_to_nw_genes_ratio = len(disease_nw_genes) / float(len(disease_genes))
        abs_nw_genes = len(disease_nw_genes)
        sys.stdout.write("ratio={}, absolute n={}\n".format(total_to_nw_genes_ratio, abs_nw_genes))
        if total_to_nw_genes_ratio < TOTAL_TO_NW_RATIO or abs_nw_genes < ABS_NW_GENES:
            sys.stdout.write("disease genes in network is too small\n")
        else:
            print "about to start analyze disease: {}".format(cur_disease_name)
            disease_counter += 1
            dataset_name = "_".join(["DISGENET", ts, str(disease_counter)])
            create_ds_folders(dataset_name)
            disease_nw_genes_selected = random.sample(set(disease_nw_genes), len(set(disease_nw_genes)) / 2)
            disease_nw_genes_unselected = set(disease_nw_genes) - set(disease_nw_genes_selected)
            create_ge(disease_nw_genes_selected, bg_genes, dataset_name)
            run_dataset(dataset_name, disease_nw_genes_unselected, cur_disease_name, constants.DEG_EDGER)

        if disease_counter == MAX_DISEASES:
            break

    aggregate_reports.aggregate_disease("_".join(["DISGENET", ts]))
    sys.stdout.write("total tested diseases: {}\n".format(disease_counter))


if __name__ == "__main__":
    constants.update_dirs(DATASET_NAME_u="DISGENET")
    main()
    # aggregate_reports.aggregate_disease("DISGENET_1542711137.18")
