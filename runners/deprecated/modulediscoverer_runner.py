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
import random
from utils.r_runner import run_rscript

from utils.server import get_parameters

import infra

import utils.go


def prepare_input(method=constants.DEG_EDGER, network_name="dip"):
    deg_file_name = os.path.join(constants.CACHE_DIR, "deg_{}.tsv".format(method))
    network_file_name = os.path.join(constants.NETWORKS_DIR, "{}.sif".format(network_name))

    network_df = pd.read_csv(network_file_name, sep="\t")
    src = np.array(network_df["ID_interactor_A"])
    dst = np.array(network_df["ID_interactor_B"])

    vertices = list(set(np.append(src,dst)))
    A = np.zeros((len(vertices), len(vertices)))

    v_list_data = np.c_[np.array([i+1 for i in range(len(vertices))]), np.ones(len(vertices)), np.zeros(len(vertices))]
    vlist = pd.DataFrame(v_list_data,index = [i for i in range(len(vertices))], columns=['content', 'weight', 'degree'], dtype=np.int)


    for i, cur_r in network_df.iterrows():
        A[vertices.index(cur_r["ID_interactor_A"]),vertices.index(cur_r["ID_interactor_B"])]=1
        A[vertices.index(cur_r["ID_interactor_B"]), vertices.index(cur_r["ID_interactor_A"])] = 1
        vlist.loc[list(set([vertices.index(cur_r["ID_interactor_B"]), vertices.index(cur_r["ID_interactor_A"])])),["degree"]]+=1


    deg = infra.load_gene_expression_profile_by_genes(gene_expression_path=deg_file_name)
    h_rows, h_cols, deg_data = infra.separate_headers(deg)
    ind = np.where(h_cols=="qval")[0][0]
    ordered_ind = np.argsort(deg_data[:,ind])
    deg_data=deg_data[ordered_ind,:]
    h_rows=h_rows[ordered_ind]
    sig_last_index = np.where(deg_data[:,np.where(h_cols=="qval")[0][0]]>0.05)[0][0]

    degs = list(set(h_rows[:sig_last_index]).intersection(vertices))
    background = list(set(vertices).intersection(set(h_rows)))
    random_sets = [random.sample(vertices, len(degs)) for x in range(10000)]

    bg_genes = vertices
    bg_genes_file_name=os.path.join(constants.OUTPUT_DIR, "keypathwayminer_bg_genes.txt")
    file(os.path.join(constants.OUTPUT_DIR, bg_genes_file_name), "w+").write("\n".join(bg_genes))


    file(os.path.join(constants.OUTPUT_DIR,"A"),"w+").write("\n".join(["\t".join([str(y) for y in x])  for x in A]))
    file(os.path.join(constants.OUTPUT_DIR, "degs"),"w+").write("\n".join(degs))
    file(os.path.join(constants.OUTPUT_DIR, "proteins"),"w+").write("\n".join(vertices))
    file(os.path.join(constants.OUTPUT_DIR, "vlist"),"w+").write("\t".join(vlist.columns)+"\n"+"\n".join(["\t".join([str(y) for y in x]) for i,x in vlist.iterrows()]))
    file(os.path.join(constants.OUTPUT_DIR, "background"),"w+").write("\n".join(background))
    file(os.path.join(constants.OUTPUT_DIR, "random_sets"),"w+").write("\n".join(["\t".join(x) for x in random_sets]))

    return network_file_name, bg_genes, \
           os.path.join(constants.OUTPUT_DIR, "A"), os.path.join(constants.OUTPUT_DIR,"degs"), \
           os.path.join(constants.OUTPUT_DIR, "proteins"), os.path.join(constants.OUTPUT_DIR,"vlist"), \
           os.path.join(constants.OUTPUT_DIR, "background"), os.path.join(constants.OUTPUT_DIR,"random_sets"),\
           0.05, os.path.join(constants.OUTPUT_DIR, "moduledicoverer_output.txt")


def run_modulediscoverer(A_file, degs_file, proteins_file, vlist_file, background_file, random_sets_file, p_value, output_file):
    script = file("../r/scripts/ModuleDiscovererMain.r").read()
    return run_rscript(script=script, A_file=A_file, degs_file=degs_file, proteins_file=proteins_file, vlist_file=vlist_file, background_file=background_file, random_sets_file=random_sets_file, p_value=p_value, output_file=output_file)

if __name__ == "__main__":
    params = get_parameters()
    if params != None:
        args, NETWORK_NAME, dataset_name = params

    score_method=constants.DEG_EDGER
    network_file_name, bg_genes, A_file, degs_file, proteins_file, vlist_file, background_file, random_sets_file, p_value, output_file = prepare_input(method=score_method)
    run_modulediscoverer(A_file, degs_file, proteins_file, vlist_file, background_file, random_sets_file, p_value, output_file)

    results = file(output_file).readlines()
    results = [y for x in results for y in x.strip().split()]
    module_genes = list(set(results))
    file(os.path.join(constants.OUTPUT_DIR,"modulediscoverer_module_genes.txt"), "w+").write("\n".join(module_genes))
    utils.go.check_group_enrichment(module_genes, bg_genes)






