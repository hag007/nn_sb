import sys
sys.path.insert(0, '../')

import argparse
import constants
import os
import time
import shutil
import pandas as pd
from runners import DEG_runner
from utils.ensembl2gene_symbol import e2g_convertor
import numpy as np
import json
from network import get_network_genes

def get_parameters():
    return None

# def get_parameters():
#     parser = argparse.ArgumentParser(description='args from wev client')
#     parser.add_argument('--gene_expression', dest='ge', default="")
#     parser.add_argument('--score', dest='score', default="")
#     parser.add_argument('--classes', dest='classes', default="")
#     parser.add_argument('--network', dest='nw', default="")
#     parser.add_argument('--reports', dest='reports', default=False)
#     parser.add_argument('--server_mode', dest='server_mode', default=False)
#     parser.add_argument('--mode', dest='mode', default=0)
#     args = parser.parse_args()
#     if args.mode==0:
#         return None
#     NETWORK_NAME = os.path.splitext(os.path.basename(args.nw))[0]
#     dataset_name = "user" + str(time.time())
#     constants.update_dirs(DATASET_NAME_u=dataset_name)
#     os.makedirs(os.path.join(constants.DATA_DIR, "data"))
#     os.makedirs(os.path.join(constants.OUTPUT_DIR, "output"))
#     os.makedirs(os.path.join(constants.CACHE_DIR, "cache"))
#
#     if args.ge != "" and args.ge != os.path.join(constants.DATA_DIR, os.path.basename(args.ge)):
#         shutil.copy(args.ge, os.path.join(constants.DATA_DIR, "ge.tsv"))
# 
#     if args.score != "" and  args.score != os.path.join(constants.DATA_DIR, os.path.basename(args.score)):
#         shutil.copy(args.score, os.path.join(constants.DATA_DIR, "score.tsv"))
# 
#     if args.classes != "" and args.classes != os.path.join(constants.DATA_DIR, os.path.basename(args.classes)):
#         shutil.copy(args.classes, os.path.join(constants.DATA_DIR,  "classes.tsv"))
# 
#     if args.nw != "" and args.nw != os.path.join(constants.NETWORKS_DIR, NETWORK_NAME + ".sif"):
#         NETWORK_NAME = os.path.splitext(os.path.basename(args.nw))[0]
#         shutil.copy(args.nw, os.path.join(constants.NETWORKS_DIR, NETWORK_NAME + ".sif"))
# 
#     if args.reports == "true":
#         constants.REPORTS = True
#         constants.HG_MODE = True
# 
#     if args.server_mode == "true":
#         constants.SERVER_MODE = True
# 
# 
#     return args, NETWORK_NAME, dataset_name

def get_score(score_method):
    if score_method != constants.PREDEFINED_SCORE:
        score_file_name = os.path.join(constants.CACHE_DIR, "deg_{}.tsv".format(score_method))
        if not os.path.exists(score_file_name):
            DEG_runner.main(method=score_method)
    else:
        score_file_name = os.path.join(constants.DATA_DIR, "score.tsv".format(score_method))
        df = pd.read_csv(score_file_name, sep="\t")
        if "pval" in df and "qval" in df:
            pass
        elif "score" in df:
            df.rename(columns={"score": "pval"})
            df["qval"] = df["pval"]
            df.to_csv(score_file_name, sep="\t", index=False)
            constants.IS_PVAL_SCORES = False
        else:
            raise Exception("expected cols: pval, qval or score. got {}".format(df.columns))

    return score_file_name

def init_common_params(NETWORK_NAME, score_method = constants.DEG_EDGER):
    params = get_parameters()
    if params != None:
        args, NETWORK_NAME, dataset_name = params
        if (args.score == ""):
            score_method = constants.DEG_EDGER
        else:
            score_method = constants.PREDEFINED_SCORE

    network_file_name = os.path.join(constants.NETWORKS_DIR, NETWORK_NAME)
    score_file_name = None
    if score_method is not None:
        score_file_name = get_score(score_method)
    network_genes = get_network_genes()
    return network_file_name, score_file_name, score_method, network_genes

