import json
import sys
sys.path.insert(0, "../")
import matplotlib as mpl 
mpl.use("Agg")
import os
import pandas as pd
import seaborn as sns
sns.set(color_codes=True)
import logging
sh = logging.StreamHandler()
logger = logging.getLogger("log")
logger.addHandler(sh)
from constants import *
from infra import *
import matplotlib.pyplot as plt
import time


ax = plt.subplot(111)
total_data = []
for cur_cancer_type in ["ICGC_LICA_FR", "ICGC_LIHC_US"]: # , "ICGC_LIRI_JP", "ICGC_PAEN_AU"]: # ["ICGC_PRAD_US", "ICGC_PAAD_US", "ICGC_LUSC_US", "ICGC_COAD_US"]:
    file_path=os.path.join(constants.DATASETS_DIR, cur_cancer_type, "data", "{}.htseq_rsem.tsv".format(cur_cancer_type))
    # file_path=os.path.join(constants.DATASETS_DIR, cur_cancer_type, "data", "exp_seq.{}-{}.tsv".format(cur_cancer_type.split("_")[1], cur_cancer_type.split("_")[2]))
    df_ge=pd.read_csv(file_path, sep='\t', index_col=0)
    group_conditions = None 
    flatten_data=np.power([2 for a in df_ge.values.flatten()], df_ge.values.flatten()) # df_ge["normalized_read_count"].values  
    flatten_data=flatten_data#-0.99
    flatten_data=np.log2(flatten_data)
    sns.distplot(flatten_data)
    plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, cur_cancer_type+"_dist.png"))
    plt.cla()
