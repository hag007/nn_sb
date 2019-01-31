import json
from matplotlib import style
from pandas._libs.parsers import k

style.use("ggplot")
import seaborn as sns
sns.set(color_codes=True)
import logging
sh = logging.StreamHandler()
logger = logging.getLogger("log")
logger.addHandler(sh)
from constants import *
from infra import *
from param_builder import build_gdc_params
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np




if __name__ == "__main__":
    # for cur_ds in ["GWAS_2hr_glucose", "GWAS_adhd", "GWAS_alzheimers"]:
    #     pval = pd.read_csv(os.path.join(constants.DATASETS_DIR, cur_ds, "data", "score.tsv"), sep='\t')["pval"]
    #     sns.distplot(pval, kde=False)
    #     plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "pval_dist_gwas_pval.png"))

    # for cur_ds in ["GE_SHERA"]:
    #     pval = pd.read_csv(os.path.join(constants.DATASETS_DIR, cur_ds, "cache", "deg_edger.tsv"), sep='\t')["pval"]
    #     sns.distplot(pval, kde=False)
    #     plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "pval_dist_rnaseq_{}.png".format(cur_ds)))

    for cur_ds in ["GE_IEM"]:
        pval = pd.read_csv(os.path.join(constants.DATASETS_DIR, cur_ds, "cache", "deg_t.tsv"), sep='\t')["pval"]
        sns.distplot(pval, kde=False)
        plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "pval_dist_microarray.png"))



