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
import networkx as nx
import shutil
from pandas.errors import EmptyDataError
from utils.permute_network import EdgeSwapGraph


if __name__ == "__main__":
    prefix="GWAS"
    datasets = [name[5:] for name in os.listdir(constants.DATASETS_DIR) if
                    os.path.isdir(os.path.join(constants.DATASETS_DIR, name)) and name.startswith("GWAS_") and "random" not in name]
    df_score = pd.DataFrame(columns=['pval', "qval"])
    for cur_ds in datasets:
        if prefix=="GE":
            score_file_name="cache/deg_edger.tsv"
            if cur_ds.startswith("IE"):
                score_file_name="cache/deg_t.tsv"
        else:
            score_file_name="data/score.tsv"

        score_file_name=os.path.join(constants.DATASETS_DIR,prefix+"_"+cur_ds,  score_file_name)
        df_score = pd.read_csv(score_file_name, sep='\t').loc[:,["pval","qval"]]

        fig, ax = plt.subplots(figsize=(12, 10))
        ax.set_yscale("log")
        df_score=df_score.loc[df_score["pval"].values != 0, ["pval","qval"]]
        pval=-np.log10(df_score["pval"].values)
        if np.size(pval)==0:
            pval=np.array([1])
        sns.distplot(pval, kde=False)
        plt.title("pval dist. dataset: {}".format(cur_ds))
        real_n_term=np.size(pval)
        plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "score_pval_dist_{}.png".format(cur_ds)))
        plt.clf()

