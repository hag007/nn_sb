import constants
import infra
import pandas as pd
import numpy as np
import os
import shutil
from utils.ensembl2gene_symbol import g2e_convertor

scz_source_folder= "/media/hag007/Data/bnet/datasets/GWAS_scz_updated/data"

scz_scores = [os.path.join(scz_source_folder, name) for name in os.listdir(scz_source_folder) if not os.path.isdir(os.path.join(scz_source_folder, name))]
for cur_scz_score in scz_scores:
    if os.path.basename(cur_scz_score) == "_README.txt":
        continue
    print "current file: {}".format(os.path.basename(cur_scz_score))
    df = pd.read_csv(cur_scz_score, sep='\t')
    df["id"] = pd.Series([g2e_convertor([cur])[0] if len(g2e_convertor([cur])) > 0 else np.NaN for cur in df["gene_symbol"].values])
    df=df.dropna()
    df.index = df["id"]
    df=df[["pvalue"]]
    df=df.rename(columns={"pvalue": "pval"})
    df["qval"] = df["pval"]
    df = df[~df.index.duplicated(keep='first')]
    new_ds_dir = os.path.join(constants.DATASETS_DIR, "GWAS_" + os.path.basename(cur_scz_score)[3:].split(".")[0])
    # os.makedirs(os.path.join(new_ds_dir,"data"))
    # os.makedirs(os.path.join(new_ds_dir, "output"))
    # os.makedirs(os.path.join(new_ds_dir, "cache"))

    df.to_csv(os.path.join(os.path.join(new_ds_dir,"data"), "score.tsv"),sep='\t')

