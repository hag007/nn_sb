import os
import shutil
import constants
from utils.network import SH_MODULE_NAME, SH_NUM_GENES
import pandas as pd
import numpy as np

def remove_false_dirs():
    base_folder = "/media/hag007/Data/bnet/output"
    gwas_folders = [os.path.join(base_folder, name) for name in os.listdir(base_folder) if
                      os.path.isdir(os.path.join(base_folder, name)) and name.startswith("GWAS_") ]

    for cur in gwas_folders:
        shutil.rmtree(cur+"/output")
        shutil.rmtree(cur + "/data")
        shutil.rmtree(cur + "/cache")


def modules_report(ds_name, algo_name):
    base_folder = os.path.join(constants.DATASETS_DIR, ds_name, "output")
    module_index = 0
    real_modules = 0
    module_file = os.path.join(base_folder, "{}_module_genes_{}.txt".format(algo_name, module_index))
    modules_summary=[]
    while os.path.exists(module_file):
        num_of_genes=len(file(os.path.join(constants.DATASETS_DIR, ds_name, "output", "{}_module_genes_{}.txt".format(algo_name, module_index))).readlines())
        if  num_of_genes > 3 and os.path.exists(os.path.join(constants.OUTPUT_GLOBAL_DIR,ds_name,algo_name, "module_{}_separated_modules_hg_samples.tsv".format(real_modules))):
            modules_summary.append({SH_MODULE_NAME: real_modules, SH_NUM_GENES: num_of_genes})
            real_modules+=1
        module_index+=1
        module_file = os.path.join(base_folder, "{}_module_genes_{}.txt".format(algo_name, module_index))

    df_summary = pd.DataFrame([], columns=[SH_MODULE_NAME, SH_NUM_GENES])
    if real_modules > 0:
        df_summary = pd.DataFrame(modules_summary).set_index("module")
    df_summary.to_csv(
        os.path.join(constants.OUTPUT_GLOBAL_DIR, ds_name, algo_name, "modules_summary.tsv"), sep="\t")




