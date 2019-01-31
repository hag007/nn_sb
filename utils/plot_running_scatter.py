import os
import constants
import shutil
import json
import numpy as np
import pandas as pd
import ebm
import network
from utils.scripts import format_script
from pandas.io.common import EmptyDataError
from scipy.stats import mannwhitneyu
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.colors as ml_colors
hg_report_header = ["GO id", "value", "pval"]
emb_headers= ["algo", "empirical_brown_score", "enriched_terms"]


if __name__ == "__main__":




    df_results = pd.DataFrame()
    datasets=["GWAS_fasting_insulin", "GWAS_2hr_glucose","GWAS_adhd" ,"GWAS_alzheimers", "GWAS_anorexia", "GWAS_autism", "GWAS_beta-cell_function", "GWAS_bipolar_disorder", "GWAS_blood_pressure_systolic"]
    for cur_ds in datasets:
        constants.update_dirs(DATASET_NAME_u=cur_ds)
        root_path = os.path.join(constants.OUTPUT_GLOBAL_DIR, constants.DATASET_NAME)

        all_algo_modules = {}
        for name in os.listdir(root_path):
            if os.path.isdir(os.path.join(root_path, name)):
                modules_summary = pd.read_csv(os.path.join(root_path, name, "modules_summary.tsv"), sep="\t").set_index(
                    "module")
                all_algo_modules[name] = np.array(
                    modules_summary.index)

        df_scale_col= 'num_of_genes'
        df_score_col = 'score'
        df_module_signal = pd.read_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, constants.DATASET_NAME, "GO_terms_per_module_ratio_best.tsv"), sep='\t')



        df_module_signal = df_module_signal.iloc[df_module_signal[df_score_col].values.argsort(), :]
        algos = np.unique(df_module_signal['algo'].values)
        algos_rank_as_index = df_module_signal['algo'].drop(["bionet"], axis=0).apply(lambda x: np.where(algos == x)[0][0]).values
        cmap = plt.cm.jet
        norm = plt.Normalize(vmin=algos_rank_as_index.min(), vmax=algos_rank_as_index.max())

        fig, ax = plt.subplots(figsize=(15, 5))
        # x = -np.log10(df_mod ule_signal[df_module_col].values)
        x = df_module_signal[df_score_col].drop(["bionet"], axis=0).values
        y = algos_rank_as_index
        s = df_module_signal[df_scale_col].drop(["bionet"], axis=0).values
        ax.scatter(x,y, s, c=y,cmap=cmap)

        colorlist = [ml_colors.rgb2hex(cmap(i)) for i in norm(np.arange(np.size(algos)))]
        patches = [Line2D([0], [0], marker='o', color='gray', label=a,
                          markerfacecolor=c) for i, a, c in zip(list(range(len(algos))), algos, colorlist)]
        ax.legend(handles=list(reversed(patches)), loc='upper left')
        ax.set_xlabel("score")
        ax.set_ylabel("algos")
        ax.set_yticks([])
        plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR,cur_ds, "fraction_score_best_plt.png"))

        plt.cla()


    df_results.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR,"FRACTION_BEST_SCORE.tsv"),sep='\t')


