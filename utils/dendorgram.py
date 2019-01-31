import sys
sys.path.insert(0, '../')

import pandas as pd

from fastsemsim.SemSim import *

import matplotlib.pyplot as plt

from rpy2.robjects import pandas2ri
pandas2ri.activate()

import constants

from scipy.cluster import hierarchy
import scipy.cluster.hierarchy as hcl
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform


import utils.go
import utils.go_hierarcies
import math
import random

from matplotlib.lines import Line2D
import matplotlib.colors as ml_colors
import pylab

from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size

ENABLE_GO_GRAPH = False
IS_GO_GRAPH_ONLY=False
GO_PCA=False
TERMS_SIMILARITY_TO_NUM_OF_TERMS=True
RATIO_TO_GO_TERM = True
QVAL_TH = 0.1

from sklearn.decomposition import PCA

import fastsemsim
from fastsemsim.SemSim import *
from fastsemsim.Ontology import ontologies
from fastsemsim.Ontology import AnnotationCorpus




if __name__ == "__main__":

    my_color = {"2hr_glucose" : 0,
                "adhd": 1,
               "alzheimers": 1,
               "anorexia": 1,
               "autism": 1,
               "beta-cell_function": 2,
               "bipolar_disorder": 1,
               "blood_pressure_systolic": 0,
               "body_mass_index": 0,
               "coronary_artery_disease": 0,
               "crohns_disease": 2,
               "cross_disorder": 1,
               "depressive_disorder": 1,
               "fasting_glucose": 0,
               "fasting_insulin": 0,
               "fasting_proinsulin": 0,
               "glycated_hemoglobin": 0,
               "hdl_cholesterol": 0,
               "height": 4,
               "hepatitis_c_resolvers": 4,
               "inflammatory_bowel_disease": 2,
               "insulin_resistance": 0,
               "insulin_secretion" : 0,
               "ldl_cholesterol": 0,
               "macular_degeneration_dry": 3,
               "macular_degeneration_neovascular": 3,
               "multiple_sclerosis": 2,
               "narcolepsy": 1,
               "osteoporosis": 4,
               "parkinsons_disease": 1,
               "rheumatoid_arthritis": 2,
               "schizophrenia": 1,
                "total_cholesterol": 0,
                "triglycerides": 0,
                "type_1_diabetes": 0,
                "type_2_diabetes": 0,
                "ulcerative_colitis": 2,
                "cancer_1001": 5,
                "cancer_1002": 5,
                "cancer_1019": 5,
                "cancer_1044": 5,
                "cancer_1059": 5,
                "cancer_1061": 5,
                "cancer_1062": 5}





    datasets = [name for name in os.listdir(constants.DATASETS_DIR) if
                os.path.isdir(os.path.join(constants.DATASETS_DIR, name)) and name.startswith("GWAS_")]

    results = pd.DataFrame()
    all_terms_names = []
    # selected_ds = [["depressive_disorder", "fasting_insulin", "rheumatoid_arthritis", "inflammatory_bowel_disease",
    #                "ulcerative_colitis", "fasting_glucose", "crohns_disease"],
    # ["cross_disorder", "hdl_cholesterol", "alzheimers", "blood_pressure_systolic", "parkinsons_disease",
    #  "hepatitis_c_resolvers",
    #  "schizophrenia", "autism", "anorexia", "insulin_resistance", "insulin_secretion", "body_mass_index",
    #  "bipolar_disorder"],
    # ["type_1_diabetes", "type_2_diabetes", "osteoporosis", "triglycerides", "height"]]
    pvals_summaries = []
    # bar_colors=plt.cm.get_cmap("jet",len(selected_ds))
    pval_scores=pd.DataFrame()
    for cur_ds in datasets:

        cur_csv=pd.read_csv(os.path.join(constants.DATASETS_DIR, cur_ds, "data", "score.tsv"), sep='\t', index_col=0)["pval"].to_frame().rename(columns={"pval" : cur_ds[5:]})
        pval_scores=pd.concat((pval_scores, cur_csv), axis=1, join='outer')
        pval_scores[pval_scores.isna()]=1
    pval_scores=pval_scores.apply(lambda x: -np.log10(x))


    fig, ax = plt.subplots(figsize=(20, 10))
    dn = hierarchy.dendrogram(linkage(pval_scores.values.T, metric='correlation' ,method="average"), labels=pval_scores.columns.values)


    # Apply the right color to each label
    my_palette = plt.cm.get_cmap("jet", 6)
    patches = [Line2D([0], [0], marker='o', color='gray', label=l,
                      markerfacecolor=my_palette(i)) for i, l in enumerate(["metabolic", "brain", "inflammatory", "eyes", "others", 'cancer'])]
    ax.legend(handles=patches)

    ax = plt.gca()
    xlbls = ax.get_xmajorticklabels()
    num = -1
    for lbl in xlbls:
        lbl.set_color(my_palette(my_color[lbl.get_text()]))

    plt.tight_layout()
    plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR,"traits_dendogram_by_pval.png"))

    plt.clf()
    plt.cla()
    pca = PCA(n_components=2)
    reduced_df = pca.fit_transform(pval_scores.values.T)
    reduced_df = pd.DataFrame(reduced_df, index=pval_scores.columns.values)

    fig, ax = plt.subplots(figsize=(15, 15))

    sc = ax.scatter(reduced_df.iloc[:, 0].values, reduced_df.iloc[:, 1].values, 50, c=np.array([my_palette(my_color[a]) for a in pval_scores.columns.values]), lw=2)
    # plt.ylim((-1000, 1000))
    # plt.xlim((-1000, 1000))
    ax.set_xscale('symlog')
    ax.set_yscale('symlog')
    ax.set_xlabel("PC_1")
    ax.set_ylabel("PC_2")
    # ax.loglog()
    for k, v in reduced_df.iterrows():
        ax.annotate("{}{}".format(k, "".join(["\n" for x in range(int(random.random() * 5))])), (v.iloc[0], v.iloc[1]))

    patches = [Line2D([0], [0], marker='o', color='gray', label=l,
                      markerfacecolor=my_palette(i)) for i, l in enumerate(["metabolic", "brain", "inflammatory", "eyes", "others", 'cancer'])]
    ax.legend(handles=patches)
    plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "traits_pca_by_pval.png"))

