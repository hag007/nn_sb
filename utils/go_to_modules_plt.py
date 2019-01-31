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
from scipy.spatial.distance import squareform


import utils.go
import utils.go_hierarcies
import math
import random

from matplotlib.lines import Line2D
import matplotlib.colors as ml_colors

from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size

ENABLE_GO_GRAPH = False
IS_GO_GRAPH_ONLY=False
GO_PCA=False
TERMS_SIMILARITY_TO_NUM_OF_TERMS=True
RATIO_TO_GO_TERM = True
QVAL_TH = 0.01
SIM_TH= 0.4

from sklearn.decomposition import PCA

import fastsemsim
from fastsemsim.SemSim import *
from fastsemsim.Ontology import ontologies
from fastsemsim.Ontology import AnnotationCorpus


##########################
if TERMS_SIMILARITY_TO_NUM_OF_TERMS:
    ontology_type = 'GeneOntology'
    ignore_parameters = {'ignore': {}}
    source_type = 'obo'
    source = os.path.join(os.path.join(constants.GO_DIR, constants.GO_FILE_NAME))

    print "\n######################"
    print "# Loading ontology... #"
    print "######################\n"

    ontology = ontologies.load(source=source, source_type=source_type, ontology_type=ontology_type,
                               parameters=ignore_parameters)

    print "\n######################"
    print "# Loading Annotation Corpus... #"
    print "######################\n"
    ac = AnnotationCorpus.AnnotationCorpus(ontology)
    ac.parse(os.path.join(constants.GO_DIR, "goa_human.gaf"), "gaf-2.0")
    ac.isConsistent()

    print "\n#################################"
    print "# Annotation corpus successfully loaded."
    print "#################################\n"

    semsim = GSESAMESemSim(ontology, ac)  # maxSemSim(ontology, ac) #


#################


def main(datasets, algos, roots):
    dict_result, go2geneids, geneids2go, entrez2ensembl = utils.go_hierarcies.build_hierarcy(roots=roots) #  'GO:0006935', 'GO:0002376','GO:0003008', 'GO:0007610' 'GO:0032502'
    dict_go_terms = {}
    terms_for_traits=[]
    for x in dict_result.values():
        dict_go_terms.update(x['vertices'])

    all_hg_score_modules = []
    dict_terms_by_traits = {}
    df_go_metadata = pd.DataFrame()
    all_df_pvals = []

    for cur_ds in datasets:
        df_all_hg_pval = pd.DataFrame()
        df_all_hg_size = pd.DataFrame()
        cur_ds_terms = pd.DataFrame()
        terms_count = 0
        for cur_algo in algos:

            summary_file_name = os.path.join(constants.OUTPUT_GLOBAL_DIR, cur_ds, cur_algo, "modules_summary.tsv")
            if not os.path.exists(summary_file_name):
                continue
            df_modules_summary = pd.read_csv(
                summary_file_name,
                sep='\t')
            i=-1
            for i in range(len(df_modules_summary.index)):

                hg_file_name = os.path.join(constants.OUTPUT_GLOBAL_DIR, cur_ds, cur_algo,
                                            "module_{}_separated_modules_hg_samples.tsv".format(i))

                print "reading module: {} {} {} from file".format(cur_ds, cur_algo, i)
                if os.path.getsize(hg_file_name) < 2:
                    continue


                else:
                    df_hg_output = pd.read_csv(hg_file_name, sep="\t")
                    df_go_metadata = pd.concat([df_go_metadata, df_hg_output[["GO id", "GO name"]]], axis=0)

                    df_hg_output.index = df_hg_output["GO name"]
                    df_hg_output=df_hg_output[df_hg_output["GO id"].isin(dict_go_terms.keys())]
                    df_hg_output=df_hg_output[df_hg_output["qval"].values <= QVAL_TH]
                    df_hg_output.loc[df_hg_output["GO name"].apply(
                        lambda x: not ("rna" in x.lower() or "dna" in x.lower() or "actine" in x.lower())).values]
                    df_hg_output = df_hg_output[df_hg_output["GO id"].apply(lambda x :(dict_go_terms[x]['n_children'] == 0) or
                                                                                  (dict_go_terms[x]['D'] > 7
                                                                                       ))]
                    if df_hg_output.empty: continue
                    df_hg_output = df_hg_output[df_hg_output["GO id"].apply(lambda x : x in go2geneids and len(go2geneids[x])<=500)]
                    if df_hg_output.empty: continue

                    cur_series_size = df_hg_output["value"].to_frame() \
                        .rename(columns={"value": "{}_{}".format(cur_algo, i)})
                    cur_series_size=cur_series_size/float(df_modules_summary.loc[i, "#_genes"])


                    cur_series_pval = df_hg_output["pval"].apply(lambda x: -np.log10(x)).to_frame() \
                        .rename(columns={"pval": "{}_{}".format(cur_algo, i)})

                    df_all_hg_size = pd.concat([df_all_hg_size,
                                                cur_series_size],
                                               join='outer', axis=1)

                    df_all_hg_pval = pd.concat([df_all_hg_pval,
                                                cur_series_pval],
                                               join='outer', axis=1)

                    df_all_hg_pval[pd.isna(df_all_hg_pval)] = 0

        df_go_metadata.index = df_go_metadata["GO name"]
        df_go_metadata = df_go_metadata[~df_go_metadata.index.duplicated(keep='first')]
        dict_terms_by_traits[cur_ds] = df_go_metadata.loc[df_all_hg_pval.index.values, "GO id"].to_frame()
        dict_terms_by_traits[cur_ds].index.name="GO name"
        dict_terms_by_traits[cur_ds]["pval"]=df_all_hg_pval.max(axis=1)
        dict_terms_by_traits[cur_ds]=dict_terms_by_traits[cur_ds].set_index("GO id")

        terms_for_traits.append({"ds_name": cur_ds, "terms_count": np.sum(df_all_hg_pval.values!=0)})
        if df_all_hg_pval.empty: continue

        sort_GO_terms=df_all_hg_pval.sum(axis=1).sort_values().index.values
        df_all_hg_pval=df_all_hg_pval.loc[sort_GO_terms,:]
        df_all_hg_size=df_all_hg_size.loc[sort_GO_terms, :]

        # plot_single_gwas_scatter(cur_ds, df_all_hg_pval, df_all_hg_size)

        all_df_pvals.append(df_all_hg_pval)
    print "# GO terms: {}, # of modules: {}".format(len(df_all_hg_pval.index),i+1)
    return terms_for_traits, dict_terms_by_traits, [dict_go_terms[cur_go]['name'] for cur_go in roots], all_df_pvals


def plot_single_gwas_scatter(cur_ds, df_all_hg_pval, df_all_hg_size):
    fig, ax = plt.subplots(figsize=(15, 15))

    x = np.arange(len(df_all_hg_pval.columns))
    y = np.arange(len(df_all_hg_pval.index))
    xx, yy = zip(*[(a, b) for a in x for b in y if df_all_hg_pval.iloc[b, a] != 0])
    c = [df_all_hg_pval.iloc[b, a] for a in x for b in y if df_all_hg_pval.iloc[b, a] != 0]
    s = [df_all_hg_size.iloc[b, a] for a in x for b in y if df_all_hg_pval.iloc[b, a] != 0]
    s = np.array(s)
    # im = plt.imshow(np.array(c).reshape(len(c), 1), cmap='bwr')
    # im.remove()
    sc = ax.scatter(xx, yy, s / float(max(s)) * 1000, c=c, cmap='bwr', vmin=np.percentile(c, 10),
                    vmax=np.percentile(c, 90))
    for s_i, a in enumerate(s):
        ax.annotate("%.2f" % a, (xx[s_i], yy[s_i]))
    ax.legend(loc='upper left')
    ax.margins(0.03, 0.03)
    # ax.locator_params(axis='x', nbins=len(df_all_hg_pval.columns))
    # ax.locator_params(axis='y', nbins=len(df_all_hg_pval.index))
    ax.set_xlabel("modules")
    plt.subplots_adjust(left=0.25, right=0.99, top=0.99, bottom=0.05)
    plt.xticks(np.arange(len(df_all_hg_pval.columns)), tuple(list(df_all_hg_pval.columns.values)), rotation='vertical')
    ax.set_ylabel("GO terms")
    plt.yticks(np.arange(len(df_all_hg_pval.index)), tuple(list(df_all_hg_pval.index.values)))
    ax_ = plt.gca()
    aspect = 20
    pad_fraction = 0.5
    divider = make_axes_locatable(ax_)
    width = axes_size.AxesY(ax_, aspect=1. / aspect)
    pad = axes_size.Fraction(pad_fraction, width)
    cax = divider.append_axes("right", size=0.3, pad=0.4)
    plt.colorbar(mappable=sc, cax=cax)
    plt.tight_layout()
    plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "go_to_modules_{}.png".format(cur_ds)))
    plt.cla()


def plot_sumary_scatter(df, all_terms_names):

    reduced_df = df
    if len(df.columns) > 2:
        pca = PCA(n_components=2)
        reduced_df = pca.fit_transform(df.values)
        reduced_df = pd.DataFrame(reduced_df, index=df.index)


    fig, ax = plt.subplots(figsize=(15, 50))
    sc = ax.scatter(reduced_df.iloc[:,0].values, reduced_df.iloc[:,1].values, 50)
    for k, v in reduced_df.iterrows():
        ax.annotate("{}{}".format(k,"".join(["\n" for x in range(int(random.random()*5))])), (v.iloc[0], v.iloc[1]))
    ax.legend(loc='upper left')
    ax.margins(0.12, 0.12)
    # ax.locator_params(axis='x', nbins=len(df_all_hg_pval.columns))
    # ax.locator_params(axis='y', nbins=len(df_all_hg_pval.index))
    ax.set_xlabel("development")
    ax.set_ylabel("immune")

    if len(df.columns) > 2:
        coeff = np.transpose(pca.components_[0:2, :]) * 150
        labels = df.columns.values
        top_pca_x_go_terms = labels[np.abs(coeff[:, 0]).argsort()[::-1][:5]]
        top_pca_x_details = ["{}: {}".format(x, all_terms_names[np.where(df.columns.values==x)[0][0]]) for x in top_pca_x_go_terms] # df_go_metadata.loc[x, "GO name"]
        top_pca_y_go_terms = labels[np.abs(coeff[:, 1]).argsort()[::-1][:5]]
        top_pca_y_details = ["{}: {}".format(x, all_terms_names[np.where(df.columns.values==x)[0][0]]) for x in top_pca_y_go_terms] # df_go_metadata.loc[x, "GO name"]
        props = dict(boxstyle='round', facecolor='gray', alpha=0.5)
        fig.text(0.01, 0.01,
                 "\n".join(["PC1 top GO terms:"] + top_pca_x_details + ["PC2 top GO terms:"] + top_pca_y_details),
                 bbox=props)

    plt.tight_layout()
    plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "terms_to_trait_summary_scatter.png"))
    plt.cla()


def single_similarity(linkage_matrix, dict_terms_by_traits, cur_row, cur_col):
    linkage_matrix.loc[cur_row, cur_col] = 0
    total_terms_score = 0
    for i, cur_term_row in enumerate(dict_terms_by_traits[cur_row]["GO id"].values):
        print "cur term index: {}".format(i)
        cur_term_score = 0
        for j, cur_term_col in enumerate(dict_terms_by_traits[cur_col]["GO id"].values):
            cur_term_score = max(cur_term_score,
                                 semsim.SemSim(cur_term_row, cur_term_col))
        total_terms_score += cur_term_score
    linkage_matrix.loc[cur_row, cur_col] = 1 - (total_terms_score / max(1, len(dict_terms_by_traits[cur_row]["GO id"].values)))


def non_weighted_average_similarity(linkage_matrix, dict_terms_by_traits, cur_row, cur_col):

    linkage_matrix.loc[cur_row, cur_col] = 0
    cur_term_score = 0
    for i, cur_term_row in enumerate(dict_terms_by_traits[cur_row].index.values):
        print "cur term index: {}".format(i)
        for j, cur_term_col in enumerate(dict_terms_by_traits[cur_col].index.values):
            cur_term_score += semsim.SemSim(cur_term_row, cur_term_col)
    linkage_matrix.loc[cur_row, cur_col] = cur_term_score / max(
        len(dict_terms_by_traits[cur_row].index) * len(dict_terms_by_traits[cur_col].index), 1)
    linkage_matrix.loc[cur_row, cur_col] = 1 - linkage_matrix.loc[cur_row, cur_col]


def weighted_maximal_similarity(linkage_matrix, dict_terms_by_traits, cur_row, cur_col):

    original_row=cur_row
    original_col=cur_col
    if len(dict_terms_by_traits[cur_col].index) > len(dict_terms_by_traits[cur_row].index):
        temp=cur_col
        cur_col=cur_row
        cur_row=temp

    linkage_matrix.loc[cur_row, cur_col] = 0
    total_terms_score = 0
    total_weights = 0
    for i, cur_term_row in enumerate(dict_terms_by_traits[cur_row].index.values):
        print "cur term index: {}".format(i)
        max_term_score=0
        max_weight=0

        for j, cur_term_col in enumerate(dict_terms_by_traits[cur_col].index.values):
            cur_weight=dict_terms_by_traits[cur_row].loc[cur_term_row, "pval"] * \
            dict_terms_by_traits[cur_col].loc[cur_term_col, "pval"]
            sim_score=semsim.SemSim(cur_term_row, cur_term_col)
            cur_term_score = (sim_score ** 3) * cur_weight

            if cur_term_score > max_term_score or True:
                max_term_score += round(cur_term_score,4)
                max_weight += round(cur_weight,4)
                if max_term_score/max_weight >1:
                    print "invalid values: {}, {}, {}, {}, {}, {}".format(cur_term_row, cur_term_col, max_term_score, max_weight, max_term_score/max_weight, sim_score)
                    exit(1)

        total_terms_score += max_term_score
        total_weights += max_weight
    linkage_matrix.loc[original_row, original_col] = 1-(total_terms_score/float(total_weights) if total_weights > 0 else 0)



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
                "cancer_1062": 5

                }


    # go_terms_level_2 = ['GO:0001909', 'GO:0030447', 'GO:0050817', 'GO:0007340', 'GO:0031503', 'GO:0001775', 'GO:0008219', 'GO:0055114', 'GO:0044851', 'GO:0061919', 'GO:0022601', 'GO:0022602', 'GO:0022608', 'GO:0007155', 'GO:0007154', 'GO:0006276', 'GO:0022413', 'GO:0022412', 'GO:0070341', 'GO:0045494', 'GO:0044033', 'GO:0045103', 'GO:0002200', 'GO:0048647', 'GO:0065008', 'GO:0140029', 'GO:0043696', 'GO:0071554', 'GO:0035074', 'GO:0035073', 'GO:0043335', 'GO:0072089', 'GO:0051821', 'GO:0140253', 'GO:0007059', 'GO:0030431', 'GO:0007631', 'GO:0007635', 'GO:0007638', 'GO:0014854', 'GO:0050673', 'GO:0010118', 'GO:0032898', 'GO:0035988', 'GO:0071625', 'GO:0044085', 'GO:0002174', 'GO:0061744', 'GO:0010312', 'GO:0009719', 'GO:0010014', 'GO:0008340', 'GO:0045730', 'GO:0072111', 'GO:0060033', 'GO:0061842', 'GO:0030588', 'GO:0009791', 'GO:0051410', 'GO:0060206', 'GO:0021700', 'GO:0002252', 'GO:0001502', 'GO:0001503', 'GO:0032963', 'GO:0045058', 'GO:0051780', 'GO:0071335', 'GO:0002941', 'GO:0044366', 'GO:0048134', 'GO:0051703', 'GO:0051705', 'GO:0002532', 'GO:0007568', 'GO:0007566', 'GO:0070661', 'GO:0048856', 'GO:0006807', 'GO:0010022', 'GO:0052192', 'GO:0032196', 'GO:0022004', 'GO:0022005', 'GO:0033002', 'GO:0002404', 'GO:0060722', 'GO:0048589', 'GO:0009845', 'GO:0009847', 'GO:0009846', 'GO:0051674', 'GO:1990110', 'GO:1990654', 'GO:0009561', 'GO:0009566', 'GO:0019882', 'GO:0010073', 'GO:0007272', 'GO:0009058', 'GO:0009056', 'GO:0044281', 'GO:0060352', 'GO:0097242', 'GO:0034337', 'GO:0022406', 'GO:0022404', 'GO:0022403', 'GO:0070085', 'GO:0071684', 'GO:0033060', 'GO:0097528', 'GO:0070988', 'GO:0032259', 'GO:1901902', 'GO:0090713', 'GO:0033059', 'GO:0033058', 'GO:0031987', 'GO:0006457', 'GO:0060756', 'GO:0003006', 'GO:0003008', 'GO:0032537', 'GO:0072690', 'GO:0006950', 'GO:0006955', 'GO:0007389', 'GO:0051641', 'GO:0060384', 'GO:0098727', 'GO:0080189', 'GO:0034381', 'GO:0002339', 'GO:0050789', 'GO:0009838', 'GO:0048532', 'GO:0061323', 'GO:0048646', 'GO:0042630', 'GO:0016203', 'GO:0007117', 'GO:0065009', 'GO:0030534', 'GO:0030537', 'GO:0019954', 'GO:0042221', 'GO:0036268', 'GO:0014009', 'GO:0048609', 'GO:0036093', 'GO:0032504', 'GO:0032505', 'GO:0044419', 'GO:0007017', 'GO:0009628', 'GO:0044237', 'GO:0042303', 'GO:0007586', 'GO:0007585', 'GO:1902579', 'GO:0044238', 'GO:0050879', 'GO:0003419', 'GO:0051867', 'GO:0035172', 'GO:0060273', 'GO:0036166', 'GO:2000793', 'GO:0006903', 'GO:0022402', 'GO:0090255', 'GO:0051716', 'GO:0044764', 'GO:0002507', 'GO:1901275', 'GO:0016049', 'GO:0016043', 'GO:0001816', 'GO:0022611', 'GO:0048144', 'GO:0007163', 'GO:0097278', 'GO:0098630', 'GO:0061887', 'GO:1990845', 'GO:0036363', 'GO:0048266', 'GO:0051450', 'GO:0090675', 'GO:1990748', 'GO:0071704', 'GO:0042330', 'GO:0033036', 'GO:0000920', 'GO:0019835', 'GO:0009856', 'GO:0042440', 'GO:0007049', 'GO:0060242', 'GO:0009653', 'GO:0007626', 'GO:0007625', 'GO:0007624', 'GO:0007623', 'GO:0007622', 'GO:0090130', 'GO:0035637', 'GO:0035638', 'GO:0010127', 'GO:0016037', 'GO:0060361', 'GO:0044663', 'GO:0008037', 'GO:0030029', 'GO:0035889', 'GO:0035314', 'GO:0010463', 'GO:0048066', 'GO:0009607', 'GO:0009605', 'GO:0033687', 'GO:0014874', 'GO:0044110', 'GO:0044111', 'GO:0051234', 'GO:0097360', 'GO:0080190', 'GO:0002440', 'GO:0043934', 'GO:0044703', 'GO:0044706', 'GO:0061687', 'GO:0006928', 'GO:0061948', 'GO:0071838', 'GO:0035640', 'GO:0051775', 'GO:1902421', 'GO:0055127', 'GO:0061351', 'GO:0090618', 'GO:0051301', 'GO:0007571', 'GO:0019748', 'GO:0001833', 'GO:0001834', 'GO:0043627', 'GO:0097737', 'GO:0097194', 'GO:0035726', 'GO:0044406', 'GO:0035187', 'GO:0048869', 'GO:0031424', 'GO:0071722', 'GO:0051606', 'GO:0014823', 'GO:0055046', 'GO:0035732', 'GO:0035736', 'GO:0048771', 'GO:0031128', 'GO:0051816']
    # go_terms_sets = [[x] for x in go_terms_level_2]
    # go_terms_sets = [['GO:0006935', 'GO:0002376', 'GO:0003008', 'GO:0032502']] # ,[]] #
    # go_terms_sets = [["GO:0051819", "GO:0010804", "GO:0010803", "GO:0071706", "GO:0071228", "GO:1903265", "GO:1990774", "GO:0080034", "GO:0042533", "GO:0042534", "GO:0042535", "GO:0042536", "GO:0032640", "GO:0032680", "GO:0032760", "GO:0032720", "GO:0034612", "GO:0002834", "GO:0002835", "GO:0002836", "GO:0002837", "GO:0002838", "GO:0002839", "GO:0002840", "GO:0002841", "GO:0002842", "GO:0002843", "GO:0002844", "GO:0002845", "GO:0002846", "GO:0002847", "GO:0002848", "GO:0002852", "GO:0002853", "GO:0002854", "GO:0002855", "GO:0002856", "GO:0002857", "GO:0002858", "GO:0002859", "GO:0002860", "GO:0002411", "GO:0002413", "GO:0002418", "GO:0002419", "GO:0002420", "GO:0002423", "GO:0002424", "GO:0002347", "GO:0002355", "GO:0002357", "GO:0033209", "GO:0071356", "GO:0072535", "GO:1903555", "GO:1903556", "GO:1903557", "GO:2000309", "GO:2000308", "GO:2000307", "GO:1904469", "GO:1904468", "GO:1904467", "GO:0051726", "GO:0036462", "GO:0039527", "GO:0039547", "GO:0036337", "GO:0038008", "GO:0071847", "GO:1903984", "GO:0043017", "GO:0043016", "GO:0043018", "GO:0045558", "GO:0045559", "GO:0045553", "GO:0042109", "GO:0032641", "GO:0032681", "GO:0032761", "GO:0032721", "GO:1903122", "GO:1903121"]]# "GO:0050789", "GO:0044237", "GO:0071704", "GO:0044238", "GO:0006807"]]
    # go_terms_sets = [['GO:0006935']]
    go_terms_sets = [['GO:0002376']]
    algos = ["keypathwayminer_INES_GREEDY", "netbox", "hotnet2", "jactivemodules_greedy", "bionet",
              "jactivemodules_sa"]  # "matisse", "reactomefi"  # "keypathwayminer_INES_OPTIMAL", "keypathwayminer_INES_ACO"
    datasets = [name for name in os.listdir(constants.DATASETS_DIR) if
                os.path.isdir(os.path.join(constants.DATASETS_DIR, name)) and name.startswith("GWAS_")]



    results = pd.DataFrame()
    all_terms_names = []
    selected_ds = [["depressive_disorder", "fasting_insulin", "rheumatoid_arthritis", "inflammatory_bowel_disease",
                   "ulcerative_colitis", "fasting_glucose", "crohns_disease"],
    ["cross_disorder", "hdl_cholesterol", "alzheimers", "blood_pressure_systolic", "parkinsons_disease",
     "hepatitis_c_resolvers",
     "schizophrenia", "autism", "anorexia", "insulin_resistance", "insulin_secretion", "body_mass_index",
     "bipolar_disorder"],
    ["type_1_diabetes", "type_2_diabetes", "osteoporosis", "triglycerides", "height"],
   ["cancer_1001", "cancer_1002", "cancer_1019", "cancer_1044", "cancer_1059", "cancer_1061", "cancer_1062"]]
    # selected_ds = [[y for x in selected_ds for y in x]]
    selected_ds = [[x[5:] for x in datasets]]
    pvals_summaries = []
    bar_colors=plt.cm.get_cmap("jet",len(selected_ds))

    for cur_selected_ds in selected_ds:
        for i, cur_roots in enumerate(go_terms_sets):

            cur_result, dict_terms_by_traits, terms_names, all_df_pvals = main([cur for cur in datasets if cur[5:] in cur_selected_ds], algos, cur_roots)

            results =pd.concat((results, pd.DataFrame(cur_result).set_index("ds_name")[["terms_count"]].rename(columns={"terms_count": " ".join(cur_roots)})), axis=1)
            if len(all_df_pvals)!=0:
                pvals_summary=pd.concat(tuple(all_df_pvals),join='outer', axis=1).max(axis=1)
                pvals_summaries.append(pvals_summary)

    w = 0.15
    fig, ax = plt.subplots(figsize=(20, 10))
    y_axis_go_terms = list(set([y for cur in pvals_summaries for y in list(cur.index)]))
    plt.yticks(np.arange(len(y_axis_go_terms))+w, y_axis_go_terms)
    for i, cur_pval_summary in enumerate(pvals_summaries):
        plt.barh(np.arange(len(y_axis_go_terms))+w*i, [cur_pval_summary.loc[cur] if cur in list(cur_pval_summary.index) else 0 for cur in y_axis_go_terms], height=w , align='center', color=bar_colors(i))

    ax.set_xlabel("-log10(p-val)")
    patches = [Line2D([0], [0], marker='o', color='gray', label=l,
                      markerfacecolor=bar_colors(i)) for i, l in
               enumerate(["c_inflamatory","c_brain","c_other", "c_cancer"])]
    plt.legend(handles=patches)
    plt.tight_layout()
    # plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR,"go_pval_summary_ds_agg.png"))

    # exit(0)
    linkage_matrix=pd.DataFrame(columns=dict_terms_by_traits.keys(), index=dict_terms_by_traits.keys())
    linkage_matrix.loc[:,:]=np.inf
    for row_i, cur_row in enumerate(dict_terms_by_traits.keys()):
        for col_i, cur_col in enumerate(dict_terms_by_traits.keys()):

            print "calc similaritry between {} (i={}, n={}), {} (i={}, n={})".format(cur_row, row_i, len(dict_terms_by_traits[cur_row].index), cur_col, col_i, len(dict_terms_by_traits[cur_col].index))
            if cur_col == cur_row:
                print "same trait. distance is 0"
                linkage_matrix.loc[cur_col, cur_row] = 0

            elif linkage_matrix.loc[cur_col, cur_row] != np.inf:
                print "already calculated similarity for {}, {}".format(cur_row, cur_col)
                linkage_matrix.loc[cur_row, cur_col] = linkage_matrix.loc[cur_col, cur_row]

            else:
                # single_similarity(linkage_matrix, dict_terms_by_traits, cur_row, cur_col)
                weighted_maximal_similarity(linkage_matrix, dict_terms_by_traits, cur_row, cur_col)
                # non_weighted_average_similarity(linkage_matrix, dict_terms_by_traits, cur_row, cur_col)


    print linkage_matrix
    linkage_matrix.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "linkage_matrix.tsv"), sep='\t')

    fig, ax = plt.subplots(figsize=(20, 10))
    dn = hierarchy.dendrogram(hcl.linkage(squareform(linkage_matrix.values)), labels=linkage_matrix.index.values)

    # Apply the right color to each label
    my_palette = plt.cm.get_cmap("jet", 6)
    patches = [Line2D([0], [0], marker='o', color='gray', label=l,
                      markerfacecolor=my_palette(i)) for i, l in enumerate(["metabolic", "brain", "inflammatory", "eyes", "others", "cancer"])]
    ax.legend(handles=patches)

    ax = plt.gca()
    xlbls = ax.get_xmajorticklabels()
    num = -1
    for lbl in xlbls:
        lbl.set_color(my_palette(my_color[lbl.get_text()[5:]]))



    plt.tight_layout()
    plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR,"traits_dendogram.png"))

    results.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "terms_by_traits.tsv"), sep='\t')
    # plot_sumary_scatter(results, all_terms_names)

lbl.get_text()