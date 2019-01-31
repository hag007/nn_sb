import json
from matplotlib import style
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
import time

def load_gene_expression_by_label(tested_gene_file_name, expression_profile_file_name, phenotype_file_name, gene_filter_file_name=None, group_conditions=None):
    dist_data = []
    flatten_data = []
    if group_conditions is not None:
        group_conditions = json.load(file("../groups/{}.json".format(group_conditions)))
    gene_expression_dataset = np.array(load_gene_expression_profile_by_patients(tested_gene_file_name, expression_profile_file_name))

    groups = []
    if group_conditions is not None:
        patients_groups = divided_patient_ids_by_label(phenotype_file_name, groups=group_conditions)
        for cur_p_group in patients_groups:
            ds = gene_expression_dataset[np.in1d(gene_expression_dataset[:,0], cur_p_group),1:]
            groups.append(ds)
    else:
        groups = [gene_expression_dataset[1:][:,1:]]

    for cur in groups:
        dist_data.append(None)
    total_var=0
    for i, cur_group in enumerate(groups):
        if len(cur_group)==0 :continue

        for cur in cur_group:
            if dist_data[i] is None:
                dist_data[i] = np.array(cur)
            else:
                dist_data[i] = np.vstack((dist_data[i],cur))
        data_primary = dist_data[i].flatten().astype(np.float64)
        flatten_data.append(data_primary)
        sns.distplot(data_primary)
        total_var += np.var(data_primary)
        plt.savefig(os.path.join(constants.BASE_PROFILE, "output" ,"dist_comparison_{}_group_{}".format(expression_profile_file_name.split("_")[1].split('.')[0],i)))
        plt.cla()

    ax = plt.subplot(111)
    positions = np.arange(len(flatten_data)) + 1
    bp = ax.boxplot(flatten_data, positions=positions, showmeans=True)
    ax.set_title("genes_statistic_{}_{}_averaged var:{}".format(constants.CANCER_TYPE, expression_profile_file_name.split(".")[0], '%.3f' % (total_var/len(groups))))
    plt.savefig(os.path.join(constants.BASE_PROFILE, "output",
                             "genes_statistic_{}_{}_{}.png".format(constants.CANCER_TYPE,
                                                                   expression_profile_file_name.split(".")[0],
                                                                   time.time())))
    return flatten_data








ax = plt.subplot(111)
total_data = []
for cur_cancer_type in ["PAAD", "STAD", "UVM", "BRCA"]:
    dataset=cur_cancer_type
    data_normalizaton = "fpkm-uq"
    constants.update_dirs(CANCER_TYPE_u=dataset)
    gene_expression_file_name, phenotype_file_name, survival_file_name, mutation_file_name, mirna_file_name, pval_preprocessing_file_name = build_gdc_params(
        dataset=dataset, data_normalizaton=data_normalizaton)
    group_conditions = None # "cancer_types"
    flatten_data = load_gene_expression_by_label("protein_coding.txt", gene_expression_file_name, phenotype_file_name, group_conditions=group_conditions)
    total_data.append(flatten_data)
positions = np.arange(len(total_data)) + 1
plt.cla()
bp = ax.boxplot(total_data, positions=positions, showmeans=True)
ax.set_title("genes_statistic")
plt.savefig(os.path.join(constants.BASE_PROFILE, "output",
                         "genes_statistic_{}_{}.png".format(constants.CANCER_TYPE,
                                                               time.time())))




