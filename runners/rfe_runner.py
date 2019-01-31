from utils.omic_svm import DISTANCE, LOGISTIC_REGRESSION, apply_svm, svm_linear_default
from utils.rfe import RANDOMIZED, REVERSED, NORMAL, RFE, print_to_excel
from utils.param_builder import build_gdc_params
import constants
import json
from utils.pca import plot_detailed_pca
from utils.groups_generator import generate_random_set
from infra import load_integrated_ge_data
import numpy as np
import matplotlib.pyplot as plt
import os

dataset = "BRCA"
constants.update_dirs(DATASET_NAME_u=dataset)


data_normalizaton = "counts_normalized_by_genes_standardization"
gene_expression_file_name, phenotype_file_name, survival_file_name, mutation_file_name, mirna_file_name, pval_preprocessing_file_name = build_gdc_params(dataset=dataset, data_normalizaton=data_normalizaton)
constants.PHENOTYPE_FORMAT = "TCGA"
phenotype_file_name = 'BRCA_clinicalMatrix'

cur_json = "brca_pam53"
meta_groups = None
filter_expression = None
meta_groups = [json.load(file("../groups/{}.json".format(cur_json)))]
filter_expression = json.load(file("../filters/{}.json".format(cur_json)))
var_th_index = None

gene_list_file_name = "dip_bg.txt"
rounds=1
rank_method = DISTANCE
recursion_number_of_steps=20
recursion_step_size= 50
permutation=RANDOMIZED

num_of_genes=0
prs = []
reduced_prs=[]

for cur in range(recursion_number_of_steps):
    num_of_genes+=recursion_step_size
    random_file_name=generate_random_set(random_size=num_of_genes, meta_gene_set="dip_bg.txt")
    X, y, pr, algo_roc = plot_detailed_pca(tested_gene_list_file_name=random_file_name,
                                                total_gene_list_file_name="dip_bg.txt",
                                                gene_expression_file_name=gene_expression_file_name,
                                                phenotype_file_name=phenotype_file_name,
                                                survival_file_name=survival_file_name,
                                                filter_expression=filter_expression,
                                                meta_groups=meta_groups,
                                                var_th_index=var_th_index,
                                                algo_name="")

    reduced_prs.append(pr)

    data = load_integrated_ge_data(tested_gene_list_file_name=random_file_name,
                                   total_gene_list_file_name="dip_bg.txt",
                                   gene_expression_file_name=gene_expression_file_name,
                                   phenotype_file_name=phenotype_file_name, survival_file_name=survival_file_name,
                                   var_th_index=var_th_index, meta_groups=meta_groups,
                                   filter_expression=filter_expression)
    if data is None:
        print "insufficient data"
    values, h_rows, h_cols, labels_assignment, survival_dataset = data

    X = values
    feature_ids = np.array([x.split('.')[0] for x in h_cols])
    labels = labels_assignment[0]
    clf_method = svm_linear_default({'C': [10], 'kernel': ['linear']})
    y=[a-1 for a in y]
    pr, roc, clf = apply_svm(clf_method, X, y, X, y, DISTANCE)
    prs.append(pr)

print "prs: {}".format(" ".join([str(a) for a in prs]))
print "reduced prs: {}".format(" ".join([str(a) for a in reduced_prs]))
plt.cla()
ax = plt.subplot()
ax.plot(np.arange(recursion_number_of_steps)*recursion_step_size+recursion_step_size, prs, label="without reduction", color='blue')
ax.plot(np.arange(recursion_number_of_steps)*recursion_step_size+recursion_step_size, reduced_prs, label="with reduction", color='red')
ax.set_xlabel("number of random genes")
ax.set_ylabel("PR score")
plt.legend()
plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR,"PCA_VS_NORMAL_DISCRIMINATION.png"))


