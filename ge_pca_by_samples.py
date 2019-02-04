from utils.ensembl2gene_symbol import e2g_convertor
from matplotlib import style

style.use("ggplot")
from scipy.stats import zscore
import scipy
import logging

sh = logging.StreamHandler()
logger = logging.getLogger("log")
logger.addHandler(sh)
from constants import *
from infra import *
from utils.pca import plot_pca_by_samples
from utils.pca import plot_pca
from utils.tsne import plot_tsne
from utils.param_builder import build_gdc_params, build_tcga_params

############################ () cluster and enrichment #############################


if __name__ == "__main__":

    for dataset in ["SKCM"]:
        constants.update_dirs(DATASET_NAME_u=dataset)
        meta_groups = [json.load(file("groups/temp.json"))]
        constants.update_dirs(CANCER_TYPE_u=dataset)
        data_normalizaton = "fpkm"
        gene_expression_file_name, phenotype_file_name, survival_file_name, mutation_file_name, mirna_file_name, pval_preprocessing_file_name = build_gdc_params(
            dataset=dataset, data_normalizaton=data_normalizaton)
        tested_gene_list_file_name = "protein_coding_long.txt"  # "mir_total.txt" #
        total_gene_list_file_name = "protein_coding_long.txt"  # "mir_total.txt" #
        var_th_index = 2000
        is_unsupervised = True
        start_k = 2
        end_k = 2
        # [{"gender.demographic": {"type": "string", "value": ["male"]}},
        # {"gender.demographic": {"type": "string", "value": ["female"]}}]
        filter_expression = None

        d_codes = []
        # for cur in meta_groups[0]:
        #     d_codes.append(cur["disease_code"]["value"][0])

        # , "person_neoplasm_cancer_status" : {"type": "string", "value" : ["WITH TUMOR"]}

        # for cur_tt in ["Primary Tumor"]:
        filter_expression = None  # [{"sample_type.samples": {"type": "string", "value": ["Primary Tumor","Metastatic"]},
        # "person_neoplasm_cancer_status" : {"type": "string", "value" : ["WITH TUMOR"]} ,
        # "disease_code": {"type": "string",
        #              "value": d_codes
        #              }
        # }]
        print "process {}".format(dataset)

        data = load_integrated_ge_data(tested_gene_list_file_name=tested_gene_list_file_name,
                                       total_gene_list_file_name=total_gene_list_file_name,
                                       gene_expression_file_name=gene_expression_file_name,
                                       phenotype_file_name=phenotype_file_name, survival_file_name=survival_file_name,
                                       var_th_index=var_th_index, meta_groups=meta_groups,
                                       filter_expression=filter_expression)

        gene_expression_top_var, gene_expression_top_var_headers_rows, gene_expression_top_var_headers_columns, labels_assignment, survival_dataset = data

        # data = load_integrated_mutation_data(mutation_file_name=mutation_file_name,
        #                                phenotype_file_name=phenotype_file_name, survival_file_name=survival_file_name,
        #                                var_th_index=var_th_index, meta_groups=meta_groups,
        #                                filter_expression=filter_expression)
        #
        #
        # if data is None:
        #     print "insufficient data"
        #     exit(1)
        #
        #
        # mutation_dataset, mutations_headers_rows, mutations_headers_columns, labels_assignment, survival_dataset, phenotype_heatmap = data
        #
        # min_ratio = 0.001
        # included_mutation_gene_list = None
        # excluded_mutation_gene_list = None
        #
        #
        # all_patients = np.unique(mutations_headers_rows).flatten()
        # all_mutated_genes = np.unique(mutation_dataset[:, 0]).flatten()
        # # mis_mutated_genes = np.unique(mu_data[np.where(np.core.defchararray.find(mu_data[1:, 8], "missense")!=-1), 1]).flatten()
        #
        # all_mutated_vectors = np.zeros((len(all_patients), len(all_mutated_genes)))
        # # mis_mutated_vectors = np.array([[0 for y in mis_mutated_genes] for x in range(len(all_patients))])
        #
        # print "build vectors from {} entries".format(len(mutation_dataset))
        #
        # stopwatch = Stopwatch()
        # stopwatch.start()
        # a = list(all_patients)
        # b = list(all_mutated_genes)
        # for i, x in enumerate(mutation_dataset):
        #     all_mutated_vectors[a.index(mutations_headers_rows[i])][b.index(x[0])] += 1
        # print stopwatch.stop("end mut")
        # all_mutated_vectors[all_mutated_vectors > 5] = 5
        #
        # if included_mutation_gene_list is not None:
        #     included_mutation_gene = load_gene_list(included_mutation_gene_list)
        #     all_mutated_vectors = all_mutated_vectors[:, np.in1d(all_mutated_genes, included_mutation_gene)]
        #     all_mutated_genes = all_mutated_genes[np.in1d(all_mutated_genes, included_mutation_gene)]
        #
        # if excluded_mutation_gene_list is not None:
        #     excluded_mutation_gene = load_gene_list(excluded_mutation_gene_list)
        #     for cur in excluded_mutation_gene:
        #         all_mutated_vectors = all_mutated_vectors[:, all_mutated_genes != cur]
        #         all_mutated_genes = all_mutated_genes[all_mutated_genes != cur]
        #
        # all_mutated_vectors[all_mutated_vectors > 5] = 5
        # all_mutated_genes = all_mutated_genes[
        #     (all_mutated_vectors != 0).sum(axis=0) > np.shape(all_mutated_vectors)[0] * min_ratio]
        # all_mutated_vectors = all_mutated_vectors[:,
        #                       (all_mutated_vectors != 0).sum(axis=0) > np.shape(all_mutated_vectors)[0] * min_ratio]
        # print "all_mutated_vectors after filter sparse: {}".format(np.shape(all_mutated_vectors))
        #
        # if np.size(all_mutated_genes) == 0:
        #     exit(0)
        #
        # pheno=divided_patient_ids_by_label(phenotype_file_name, groups=meta_groups[0])
        # labels_assignment=[np.array([1 if x in pheno[0] else 2 if x in pheno[1]  else 3 if x in pheno[2]  else 4 if x in pheno[3] else 0 for x in all_patients])]

        plot_pca_by_samples(gene_expression_top_var, labels_assignment, meta_groups, n_components=2)
        # from sklearn.decomposition import PCA
        # n_components=2
        # import matplotlib.pyplot as plt
        # import time
        #
        # X = np.array(gene_expression_top_var)
        # y = np.array([0 for x in gene_expression_top_var])
        # pca = PCA(n_components=n_components)
        # pca.fit_transform(X)
        #
        # fig = plt.figure(1, figsize=(20, 20))
        # plt.clf()
        #
        # if n_components == 2:
        #     ax = fig.add_subplot(111)
        #     ax.scatter(X[:, 0], X[:, 1], c=y, cmap='jet')
        # plt.savefig(os.path.join(constants.BASE_PROFILE, "output", "PCA_by_samples_{}_{}_{}_{}.png").format(
        #     constants.CANCER_TYPE,
        #     tested_gene_list_file_name.split(".")[0] if tested_gene_list_file_name is not None else "none",
        #     n_components, time.time()))


        plot_tsne(gene_expression_top_var, labels_assignment, meta_groups, n_components=2)
