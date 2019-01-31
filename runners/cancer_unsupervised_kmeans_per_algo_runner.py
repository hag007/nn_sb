import constants
import os
import infra
import json
from datasets_multithread_runner import run_dataset
from utils.omic_svm import prediction_by_gene_expression
from utils.omic_svm import DISTANCE, LOGISTIC_REGRESSION
from utils.param_builder import build_gdc_params
from utils.patients_clustering import find_clusters_and_survival
import shutil
import pandas as pd
from utils.groups_generator import generate_random_set
import numpy as np


RAND_TIMES = 10

def main(dataset="BRCA"):
    constants.update_dirs(DATASET_NAME_u=dataset)
    data_normalizaton = "counts_normalized_by_genes_standardization"
    cur_json = "brca_pam53"
    meta_groups = None
    filter_expression = None
    meta_groups = [json.load(file("../groups/{}.json".format(cur_json)))]
    filter_expression = json.load(file("../filters/{}.json".format(cur_json)))

    gene_expression_file_name, phenotype_file_name, survival_file_name, mutation_file_name, mirna_file_name, pval_preprocessing_file_name = build_gdc_params(
        dataset=dataset, data_normalizaton=data_normalizaton)
    phenotype_file_name='BRCA_clinicalMatrix'
    tested_gene_expression, h_rows, h_cols, labels_assignment, survival_dataset = infra.load_integrated_ge_data("dip_bg.txt", "dip_bg.txt", gene_expression_file_name,
                            survival_file_name, phenotype_file_name, gene_filter_file_name=None, filter_expression=filter_expression,
                            meta_groups=meta_groups, var_th_index=None)
    file(os.path.join(constants.DATA_DIR, "classes.tsv"),'w+').write('\t'.join([str(x) for x in labels_assignment[0]]))
    h_cols= [x.split('.')[0] for x in h_cols]
    df_data = pd.DataFrame(index=h_rows, columns=h_cols,data=tested_gene_expression).T
    df_data.to_csv(os.path.join(constants.DATA_DIR, 'ge.tsv'), index_label="eid", sep="\t")
    var_th_index = None
    start_k = 2
    end_k = 2

    # algos = ["matisse", "keypathwayminer_INES_GREEDY", "netbox", "hotnet2", "bionet", "jactivemodules_greedy",
    #          "jactivemodules_sa", "reactomefi"]
    algos=["netbox"]
    run_dataset(dataset, score_method=constants.DEG_EDGER)
    gene_list_file_names = []
    generate_plot = True
    clustering_algorithm = "correlation"
    for cur_algo in algos:
        algo_output = json.loads(file(os.path.join(constants.OUTPUT_DIR,"{}_client_output.txt".format(cur_algo))).read().split("\n")[1])
        i=0
        algo_pvals = []
        random_pvals = []
        df_mean = pd.DataFrame()
        all_algo_genes_flatted = []
        gene_2_module = {}
        while True:
            algo_genes_flatted = [x['eid'] for x in algo_output if i in x['modules']]
            for cur in algo_genes_flatted:
                gene_2_module[cur] = i
            all_algo_genes_flatted+=algo_genes_flatted
            if len(algo_genes_flatted)==0 and i>0: break
            if len(algo_genes_flatted) < 4:
                i+=1
                continue
            gene_list_file_names.append(os.path.join(constants.LIST_DIR, cur_algo + ".txt"))
            file(gene_list_file_names[-1],'w+').write("\n".join(algo_genes_flatted))
            df_mean = pd.concat((df_mean, df_data[df_data.index.isin(algo_genes_flatted)].mean().to_frame().T))

            i+=1

        all_genes_file_name = os.path.join(constants.LIST_DIR, "{}_all_genes.txt".format(cur_algo))
        file(all_genes_file_name, 'w+').write('\n'.join(all_algo_genes_flatted))

        mean_file_name=os.path.join(constants.DATA_DIR, "mean.tsv")
        df_mean.index = np.arange(df_mean.shape[0])
        df_mean.to_csv(mean_file_name, sep="\t", index_label="eid")
        index_file_name = os.path.join(constants.LIST_DIR, "{}_indices.txt".format(cur_algo))
        file(index_file_name,'w+').write('\n'.join([str(x) for x in df_mean.index.values]))

        algo_pvals.append(find_clusters_and_survival(tested_gene_list_file_name="{}_indices.txt".format(cur_algo),
                                                     total_gene_list_file_name="protein_coding.txt",
                                                     gene_expression_file_name=mean_file_name,
                                                     phenotype_file_name=phenotype_file_name,
                                                     survival_file_name=survival_file_name,
                                                     var_th_index=var_th_index, is_unsupervised=True, start_k=start_k,
                                                     end_k=end_k, filter_expression=filter_expression,
                                                     meta_groups=meta_groups, clustering_algorithm=clustering_algorithm,
                                                     plot=generate_plot))

        for cur in range(RAND_TIMES):
            random_set_file_name = generate_random_set(random_size=len(df_mean.index),
                                                       meta_gene_set="{}_all_genes.txt".format(cur_algo))

            random_pvals.append(find_clusters_and_survival(tested_gene_list_file_name=random_set_file_name,
                                                           total_gene_list_file_name="dip_bg.txt",
                                                           gene_expression_file_name=gene_expression_file_name,
                                                           phenotype_file_name=phenotype_file_name,
                                                           survival_file_name=survival_file_name,
                                                           var_th_index=var_th_index, is_unsupervised=True, start_k=start_k,
                                                           end_k=end_k, filter_expression=filter_expression,
                                                           meta_groups=meta_groups,
                                                           clustering_algorithm=clustering_algorithm, plot=generate_plot))





        print " algo pvals"
        print algo_pvals
        print "# above TH: {}".format(len([x for x in algo_pvals if any(y < 0.001 for y in x)]))
        print " random pvals"
        print random_pvals
        print "# above TH: {}".format(len([x for x in random_pvals if any(y < 0.001 for y in x)]))
        print "# of modules better over random: {}/{}" .format(len([x for x1,x2 in zip(algo_pvals, random_pvals) if x1[0] < x2[0]]), len(algo_pvals))
if __name__=='__main__':
    main()