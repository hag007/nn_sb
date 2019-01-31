import json
import os

import numpy as np
import pandas as pd
from scipy.stats import zscore

import constants
import infra
from datasets_multithread_runner import run_dataset
from utils.cache import clear_cache
from utils.groups_generator import generate_random_set
from utils.kmeans import kmeanssample
from utils.param_builder import build_params
from utils.pca import plot_detailed_pca

RAND_TIMES = 2
KMEANS_TIMES = 2
TOP_SIG= True

data_normalizaton = "counts"
deg_method = constants.DEG_EDGER
deg_file_name = "deg_{}.tsv".format(deg_method)
use_algo_cache=False
algos = ["netbox", "hotnet2", "bionet", "jactivemodules_greedy",
                "jactivemodules_sa", "keypathwayminer_INES_GREEDY"] # "matisse", "reactomefi"
# algos=['jactivemodules_sa'] # jactivemodules_greedy
plot_svm=False

def main(dataset, cur_json, ds_types = "GDC"):
    constants.update_dirs(DATASET_NAME_u=dataset)
    meta_groups = None
    filter_expression = None
    meta_groups = [json.load(file("../filters/{}.json".format(cur_json)))]
    filter_expression = json.load(file("../filters/{}.json".format(cur_json)))

    gene_expression_file_name, phenotype_file_name, survival_file_name, mutation_file_name, mirna_file_name, pval_preprocessing_file_name = build_params(
        type=ds_types, dataset=constants.DATASET_NAME, data_normalizaton=data_normalizaton)
    gene_expression_normalized_file_name =  "ge_normalized.tsv" # gene_expression_file_name
    survival_file_name="none"
    tested_gene_expression, h_rows, h_cols, labels_assignment, survival_dataset = infra.load_integrated_ge_data("dip_bg.txt", "dip_bg.txt", gene_expression_file_name,
                            survival_file_name, phenotype_file_name, gene_filter_file_name=None, filter_expression=filter_expression,
                            meta_groups=meta_groups, var_th_index=None)
    file(os.path.join(constants.DATA_DIR, "classes.tsv"), 'w+').write('\t'.join([str(x) for x in labels_assignment[0]]))
    h_cols= [x.split('.')[0] for x in h_cols]
    df_data = pd.DataFrame(index=h_rows, columns=h_cols,data=tested_gene_expression).T
    df_data.to_csv(os.path.join(constants.DATA_DIR, "ge.tsv"), index_label="eid", sep="\t")
    var_th_index = None

    if not use_algo_cache:
        run_dataset(dataset, score_method=deg_method, algos=algos)
    # exit(0)
    gene_list_file_names = []
    prs = pd.DataFrame(columns=['algo', 'algo_pr', 'kmean_pr_avg', 'kmean_pr_std',
               'algo_kmean_ratio', 'top_sig_pr', "algo_top_sig_ratio", 'rand_pr_avg', 'rand_pr_std',
                "algo_rand_ratio", 'algo_pr_rank_from_rand', "num of modules", 'num_of_genes'])

    for cur_algo in algos:
        print "about to start running {}".format(cur_algo)
        if not os.path.exists(os.path.join(constants.OUTPUT_GLOBAL_DIR,'pca', dataset, cur_algo)):
            os.makedirs(os.path.join(constants.OUTPUT_GLOBAL_DIR,'pca', dataset, cur_algo))
        algo_output = json.loads(file(os.path.join(constants.OUTPUT_DIR,"{}_client_output.txt".format(cur_algo))).read().split("\n")[1])
        module_i=0
        algo_pvals = []
        df_mean = pd.DataFrame()
        gene_2_module = {}
        num_of_genes=0
        algo_genes_flatted = []

        while True:
            module_genes_flatted = [x['eid'] for x in algo_output if module_i in x['modules']]
            algo_genes_flatted += module_genes_flatted
            num_of_genes += len(module_genes_flatted)
            print "# of genes in module {} : {}".format(module_i,len(module_genes_flatted))
            for cur in module_genes_flatted:
                gene_2_module[cur] = module_i

            algo_genes_flatted += module_genes_flatted
            if len(module_genes_flatted)==0 and module_i>0: break
            if len(module_genes_flatted) < 4 or sum(df_data.index.isin(module_genes_flatted))==0:
                module_i+=1
                continue
            gene_list_file_names.append(os.path.join(constants.LIST_DIR, cur_algo + ".txt"))
            file(gene_list_file_names[-1],'w+').write("\n".join(module_genes_flatted))
            df_mean = pd.concat((df_mean, pd.DataFrame(zscore(df_data[df_data.index.isin(module_genes_flatted)],axis=1).mean(axis=0).reshape(1,len(df_data.columns)), columns=df_data.columns)))

            module_i+=1

        module_i = len(df_mean.index)

        if module_i < 2:
            print "not enough modules. retrieved {}".format(module_i)
            continue

        mean_file_name = os.path.join(constants.DATA_DIR, "mean.tsv")
        df_mean.index = np.arange(df_mean.shape[0])
        df_mean.to_csv(mean_file_name, sep="\t", index_label="eid")
        index_file_name = os.path.join(constants.LIST_DIR, "{}_indices.txt".format(cur_algo))
        file(index_file_name,'w+').write('\n'.join([str(x) for x in df_mean.index.values]))

        algo_genes_file_name=os.path.join(constants.LIST_DIR, "{}_all_genes.txt".format(cur_algo))
        file(algo_genes_file_name,'w+').write('\n'.join([str(x) for x in algo_genes_flatted]))

        df_algo_gene_matrix = df_data[df_data.index.isin(algo_genes_flatted)]


        results = plot_detailed_pca(tested_gene_list_file_name="{}_indices.txt".format(cur_algo),
                                                     total_gene_list_file_name="protein_coding.txt",
                                                     gene_expression_file_name=mean_file_name,
                                                     phenotype_file_name=phenotype_file_name,
                                                     survival_file_name=survival_file_name,
                                                     filter_expression=filter_expression,
                                                     meta_groups=meta_groups,
                                                     var_th_index=var_th_index,
                                                     algo_name=cur_algo,
                                                 plot_svm=plot_svm)

        if results is None:
            continue
        X, y, algo_pr, algo_roc = results

        print "results for mean: {}".format(algo_pr)
        algo_bg_pr_mean = 0
        algo_bg_pr_std = 0
        if KMEANS_TIMES > 1:
            all_algo_bg_pr = []
            for kmean_i in range(KMEANS_TIMES):
                _1, clusters,_2 = kmeanssample(X=df_algo_gene_matrix.values, k=module_i, metric="euclidean")

                bg_modules=[df_algo_gene_matrix.index.values[clusters == cur_i] for cur_i in range(module_i)]
                df_mean_bg = pd.DataFrame()
                for i, module in enumerate(bg_modules):
                    print "# of genes in background module {} : {}".format(i, len(module))

                    gene_list_file_names.append(os.path.join(constants.LIST_DIR, cur_algo + "_bg.txt"))
                    file(gene_list_file_names[-1],'w+').write("\n".join(module_genes_flatted))
                    df_mean_bg = pd.concat((df_mean_bg, pd.DataFrame(zscore(df_data[df_data.index.isin(module)],axis=1).mean(axis=0).reshape(1,len(df_data.columns)), columns=df_data.columns)))

                bg_mean_file_name = os.path.join(constants.DATA_DIR, "mean_bg.tsv")
                df_mean_bg.index = np.arange(df_mean_bg.shape[0])
                df_mean_bg.to_csv(bg_mean_file_name, sep="\t", index_label="eid")
                index_file_name = os.path.join(constants.LIST_DIR, "{}_bg_indices.txt".format(cur_algo))
                file(index_file_name,'w+').write('\n'.join([str(x) for x in df_mean_bg.index.values]))

                all_genes_file_name = os.path.join(constants.LIST_DIR, "{}_all_genes.txt".format(cur_algo))
                file(all_genes_file_name, 'w+').write('\n'.join(algo_genes_flatted))


                bg_genes = pd.read_csv(os.path.join(constants.CACHE_DIR, deg_file_name.format(deg_method)), sep='\t',
                                            index_col=0).index.values[:len(df_mean.index)]
                bg_genes_file_name = os.path.join(constants.LIST_DIR,
                                                       "{}_{}_bg_genes.txt".format(cur_algo, deg_method))
                file(bg_genes_file_name, 'w+').write('\n'.join([x.split('.')[0] for x in bg_genes]))

                X,y,algo_bg_pr,algo_bg_roc = plot_detailed_pca(tested_gene_list_file_name="{}_bg_indices.txt".format(cur_algo),
                                                             total_gene_list_file_name="protein_coding.txt",
                                                             gene_expression_file_name=bg_mean_file_name,
                                                             phenotype_file_name=phenotype_file_name,
                                                             survival_file_name=survival_file_name,
                                                             filter_expression=filter_expression,
                                                             meta_groups=meta_groups,
                                                             var_th_index=var_th_index,
                                                             algo_name=cur_algo, plot_svm=plot_svm)

                print "results for mean: {}".format(algo_bg_pr)
                all_algo_bg_pr.append(algo_bg_pr)
            all_algo_bg_pr = np.array(all_algo_bg_pr)
            algo_bg_pr_mean = all_algo_bg_pr.mean()
            algo_bg_pr_std = all_algo_bg_pr.mean()

        top_sig_pr=0
        if TOP_SIG:
            top_sig_genes = pd.read_csv(os.path.join(constants.CACHE_DIR, deg_file_name),sep='\t', index_col=0)
            top_sig_genes=top_sig_genes.index.values[:len(top_sig_genes.index)/200] # len(df_mean.index)
            top_sig_genes_file_name = os.path.join(constants.LIST_DIR, "{}_{}_top_sig_genes.txt".format(cur_algo,deg_method))
            file(top_sig_genes_file_name, 'w+').write('\n'.join([x.split('.')[0] for x in top_sig_genes]))

            X, y, top_sig_pr, top_sig_roc = plot_detailed_pca(tested_gene_list_file_name=os.path.basename(top_sig_genes_file_name),
                                     total_gene_list_file_name="protein_coding.txt",
                                     gene_expression_file_name=gene_expression_normalized_file_name,
                                     phenotype_file_name=phenotype_file_name,
                                     survival_file_name=survival_file_name,
                                     filter_expression=filter_expression,
                                     meta_groups=meta_groups,
                                     var_th_index=var_th_index,
                                     algo_name=cur_algo, plot_svm=plot_svm)
            print "results for top: {}".format(top_sig_pr)


        rand_prs_mean = 0
        rand_prs_std = 0
        rand_prs = []
        trials=0
        if RAND_TIMES > 1:
            while trials < RAND_TIMES:
                random_set_file_name = generate_random_set(random_size=len(df_mean.index), # df_mean.index
                                                           meta_gene_set="dip_bg.txt".format(cur_algo))
                print "running {} iteration for {} random bg with {} genes".format(trials, cur_algo, len(df_mean.index))
                results =plot_detailed_pca(tested_gene_list_file_name=random_set_file_name,
                                  total_gene_list_file_name="protein_coding.txt",
                                  gene_expression_file_name=gene_expression_normalized_file_name,
                                  phenotype_file_name=phenotype_file_name,
                                  survival_file_name=survival_file_name,
                                  filter_expression=filter_expression,
                                  meta_groups=meta_groups,
                                  var_th_index=var_th_index,
                                  feature_names=gene_2_module,
                                  algo_name=cur_algo, plot_svm=plot_svm)

                if results is None:
                    print "not enough genes retrieved. retry.."
                    continue
                X, y, rand_pr, rand_roc = results
                trials+=1
                rand_prs.append(rand_pr)
                print "results for random {}: {}".format(trials, rand_pr)
            rand_prs=np.array(rand_prs)
            rand_prs_mean = rand_prs.mean()
            rand_prs_std = rand_prs.std()

        row = {'algo': cur_algo, 'algo_pr': algo_pr, 'kmean_pr_avg': algo_bg_pr_mean, 'kmean_pr_std': algo_bg_pr_std,
               'algo_kmean_ratio': algo_pr/algo_bg_pr_mean, 'top_sig_pr': top_sig_pr, "algo_top_sig_ratio": algo_pr/top_sig_pr,
         'rand_pr_mean': rand_prs_mean, 'rand_pr_std': rand_prs_std, "algo_rand_ratio": algo_pr/rand_pr,
         'algo_pr_rank_from_rand': len([cur for cur in rand_prs if cur > algo_pr]), "num of modules": module_i,
         'num_of_genes': num_of_genes}
        row.update({'rand_pr' + str(i): v for i, v in enumerate(rand_prs)})
        prs=prs.append(row, ignore_index=True)


    prs=prs.set_index('algo')
    prs.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, 'pca', constants.DATASET_NAME, "pr_summary_{}_{}.tsv".format(dataset, cur_json)), sep='\t')
    print " algo pvals"
    print algo_pvals

if __name__=='__main__':
    # ds_types = "GDC"
    ds_types = "DEFAULT"

    # dataset_name = "UVM"
    # filter_jsons = ["uvm_eif1ax", "uvm_bap1", "uvm_sf3b1",
    #             "uvm_chromosome_abnormality", "uvm_eye_color", "uvm_histological_type",
    #             "uvm_tumor_shape_pathologic_clinical", "uvm_gna11", "uvm_gnaq", "uvm_extrascleral_extension", "new_tumor_event_after_initial_treatment"]
    # for cur_json in filter_jsons:
    #     print cur_json
    #     clear_cache(dataset_name)
    #     main(dataset_name, cur_json, ds_types = ds_types)
    # dataset_name = "SKCM"
    # filter_jsons = ["skcm_tss_tcga_0_1", "skcm_tss_tcga_1_2", "skcm_tss_tcga_2_3", "skcm_clark", "skcm_breslow", "new_tumor_event_after_initial_treatment"]
    # for cur_json in filter_jsons:
    #     print cur_json
    #     clear_cache(dataset_name)
    #     main(dataset_name, cur_json, ds_types = ds_types)
    # dataset_name = "BRCA"
    # ds_types = "TCGA"
    # constants.PHENOTYPE_FORMAT = "TCGA"
    # filter_jsons = ["brca_pam53"]
    # for cur_json in filter_jsons:
    #     print cur_json
    #     clear_cache(dataset_name)
    #     main(dataset_name, cur_json, ds_types = ds_types)


    dataset_name = "PRAD_2"
    filter_jsons = ["prad2_bone_metastatic_ar_noar", "prad2_bone_metastatic_crpc_other", "prad2_bone_metastatic_type"]
    for cur_json in filter_jsons:
        print cur_json
        clear_cache(dataset_name)
        main(dataset_name, cur_json, ds_types = ds_types)