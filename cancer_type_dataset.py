import numpy as np
import pandas as pd
from torch.utils.data import Dataset
import torch
import infra
import constants
import simplejson as json
from utils.param_builder import build_gdc_params

# CANCER_TYPES=['LUSC', 'LUAD' , 'MESO', 'HNSC', 'BRCA', 'PRAD', 'SKCM', 'UVM', 'KIRP', 'KICH', 'KIRC', 'GBM', 'LGG', 'STAD', 'PAAD']
CANCER_TYPES=["KIRC", "KIRP", "KICH", "LUSC", "LUAD", "COAD", "BRCA", "STAD", "LIHC", "READ", "PRAD", "BLCA", "ESCA", "HNSC", "THCA", "UCEC"]
META_GROUPS = ["groups/temp.json", "groups/temp.json", "groups/temp.json","groups/temp.json", "groups/temp.json",  "groups/temp.json", "groups/temp.json",  "groups/temp.json", "groups/temp.json",  "groups/temp.json", "groups/temp.json",  "groups/temp.json", "groups/temp.json",  "groups/temp.json", "groups/temp.json", "groups/temp.json"]

n_input_layer=2000


class CancerTypesDataset(Dataset):

    def __init__(self, dataset_names, meta_groups_files, metagroups_names):

        self.labels = np.array([])
        self.labels_unique=np.array([])
        self.samples = pd.DataFrame()
        label_counter=0
        for dataset_name, meta_groups_file, metagroups_name in zip(dataset_names, meta_groups_files, metagroups_names):
            constants.update_dirs(DATASET_NAME_u=dataset_name)
            meta_groups = [json.load(file(meta_groups_file))]
            data_normalizaton = "fpkm"
            gene_expression_file_name, phenotype_file_name, survival_file_name, mutation_file_name, mirna_file_name, pval_preprocessing_file_name = build_gdc_params(
                dataset=dataset_name, data_normalizaton=data_normalizaton)

            tested_gene_list_file_name = "protein_coding_long.txt"  #
            total_gene_list_file_name = "protein_coding_long.txt"  #
            filter_expression = None
            data = infra.load_integrated_ge_data(tested_gene_list_file_name=tested_gene_list_file_name,
                                                 total_gene_list_file_name=total_gene_list_file_name,
                                                 gene_expression_file_name=gene_expression_file_name,
                                                 phenotype_file_name=phenotype_file_name,
                                                 survival_file_name=survival_file_name,
                                                 var_th_index=None, meta_groups=meta_groups,
                                                 filter_expression=filter_expression)

            gene_expression_top_var, gene_expression_top_var_headers_rows, gene_expression_top_var_headers_columns, labels_assignment, survival_dataset = data

            labels_assignment=np.array(labels_assignment)[0]
            for cur_label in np.unique(labels_assignment):
                cur_label_name=[cur["_name"] for cur in meta_groups[0] if "_label" in cur and int(cur["_label"])==cur_label]
                cur_label_name = "{}, {}".format(metagroups_name, cur_label_name[0] if len(cur_label_name) > 0 else "unknown")
                print cur_label_name
                df_new = pd.DataFrame(data=gene_expression_top_var[labels_assignment==cur_label], index=gene_expression_top_var_headers_rows[labels_assignment==cur_label],
                                      columns=gene_expression_top_var_headers_columns)
                self.samples = pd.concat([self.samples, df_new], axis=0)
                self.labels = np.append(self.labels, [cur_label_name for x in range(len(df_new.index))])
                label_counter+=1
                self.labels_unique = np.append(self.labels_unique, [cur_label_name])

        var_th_index =n_input_layer-1
        if var_th_index is not None:
            print "filtering top vars"
            gene_expression_top_var, gene_expression_top_var_headers_rows, gene_expression_top_var_headers_columns = infra.filter_top_var_genes(
                self.samples.values.T, self.samples.index.values, self.samples.columns, var_th_index)
            self.samples = pd.DataFrame(data=gene_expression_top_var, index=gene_expression_top_var_headers_rows,
                                        columns=gene_expression_top_var_headers_columns).T

        self.samples = self.samples.dropna(axis=1) / self.samples.dropna(axis=1).max()

    def __len__(self):
        return self.samples.shape[0]

    def __getitem__(self, idx):
        if type(idx)==torch.Tensor:
            idx=idx.item()

        result = torch.tensor(self.samples.iloc[idx, :].values, dtype=torch.float), torch.tensor(
            (self.labels_unique == list(self.labels)[idx]).astype(np.int), dtype=torch.long)
        # print idx, len(self.samples.index), result[1]
        return result

    def get_labels_unique(self):

        return self.labels_unique
