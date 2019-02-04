import numpy as np
import pandas as pd
from torch.utils.data import Dataset
import torch
import infra
import constants
import simplejson as json
from utils.param_builder import build_gdc_params
class CancerTypesDataset(Dataset):

    def __init__(self, csv_files, labels=None):


        self.labels=np.array([])
        self.samples=pd.DataFrame()
        for csv_file, csv_label in zip(csv_files, labels):
            constants.update_dirs(DATASET_NAME_u=csv_label)
            meta_groups = [json.load(file("groups/temp.json"))]
            data_normalizaton = "fpkm"
            gene_expression_file_name, phenotype_file_name, survival_file_name, mutation_file_name, mirna_file_name, pval_preprocessing_file_name = build_gdc_params(
                dataset=csv_label, data_normalizaton=data_normalizaton)

            tested_gene_list_file_name = "protein_coding_long.txt"  # "mir_total.txt" #
            total_gene_list_file_name = "protein_coding_long.txt"  # "mir_total.txt" #
            filter_expression = None
            data = infra.load_integrated_ge_data(tested_gene_list_file_name=tested_gene_list_file_name,
                                           total_gene_list_file_name=total_gene_list_file_name,
                                           gene_expression_file_name=gene_expression_file_name,
                                           phenotype_file_name=phenotype_file_name,
                                           survival_file_name=survival_file_name,
                                           var_th_index=None, meta_groups=meta_groups,
                                           filter_expression=filter_expression  )

            gene_expression_top_var, gene_expression_top_var_headers_rows, gene_expression_top_var_headers_columns, labels_assignment, survival_dataset = data


            df_new=pd.DataFrame(data=gene_expression_top_var, index=gene_expression_top_var_headers_rows, columns=gene_expression_top_var_headers_columns)
            self.samples = pd.concat([self.samples, df_new], axis=0)
            self.labels = np.append(self.labels, [csv_label for x in range(df_new.shape[0])])
             
        var_th_index = 99
        if var_th_index is not None:
            print "filtering top vars"
            gene_expression_top_var, gene_expression_top_var_headers_rows, gene_expression_top_var_headers_columns = infra.filter_top_var_genes(
                self.samples.values.T, self.samples.index.values, self.samples.columns, var_th_index)
            self.samples = pd.DataFrame(data=gene_expression_top_var, index=gene_expression_top_var_headers_rows,
                                  columns=gene_expression_top_var_headers_columns).T

        
        self.samples=self.samples.dropna()/self.samples.dropna().max()



    def __len__(self):
        return self.samples.shape[0]

    def __getitem__(self, idx):
        labels=np.array(['LUSC', 'LUAD', 'MESO', 'HNSC', 'BRCA', 'PRAD', 'SKCM', 'UVM', 'KIRP', 'KICH', 'KIRC', 'GBM', 'LGG', 'STAD',  'PAAD'])
        result=torch.tensor(self.samples.iloc[idx,:].values, dtype=torch.float), torch.tensor((labels==list(self.labels)[idx]).astype(np.int), dtype=torch.long)
        # print idx, len(self.samples.index), result[1]
        return result
