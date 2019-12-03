import numpy as np
import pandas as pd
from torch.utils.data import Dataset
import torch
import infra
import constants
import simplejson as json
from utils.param_builder import build_gdc_params
import os
# CANCER_TYPES=['LUSC', 'LUAD' , 'MESO', 'HNSC', 'BRCA', 'PRAD', 'SKCM', 'UVM', 'KIRP', 'KICH', 'KIRC', 'GBM', 'LGG', 'STAD', 'PAAD']
study="sc_melanoma"
DATASETS= ["CDKexp.UACC257_tpm.txt"] # "CDKexp.A2058_tpm.txt", "CDKexp.IGR37_tpm.txt",

# META_GROUPS =  ["groups/stages.json", "groups/stages.json", "groups/stages.json", "groups/stages.json", "groups/stages.json", "groups/stages.json", "groups/stages.json", "groups/stages.json", "groups/stages.json", "groups/stages.json",
#                 "groups/stages.json", "groups/stages.json", "groups/stages.json", "groups/stages.json", "groups/stages.json", "groups/stages.json"] # ["groups/temp.json", "groups/temp.json" , "groups/temp.json","groups/temp.json", "groups/temp.json",  "groups/temp.json", "groups/temp.json",  "groups/temp.json", "groups/temp.json",  "groups/temp.json", "groups/temp.json",  "groups/temp.json", "groups/temp.json",  "groups/temp.json", "groups/temp.json", "groups/temp.json"]

META_GROUPS =  ["groups/temp.json"]

n_input_layer=2000


class SingleCellDataset(Dataset):
    
    def __init__(self, dataset_names, meta_groups_files, metagroups_names):

        self.labels = np.array([])
        self.labels_unique=np.array([])
        self.samples = pd.DataFrame()
        self.survival= pd.DataFrame()
        label_counter=0
        for dataset in dataset_names:

            df_ge=pd.read_csv(os.path.join(constants.DATASETS_DIR, study, "data", dataset), sep='\t', index_col=0)
            self.samples = pd.concat([self.samples, df_ge], axis=1)

        # print "filtering top vars"
        # print self.samples.shape
        tested_gene_expression, tested_gene_expression_headers_rows, tested_gene_expression_headers_columns = infra.filter_top_var_genes(
            self.samples.values, self.samples.columns.values,
            self.samples.index.values, n_input_layer-1)
        self.samples = pd.DataFrame(data=tested_gene_expression, index=tested_gene_expression_headers_rows,
                                    columns=tested_gene_expression_headers_columns).T
        # print self.samples.shape
        self.samples = self.samples.dropna(axis=1)/ self.samples.dropna(axis=1).max()

    def __len__(self):
        return self.samples.shape[0]

    def __getitem__(self, idx):
        if type(idx)==torch.Tensor:
            idx=idx.item()

        result = torch.tensor(self.samples.iloc[idx, :].values, dtype=torch.float), torch.tensor([0], dtype=torch.float)
        # print "get_item shape:", result[0].shape
        return result


    def get_labels_unique(self):

        return [0]


