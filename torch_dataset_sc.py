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
DATASETS= ["sc_melanoma"]

# META_GROUPS =  ["groups/stages.json", "groups/stages.json", "groups/stages.json", "groups/stages.json", "groups/stages.json", "groups/stages.json", "groups/stages.json", "groups/stages.json", "groups/stages.json", "groups/stages.json",
#                 "groups/stages.json", "groups/stages.json", "groups/stages.json", "groups/stages.json", "groups/stages.json", "groups/stages.json"] # ["groups/temp.json", "groups/temp.json" , "groups/temp.json","groups/temp.json", "groups/temp.json",  "groups/temp.json", "groups/temp.json",  "groups/temp.json", "groups/temp.json",  "groups/temp.json", "groups/temp.json",  "groups/temp.json", "groups/temp.json",  "groups/temp.json", "groups/temp.json", "groups/temp.json"]

META_GROUPS =  ["groups/temp.json"]

n_input_layer=2000


class SingleCellDataset(Dataset):

    def __init__(self, dataset_names=DATASETS, meta_groups_files=META_GROUPS, metagroups_names=DATASETS):

        self.labels = np.array([])
        self.labels_unique=np.array([])
        self.samples = pd.DataFrame()
        self.survival= pd.DataFrame()
        label_counter=0
        for dataset_name, meta_groups_file, metagroups_name in zip(dataset_names, meta_groups_files, metagroups_names):

            df_ge=pd.read_csv(os.path.join(constants.DATASETS_DIR, dataset_name, "ge.tsv"), sep='\t', index_col=0).T
            self.samples = pd.concat([self.samples, df_ge], axis=0)


        self.samples = self.samples.dropna(axis=1)

    def __len__(self):
        return self.samples.shape[0]

    def __getitem__(self, idx):
        if type(idx)==torch.Tensor:
            idx=idx.item()

        result = torch.tensor(self.samples.iloc[idx, :].values, dtype=torch.float), 0
        return result


    def get_labels_unique(self):

        return [0]
