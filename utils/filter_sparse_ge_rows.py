import constants
import infra
import pandas as pd
import numpy as np
import os

constants.update_dirs(DATASET_NAME_u="PRAD_2")


#### PREPARE DICT ###

company="illumina" # agilent

# old_key_field = "AGILENT WholeGenome 4x44k v1 probe"
# new_field = "agilent_44_v1"

old_key_field = "ILLUMINA HumanWG 6 V3 probe"
new_key_field = "illumina_WG_v3"
old_value_field= "Gene stable ID"
new_value_field= "eid"


df_dict=pd.read_csv(os.path.join(constants.DICTIONARIES_DIR, "{}_ensembl_biomart.tsv".format(company)),sep='\t')
df_dict=df_dict.dropna()
df_dict.index=df_dict[old_key_field]
df_dict=df_dict.drop([old_key_field], axis=1)
df_dict=df_dict.rename(columns={old_value_field: new_value_field})
df_dict=df_dict[~df_dict.index.duplicated(keep='first')]
df_dict.to_csv(os.path.join(constants.DICTIONARIES_DIR, "{}_ensembl.tsv".format(company)), sep='\t', index_label=new_key_field)

#####################


#####PREPARE GE #######

df_dict=pd.read_csv(os.path.join(constants.DICTIONARIES_DIR, "{}_ensembl.tsv".format(company)),sep='\t', index_col=0)
raw_ge_matrix = pd.read_csv(os.path.join(constants.DATA_DIR, "ge_raw.tsv"),sep='\t',index_col=0)
raw_ge_matrix.index=df_dict.loc[raw_ge_matrix.index]["eid"]
raw_ge_matrix=raw_ge_matrix[~raw_ge_matrix.index.isna()]


counter=0
indices_to_drop=[]

for i,cur_row in raw_ge_matrix.iterrows():
    if cur_row.isna().sum() > 3:
        indices_to_drop.append(i)
        continue
    if cur_row.isna().sum() !=0:
        row_mean=cur_row.dropna().mean()
        raw_ge_matrix.loc[i,cur_row.isna()]=row_mean

raw_ge_matrix=raw_ge_matrix.drop(indices_to_drop, axis=0)
raw_ge_matrix=raw_ge_matrix[~raw_ge_matrix.index.duplicated(keep='first')]
raw_ge_matrix.to_csv(os.path.join(constants.DATA_DIR,"ge.tsv"),sep='\t')