import sys
sys.path.insert(0, "../")
import pandas as pd
import numpy as np
import constants
import os
from multiprocessing import Process
import constants 


# cancer_codes=["KIRC", "KIRP", "LUSC", "LUAD", "COAD", "BRCA", "STAD", "LIHC", "READ", "PRAD", "BLCA", "HNSC", "THCA", "UCEC", "OV", "PAAD"]
# country_codes=["US" for a  in cancer_codes]

cancer_codes=["PACA", "LIRI"] # ,  "BRCA"] # ["PRAD", "LICA", "RECA", "OV", "PACA", "PACA", "PAEN"]

country_codes=["CA", "JP"] # , "KR"] # ["FR", "FR", "EU", "AU", "AU", "CA", "AU"]

  
def format_matrix(cancer_code, country_code):

    dataset_name="ICGC_{}_{}".format(cancer_code, country_code)
    file_name="exp_seq.{}-{}.tsv".format(cancer_code, country_code)

    gene_id_field="gene_id"
    norm_value_field="total_read_count"
    gene_value_field="normalized_read_count"
    patient_id_field="icgc_sample_id"


    df_icgc=pd.read_csv(os.path.join(constants.DATASETS_DIR, dataset_name, "data", file_name), sep='\t')
    df_matrix=pd.DataFrame()

    gene_ids=df_icgc[gene_id_field]
    gene_values=df_icgc[gene_value_field]
    # df_icgc[norm_value_field]=1.0 
    norm_values=df_icgc[norm_value_field]
    patient_ids=df_icgc[patient_id_field]


    for g, v, n, p in zip(gene_ids.values, gene_values.values, norm_values.values, patient_ids.values):
       if np.isnan(n): continue
       df_matrix.loc[g, p]= v # (v/float(n)) * 1e6 
    
    df_matrix=df_matrix.apply(lambda x: np.log2(1+x)) 
    df_matrix.to_csv(os.path.join(constants.DATASETS_DIR, dataset_name, "data", "{}.htseq_rsem.tsv".format(dataset_name)), sep='\t', index_label="id")

    print "done {}, {}".format(cancer_code, country_code)


for cancer_code, country_code in zip(cancer_codes, country_codes):
    print "start {}, {}".format(cancer_code, country_code)
    p = Process(target=format_matrix, args=(cancer_code, country_code))
    p.start()  


