import sys
sys.path.insert(0, "../")
import constants
import os
import shutil
import pandas as pd
import numpy as np 
cancer_codes=["PRAD", "BRCA", "RECA", "LIRI", "LICA"] 
country_codes=["FR", "KR", "EU", "JP", "FR"]

for cancer_code, country_code in zip(cancer_codes, country_codes):
    file_path=os.path.join(constants.DATASETS_DIR, "ICGC_{}_{}".format(cancer_code, country_code), "data", "ICGC_{}_{}.htseq_rsem.tsv".format(cancer_code, country_code))
    df_matrix=pd.read_csv(file_path, sep='\t', index_col=0)
    df_matrix=df_matrix.apply(lambda x: np.log2(x))
    df_matrix.to_csv(file_path,sep='\t')
   
   
