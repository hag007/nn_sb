import os
import pandas as pd
import numpy as np
import constants
import subprocess
import shutil
from utils.ensembl2gene_symbol import g2e_convertor

output_file_format="{}.sum.genescores.txt"
pascal_base_folder="/home/hag007/Desktop/PASCAL"

file_base_names = ["20001_1001.assoc", "20001_1002.assoc", "20001_1019.assoc"]

for gwas_file_base_name in file_base_names:
    print "current file: {}".format(gwas_file_base_name)
    gwas_file_name_format = gwas_file_base_name + "{}.tsv"
    gwas_file_name = gwas_file_name_format.format("")
    output_file_name=output_file_format.format(gwas_file_base_name)
    df=pd.read_csv(os.path.join(constants.RAW_DIR,"gwas", gwas_file_name), sep='\t')
    df=df[["rsid", "pval"]]
    df.loc[df['pval']<=0]=np.nan
    df.loc[df['pval']>=1]=np.nan
    df=df.dropna()
    df.to_csv(os.path.join(constants.RAW_DIR,"gwas", "filtered", gwas_file_name), sep='\t', index=False)
    print "done filter corrupted tuples. About to run PASCAL"
    print subprocess.Popen("./Pascal --pval='{}'".format(os.path.join(constants.RAW_DIR,"gwas","filtered", gwas_file_name)),
                           cwd=os.path.join(pascal_base_folder), shell=True, stdout=subprocess.PIPE).stdout.read()
    print "done run PASCAL. About format results"
    df_results=pd.read_csv(os.path.join(pascal_base_folder,"output",output_file_name), sep='\t')[['gene_symbol','pvalue']]
    df_results['eid']=pd.Series([g2e_convertor([cur])[0] if len(g2e_convertor([cur]))>0 else np.nan for cur in df_results['gene_symbol'].values], index=df_results.index)
    df_results=df_results.dropna()
    df_results.index=df_results['eid']
    df_results=df_results[['pvalue']].rename(columns={"pvalue": "pval"})
    print "done format results. About to save results"
    df_results.to_csv(os.path.join(constants.RAW_DIR,"output", gwas_file_name_format.format("_pval")),sep='\t', index_label='id')

