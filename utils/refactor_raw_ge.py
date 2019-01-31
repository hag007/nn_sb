import pandas as pd
import numpy as np
import constants
import os
import ensembl2gene_symbol


# # convert ensembl to gene symbols
# constants.update_dirs(DATASET_NAME_u="GE_ERS_2")
# df_raw = pd.read_csv(os.path.join(constants.DATA_DIR, "ge_raw.tsv"), sep="\t", index_col=0)
# ensg_ids = [ensembl2gene_symbol.g2e_convertor([cur])[0] if len(ensembl2gene_symbol.g2e_convertor([cur])) > 0 else "nan" for cur in df_raw.index]
# ensg_ids=[x if not x.startswith("ENSGR") else "nan" for x in ensg_ids ]
# df_raw.index=ensg_ids
# df_raw=df_raw[df_raw.index!="nan"]
# df_raw=df_raw[(df_raw.T != 0).any()]
# # df_raw=df_raw.dropna()
# df_raw= df_raw[~df_raw.index.duplicated(keep="first")]
# df_raw.to_csv(os.path.join(constants.DATA_DIR, "ge.tsv"), sep="\t", index_label="id")



# convert ensembl to gene symbols
constants.update_dirs(DATASET_NAME_u="EN_BRCA")
df_raw = pd.read_csv(os.path.join(constants.DATA_DIR, "ge_raw.tsv"), sep="\t", index_col=0)
ensg_ids = [ensembl2gene_symbol.g2e_convertor([cur])[0] if len(ensembl2gene_symbol.g2e_convertor([cur])) > 0 else "nan" for cur in df_raw.index]
ensg_ids=[x if not x.startswith("ENSGR") else "nan" for x in ensg_ids ]
df_raw.index=ensg_ids
df_raw=df_raw[df_raw.index!="nan"]
df_raw=df_raw[(df_raw.T != 0).any()]
# df_raw=df_raw.dropna()
df_raw= df_raw[~df_raw.index.duplicated(keep="first")]
df_raw.to_csv(os.path.join(constants.DATA_DIR, "score.tsv"), sep="\t", index_label="id")



# # remove ensembl "tail" dot
# constants.update_dirs(DATASET_NAME_u="GE_SHEZH_2")
# df_raw = pd.read_csv(os.path.join(constants.DATA_DIR, "ge_raw.tsv"), sep="\t", index_col=0)
# ensg_ids = [cur.split(".")[0] for cur in df_raw.index]
# ensg_ids=[x if not x.startswith("ENSGR") else "nan" for x in ensg_ids ]
# df_raw.index=ensg_ids
# df_raw=df_raw[(df_raw.T != 0).any()]
# df_raw=df_raw[df_raw.index!="nan"]
# df_raw=df_raw.dropna()
# df_raw= df_raw[~df_raw.index.duplicated(keep="first")]
# df_raw.to_csv(os.path.join(constants.DATA_DIR, "ge.tsv"), sep="\t", index_label="id")

