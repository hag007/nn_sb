import pandas as pd
import constants
import os
import numpy as np
constants.update_dirs(DATASET_NAME_u="MCF7_2")


ge_raw = pd.read_csv(os.path.join(constants.DATA_DIR, "ge_raw.tsv"), sep="\t", index_col=0)
enst2ensg = pd.read_csv(os.path.join(constants.DICTIONARIES_DIR, "enst2ensg.txt"), sep="\t", index_col=0)
result = pd.concat([enst2ensg, ge_raw], axis=1, join='inner')
result = result.drop_duplicates('id')
result = result.set_index("id")
result.to_csv(os.path.join(constants.DATA_DIR, "ge.tsv"), sep="\t")
for col, col_type in result.items():
    result[col] = result[col].astype(np.float)
result = result.apply(lambda x: np.log2(1 + x))
result.to_csv(os.path.join(constants.DATA_DIR, "ge_log2.tsv"), sep="\t")
