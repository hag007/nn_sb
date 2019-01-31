import pandas as pd
import constants
import os

constants.update_dirs(DATASET_NAME_u="IES")
df_ge = pd.read_csv(os.path.join(constants.DATA_DIR, "ge_mouse.tsv"), sep="\t")
df_ge = df_ge.set_index("id")
df_ge = df_ge[~df_ge.index.duplicated(keep='first')]
df_mouse2human = pd.read_csv(os.path.join(constants.DICTIONARIES_DIR, "mouse2human.txt"), sep="\t")
df_mouse2human = df_mouse2human.set_index("Mouse gene stable ID")
df_mouse2human = df_mouse2human[~df_mouse2human.index.duplicated(keep='first')]
df_converted_ge = pd.concat([df_ge, df_mouse2human], join="inner", axis=1)
df_converted_ge = df_converted_ge.set_index("Gene stable ID")
df_converted_ge = df_converted_ge[~df_converted_ge.index.duplicated(keep="first")]
df_converted_ge.to_csv(os.path.join(constants.DATA_DIR, "ge.tsv"), sep="\t", index_label="id")

