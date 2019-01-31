import numpy as np
import pandas as pd

def to_full_np(df, index_title):
    if type(df) is list:
        if len(df) ==0 : return []
        else: df = pd.DataFrame(df)
        if index_title in df:
            df = df.set_index(index_title)

    return np.r_[[[index_title] + list(df.columns)], np.c_[np.array(df.index), df.values]]

def to_full_list(df, index_title):

    return [list (x) for x in to_full_np(df, index_title)]