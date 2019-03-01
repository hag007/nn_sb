import os
import constants

import matplotlib
import os

import matplotlib

import constants

matplotlib.use('Agg')
import matplotlib.pyplot as plt

from sklearn.manifold import TSNE
import pandas as pd

df=pd.read_csv(os.path.join("/home/hag007", "CDKexp.A2058_tpm.txt"), sep='\t', index_col=0)

df.values

n_components=2
X = TSNE(n_components=n_components, metric="correlation", perplexity=30.0).fit_transform(df.T.values)
fig = plt.figure(1, figsize=(20, 20))
ax = fig.add_subplot(111)
ax.scatter(X[:, 0], X[:, 1], cmap='jet')
plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "CDKexp.A2058_tpm.png"))


