import os
import constants

import matplotlib
import os

import matplotlib

import constants

matplotlib.use('Agg')
import matplotlib.pyplot as plt

from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
import pandas as pd
import jarvispatrick
import numpy as np
base_folder="/specific/netapp5/gaga/hagailevi/sc/datasets"
data_folder='test' # 'sc' # 'sc_zebrafish'
file_name='GSE100866_CBMC_8K_13AB_10X-RNA_umi.csv' # 'CDKexp.A2058_tpm.txt' # 'URD_Dropseq_Expression_Log2TPM.txt'

sep=',' # '\t'

df=pd.read_csv(os.path.join(base_folder, data_folder, file_name), sep=sep, index_col=0)



n_components=2
print "start pca..."
X_pca = PCA(n_components=10).fit_transform(df.T.values)

clustering = jarvispatrick.JarvisPatrick([i for i in range(X_pca.shape[0])], lambda x,y: np.linalg.norm(X_pca[x]-X_pca[y]))
sol= clustering(1,0)
print sol
print "start tsne..."
X = TSNE(n_components=n_components, metric="euclidean", perplexity=15.0).fit_transform(X_pca)
fig = plt.figure(1, figsize=(20, 20))
ax = fig.add_subplot(111)
ax.scatter(X[:, 0], X[:, 1], cmap='jet')
plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, file_name.split('.')[0]+".png"))



