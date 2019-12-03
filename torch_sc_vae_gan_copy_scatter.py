import argparse
import numpy as np
import os
import time
import torch
import torch.utils.data
import torch.nn as nn
import torch.optim as optim
from torch.utils.data.dataset import Dataset
from torch.autograd import Variable
import torch_dataset_cancer
from torch.nn import functional as F
import constants

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as ml_colors
import torch_vae_gan_copy_model

import torch_dataset_sc

from matplotlib.lines import Line2D

from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

batch_size_train=100
batch_size_val=10
num_workers=25

ctype='average' #'average'

trainset = torch_dataset_sc.SingleCellDataset(dataset_names=torch_dataset_sc.DATASETS, meta_groups_files=torch_dataset_cancer.META_GROUPS, metagroups_names=torch_dataset_sc.DATASETS)
trainloader = torch.utils.data.DataLoader(trainset, batch_size=10,
                                          shuffle=True, num_workers=5, pin_memory=True)


# trainset = torch_dataset_cancer.CancerTypesDataset(dataset_names=torch_dataset_cancer.CANCER_TYPES, meta_groups_files=torch_dataset_cancer.META_GROUPS, metagroups_names=torch_dataset_cancer.CANCER_TYPES)
# trainloader = torch.utils.data.DataLoader(trainset, batch_size=10,
#                                           shuffle=True, num_workers=5, pin_memory=True)

testset = trainset
testloader = trainloader
n_latent_vector=2
encoder=torch_vae_gan_copy_model.Encoder(n_latent_vector=n_latent_vector)
decoder=torch_vae_gan_copy_model.Decoder(n_latent_vector=n_latent_vector)
discriminator=torch_vae_gan_copy_model.Discriminator(n_latent_vector=n_latent_vector)


correct = 0
total = 0
X = None
X_z = None
X_mu = None
X_var = None
X_survival = None
y = []

model_base_folder="/home/hag007/Desktop/nn/"
PATH_DISCRIMINATOR= model_base_folder+"VAE_DIS_mdl"# os.path.join(constants.OUTPUT_GLOBAL_DIR, "VAE_model")
PATH_ENCODER= model_base_folder+"VAE_ENC_mdl"
PATH_DECODER= model_base_folder+"VAE_DEC_mdl"
load_model=True
if load_model and os.path.exists(PATH_ENCODER):
    encoder.load_state_dict(torch.load(PATH_ENCODER))
    encoder.eval()
    decoder.load_state_dict(torch.load(PATH_DECODER))
    decoder.eval()
    discriminator.load_state_dict(torch.load(PATH_DISCRIMINATOR))
    discriminator.eval()

m_VAE = nn.Sequential(encoder,decoder)
m_GAN = nn.Sequential(decoder,discriminator)
m_FULL = nn.Sequential(encoder,decoder,discriminator)


with torch.no_grad():
    for i in range(trainset.__len__()):
        try:
            features, label = trainset.__getitem__(i) # get_full_item(i)
            # if survival.iloc[3] == '0': continue
        except:
            continue
        features=torch.stack([features, features])
        label=torch.stack([label])
        mu, var = encoder.encode(features)[0]
        mu=[mu.numpy()]
        var = [var.numpy()]
        z = [encoder(features)[0].numpy()]

        X_z = np.append(X_z, z, axis=0) if X_z is not None else z
        X_mu=np.append(X_mu, mu, axis=0) if X_mu is not None else mu
        X_var=np.append(X_var, var, axis=0) if X_var is not None else var
        # X_survival = np.append(X_survival, [survival.iloc[4]], axis=0) if X_survival is not None else [survival.iloc[4]]
        y = np.append(y, torch.argmax(label).numpy())

# X_survival=X_survival.astype(np.int)

print "len samples: {}".format(len(X_mu))
colormap = cm.jet
label_ids_unique = np.unique(y)
label_ids = y

n_components = 2
fig = plt.figure(1, figsize=(20, 20))
plt.clf()
if n_components == 3:
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(X_z[:, 0], X_z[:, 1], X_z[:, 2], c=y, cmap='jet')
if n_components == 2:
    ax = fig.add_subplot(111)
    ax.scatter(X_z[:, 0], X_z[:, 1], c=y, cmap='jet')

colorlist_unique = [ml_colors.rgb2hex(colormap(a)) for a in
                    label_ids_unique / float(max(label_ids))]
patches = [Line2D([0], [0], marker='o', color='gray', label=a,
                  markerfacecolor=c) for a, c in
           zip(trainset.get_labels_unique(), colorlist_unique)]
ax.legend(handles=patches)

plt.savefig(
    os.path.join(constants.BASE_PROFILE, "output", "AE_by_samples_sc_z.png"))

n_components = 2
fig = plt.figure(1, figsize=(20, 20))
plt.clf()
if n_components == 3:
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(X_mu[:, 0], X_mu[:, 1], X_mu[:, 2], c=y, cmap='jet')
if n_components == 2:
    ax = fig.add_subplot(111)
    # ax.scatter(X_mu[:, 0], X_mu[:, 1], c=X_survival, cmap='jet', vmin=np.percentile(X_survival,20), vmax=np.percentile(X_survival,80))
    ax.scatter(X_mu[:, 0], X_mu[:, 1], c=y, cmap='jet')

reduced_dict={}
for x1, x2, y_i in zip(X_mu[:, 0], X_mu[:, 1],y):
    if y_i not in reduced_dict:
        reduced_dict[y_i]=[]

    reduced_dict[y_i].append([x1, x2])

for k, v in reduced_dict.iteritems():
    print trainset.get_labels_unique()[int(k)], np.var(np.array(reduced_dict[k]), axis=0), np.linalg.norm(np.var(np.array(reduced_dict[k]), axis=0)), np.array(reduced_dict[k]).shape[0]

colorlist_unique = [ml_colors.rgb2hex(colormap(a)) for a in
                    label_ids_unique / float(max(label_ids))]
patches = [Line2D([0], [0], marker='o', color='gray', label=a,
                  markerfacecolor=c) for a, c in
           zip(trainset.get_labels_unique(), colorlist_unique)]
ax.legend(handles=patches)

plt.savefig(
    os.path.join(constants.BASE_PROFILE, "output", "AE_by_samples_sc_mu.png"))

n_components = 2
fig = plt.figure(1, figsize=(20, 20))
plt.clf()
if n_components == 3:
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(X_var[:, 0], X_var[:, 1], X_var[:, 2], c=y, cmap='jet')
if n_components == 2:
    ax = fig.add_subplot(111)
    ax.scatter(X_var[:, 0], X_var[:, 1], c=y, cmap='jet')

colorlist_unique = [ml_colors.rgb2hex(colormap(a)) for a in
                    label_ids_unique / float(max(label_ids))]
patches = [Line2D([0], [0], marker='o', color='gray', label=a,
                  markerfacecolor=c) for a, c in
           zip(trainset.get_labels_unique(), colorlist_unique)]
ax.legend(handles=patches)

plt.savefig(
    os.path.join(constants.BASE_PROFILE, "output", "AE_by_samples_sc_logvar.png"))

# ###########################
#
# correct = 0
# total = 0
# X = None
# X_r = None
# y = []
# with torch.no_grad():
#     for i, data in enumerate(trainloader, 0):
#         features, labels = data
#         X_r = np.append(X_r, features, axis=0) if X_r is not None else features
#         _,lbl=torch.max(labels, 1)
#         y = np.append(y, lbl)
#
# n_components=2
# X = PCA(n_components=n_components).fit_transform(X_r)
# fig = plt.figure(1, figsize=(20, 20))
# plt.clf()
# if n_components == 3:
#     ax = fig.add_subplot(111, projection='3d')
#     ax.scatter(X[:, 0], X[:, 1], X[:, 2], c=y, cmap='jet')
# if n_components == 2:
#     ax = fig.add_subplot(111)
#     ax.scatter(X[:, 0], X[:, 1], c=y, cmap='jet')
#
# colormap = cm.jet
# label_ids_unique = np.unique(y)
# label_ids = y
#
# colorlist_unique = [ml_colors.rgb2hex(colormap(a)) for a in
#                     label_ids_unique / float(max(label_ids))]
# patches = [Line2D([0], [0], marker='o', color='gray', label=a,
#                   markerfacecolor=c) for a, c in
#            zip(datasets, colorlist_unique)]
# ax.legend(handles=patches)
#
# plt.savefig(
#     os.path.join(constants.BASE_PROFILE, "output", "PCA_by_samples.png").format(constants.CANCER_TYPE))
#
#
#
# n_components=2
# X = TSNE(n_components=n_components, metric="correlation", perplexity=30.0).fit_transform(X_r)
# fig = plt.figure(1, figsize=(20, 20))
# plt.clf()
# if n_components == 3:
#     ax = fig.add_subplot(111, projection='3d')
#     ax.scatter(X[:, 0], X[:, 1], X[:, 2], c=y, cmap='jet')
# if n_components == 2:
#     ax = fig.add_subplot(111)
#     ax.scatter(X[:, 0], X[:, 1], c=y, cmap='jet')
#
# colormap = cm.jet
# label_ids_unique = np.unique(y)
# label_ids = y
#
# colorlist_unique = [ml_colors.rgb2hex(colormap(a)) for a in
#                     label_ids_unique / float(max(label_ids))]
# patches = [Line2D([0], [0], marker='o', color='gray', label=a,
#                   markerfacecolor=c) for a, c in
#            zip(datasets, colorlist_unique)]
# ax.legend(handles=patches)
#
# plt.savefig(
#     os.path.join(constants.BASE_PROFILE, "output", "TSNE_by_samples.png").format(constants.CANCER_TYPE))

