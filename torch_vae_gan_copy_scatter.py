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
import cancer_type_dataset
from torch.nn import functional as F
import constants

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as ml_colors

from matplotlib.lines import Line2D

from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

batch_size_train=100
batch_size_val=10
num_workers=25

datasets=cancer_type_dataset.CANCER_TYPES
torch_dataset=cancer_type_dataset.CancerTypesDataset(dataset_names=cancer_type_dataset.CANCER_TYPES, meta_groups_files=cancer_type_dataset.META_GROUPS, metagroups_names=["{}".format(x,i_x) for i_x, x in enumerate(cancer_type_dataset.CANCER_TYPES)])
# train_dataset,test_dataset = torch.utils.data.random_split(torch_dataset, [torch_dataset.__len__()-torch_dataset.__len__()/100, torch_dataset.__len__()/100])


class Encoder(nn.Module):
    def __init__(self, factor=0.5, n_mito_input_layer=cancer_type_dataset.n_input_layer, n_cancer_types=2,
                 n_latent_vector=100, n_reduction_layers=2):
        super(Encoder, self).__init__()
        self.n_reduction_layers = n_reduction_layers
        self.n_latent_vector = n_latent_vector

        for cur in np.arange(1, n_reduction_layers + 1):
            setattr(self, "fc_enc" + str(cur),
                    nn.Linear(int(n_mito_input_layer * factor ** (cur - 1)), int(n_mito_input_layer * factor ** cur)))
            setattr(self, "fc_bn_enc" + str(cur), nn.BatchNorm1d(int(n_mito_input_layer * factor ** cur)))

        self.fc_enc_l_mu = nn.Linear(int(n_mito_input_layer * factor ** n_reduction_layers), n_latent_vector)
        self.fc_bn_enc_l_mu = nn.BatchNorm1d(n_latent_vector)
        self.fc_enc_l_var = nn.Linear(int(n_mito_input_layer * factor ** n_reduction_layers), n_latent_vector)
        self.fc_bn_enc_l_var = nn.BatchNorm1d(n_latent_vector)

    def encode(self, x):
        h = x
        for cur in np.arange(1, self.n_reduction_layers + 1):
            h = getattr(self, "fc_bn_enc" + str(cur))(F.relu(getattr(self, "fc_enc" + str(cur))(h)))

        l_mu = getattr(self, "fc_bn_enc_l_mu")(F.relu(getattr(self, "fc_enc_l_mu")(h)))
        l_var = getattr(self, "fc_bn_enc_l_var")(F.relu(getattr(self, "fc_enc_l_var")(h)))

        return l_mu, l_var

    def reparameterize(self, mu, logvar):
        std = torch.exp(0.5 * logvar)
        eps = torch.randn_like(std)

        return eps.mul(std).add_(mu)

    def forward(self, x):
        mu, logvar = self.encode(x)
        z = self.reparameterize(mu, logvar)
        return z


class Decoder(nn.Module):

    def __init__(self, factor=0.5, n_mito_input_layer=cancer_type_dataset.n_input_layer, n_cancer_types=2,
                 n_latent_vector=100, n_reduction_layers=2):
        super(Decoder, self).__init__()
        self.n_reduction_layers = n_reduction_layers
        self.n_latent_vector = n_latent_vector

        self.fc_dec_l = nn.Linear(n_latent_vector, int(n_mito_input_layer * factor ** n_reduction_layers))
        self.fc_bn_dec_l = nn.BatchNorm1d(int(n_mito_input_layer * factor ** n_reduction_layers))

        for cur in np.arange(n_reduction_layers, 1, -1):
            setattr(self, "fc_dec" + str(cur),
                    nn.Linear(int(n_mito_input_layer * factor ** cur), int(n_mito_input_layer * factor ** (cur - 1))))
            setattr(self, "fc_bn_dec" + str(cur), nn.BatchNorm1d(int(n_mito_input_layer * factor ** (cur - 1))))
        setattr(self, "fc_dec1",
                nn.Linear(int(n_mito_input_layer * factor), int(n_mito_input_layer)))

    def decode(self, z):

        h = getattr(self, "fc_bn_dec_l")(F.relu(getattr(self, "fc_dec_l")(z)))
        for cur in np.arange(self.n_reduction_layers, 1, -1):
            h = getattr(self, "fc_bn_dec" + str(cur))(F.relu(getattr(self, "fc_dec" + str(cur))(h)))

        h = F.sigmoid(getattr(self, "fc_dec1")(h))

        return h

    def forward(self, input):
        z = input

        decoded = self.decode(z)
        return decoded


class VAE_GAN_Generator(nn.Module):
    def __init__(self, factor=0.5, n_mito_input_layer=cancer_type_dataset.n_input_layer, n_cancer_types=2,
                 n_latent_vector=100, n_reduction_layers=2):
        super(VAE_GAN_Generator, self).__init__()
        self.n_reduction_layers = n_reduction_layers
        self.n_latent_vector = n_latent_vector

        self.encoder = Encoder(factor=0.5, n_mito_input_layer=cancer_type_dataset.n_input_layer, n_cancer_types=2,
                               n_latent_vector=100, n_reduction_layers=2)
        self.decoder = Decoder(factor=0.5, n_mito_input_layer=cancer_type_dataset.n_input_layer, n_cancer_types=2,
                               n_latent_vector=100, n_reduction_layers=2)

    def forward(self, x):
        z, mu, logvar = self.encoder(x)
        rec_images = self.decoder(z)

        return mu, logvar, rec_images


class Discriminator(nn.Module):

    def __init__(self, factor=0.5, n_mito_input_layer=cancer_type_dataset.n_input_layer, n_cancer_types=2,
                 n_latent_vector=100, n_reduction_layers=2):
        super(Discriminator, self).__init__()
        # self.factor = factor
        # self.n_mito_input_layer=n_mito_input_layer
        self.n_reduction_layers = n_reduction_layers
        self.n_latent_vector = n_latent_vector

        for cur in np.arange(1, n_reduction_layers + 1):
            setattr(self, "fc_dis" + str(cur),
                    nn.Linear(int(n_mito_input_layer * factor ** (cur - 1)), int(n_mito_input_layer * factor ** cur)))
            setattr(self, "fc_bn_dis" + str(cur), nn.BatchNorm1d(int(n_mito_input_layer * factor ** cur)))

        self.fc_dis_l = nn.Linear(int(n_mito_input_layer * factor ** n_reduction_layers), n_latent_vector)
        self.fc_bn_dis_l = nn.BatchNorm1d(n_latent_vector)
        self.fc_out = nn.Linear(int(n_latent_vector), 1)

    def discriminate(self, x_hat):

        h = x_hat
        for cur in np.arange(1, self.n_reduction_layers + 1):
            h = getattr(self, "fc_bn_dis" + str(cur))(F.relu(getattr(self, "fc_dis" + str(cur))(h)))

        l = F.sigmoid(getattr(self, "fc_dis_l")(h))

        out_dis = F.sigmoid(self.fc_out(getattr(self, "fc_bn_dis_l")(l)))

        return out_dis, l

    def forward(self, input):
        encoded = input
        dis_prediction, l = self.discriminate(encoded)
        return dis_prediction, l

    def similarity(self, x):
        batch_size = x.size()[0]
        outputs, features = self.discriminate(x)
        return features


ctype='average' #'average'

datasets=cancer_type_dataset.CANCER_TYPES
trainset = cancer_type_dataset.CancerTypesDataset(dataset_names=cancer_type_dataset.CANCER_TYPES, meta_groups_files=cancer_type_dataset.META_GROUPS, metagroups_names=cancer_type_dataset.CANCER_TYPES)
trainloader = torch.utils.data.DataLoader(trainset, batch_size=10,
                                          shuffle=True, num_workers=5, pin_memory=True)
testset = trainset
testloader = trainloader
n_latent_vector=2
encoder=Encoder(n_latent_vector=n_latent_vector)
decoder=Decoder(n_latent_vector=n_latent_vector)
discriminator=Discriminator(n_latent_vector=n_latent_vector)



correct = 0
total = 0
X = None
X_z = None
X_mu = None
X_var = None
y = []

model_base_folder="/home/hag007/Desktop/nn/"
PATH_DISCRIMINATOR= model_base_folder+"GAN_DIS_mdl"# os.path.join(constants.OUTPUT_GLOBAL_DIR, "VAE_model")
PATH_ENCODER= model_base_folder+"GAN_ENC_mdl"
PATH_DECODER= model_base_folder+"GAN_DEC_mdl"
load_model=True # False
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
    for i, data in enumerate(trainloader, 0):
        features, labels = data
        _, labels = torch.max(labels, 1)
        result = m_FULL(features)
        mu, var = encoder.encode(features)
        z = encoder(features)
        decoded = decoder(z)

        if len(result) == 2:
            auth, l = m_FULL(features)
        else:
            auth, l, decoded, z, mu, var = m_FULL(features)
        X_z = np.append(X_z, z, axis=0) if X_z is not None else z
        X_mu=np.append(X_mu, mu, axis=0) if X_mu is not None else mu
        X_var=np.append(X_var, var, axis=0) if X_var is not None else var
        y = np.append(y, labels)

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
           zip(torch_dataset.get_labels_unique(), colorlist_unique)]
ax.legend(handles=patches)

plt.savefig(
    os.path.join(constants.BASE_PROFILE, "output", "AE_by_samples_z.png"))

n_components = 2
fig = plt.figure(1, figsize=(20, 20))
plt.clf()
if n_components == 3:
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(X_mu[:, 0], X_mu[:, 1], X_mu[:, 2], c=y, cmap='jet')
if n_components == 2:
    ax = fig.add_subplot(111)
    ax.scatter(X_mu[:, 0], X_mu[:, 1], c=y, cmap='jet')

colorlist_unique = [ml_colors.rgb2hex(colormap(a)) for a in
                    label_ids_unique / float(max(label_ids))]
patches = [Line2D([0], [0], marker='o', color='gray', label=a,
                  markerfacecolor=c) for a, c in
           zip(torch_dataset.get_labels_unique(), colorlist_unique)]
ax.legend(handles=patches)

plt.savefig(
    os.path.join(constants.BASE_PROFILE, "output", "AE_by_samples_mu.png"))

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
           zip(torch_dataset.get_labels_unique(), colorlist_unique)]
ax.legend(handles=patches)

plt.savefig(
    os.path.join(constants.BASE_PROFILE, "output", "AE_by_samples_logvar.png"))

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
