import os

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as ml_colors

import numpy as np
import simplejson as json
import torch
import torch.nn as nn

import torch_dataset_cancer
import infra
from torch_dataset_cancer import CancerTypesDataset

matplotlib.use('Agg')
from utils import param_builder
import pandas as pd
import torch_vae_gan_copy_model
import scipy.stats
import constants

datasets=torch_dataset_cancer.CANCER_TYPES
trainset = CancerTypesDataset(dataset_names=torch_dataset_cancer.CANCER_TYPES, meta_groups_files=torch_dataset_cancer.META_GROUPS, metagroups_names=torch_dataset_cancer.CANCER_TYPES)
trainloader = torch.utils.data.DataLoader(trainset, batch_size=10,
                                          shuffle=True, num_workers=5, pin_memory=True)
testset = trainset
testloader = trainloader
n_latent_vector=2
encoder=torch_vae_gan_copy_model.Encoder(n_latent_vector=n_latent_vector)
decoder=torch_vae_gan_copy_model.Decoder(n_latent_vector=n_latent_vector)
discriminator=torch_vae_gan_copy_model.Discriminator(n_latent_vector=n_latent_vector)

model_base_folder="/home/hag007/Desktop/nn/"
PATH_DISCRIMINATOR= model_base_folder+"GAN_DIS_mdl"
PATH_ENCODER= model_base_folder+"GAN_ENC_mdl"
PATH_DECODER= model_base_folder+"GAN_DEC_mdl"
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

correct = 0
total = 0
X = None
X_z = None
X_mu = None
X_var = None
y_survival = []
y_labels = []
epoch="inf"
m_FULL.train(False)

with torch.no_grad():
    for i in range(trainset.__len__()):
        try:
            features, labels, survival = trainset.get_full_item(i)
        except:
            continue
        features=torch.stack([features])
        labels=torch.stack([labels])
        # survival = torch.stack([survival])
        _, labels = torch.max(labels, 1)
        result=m_FULL(features)
        mu, var=encoder.encode(features)
        z = encoder(features)
        decoded = decoder(z)

        if len(result) == 2:
            auth, l = m_FULL(features)
        else:
            auth, l, decoded, z, mu, var = m_FULL(features)

        X_z = np.append(X_z, z, axis=0) if X_z is not None else z
        # mu=features.numpy()
        X_mu= np.append(X_mu, mu, axis=0) if X_mu is not None else mu
        X_var=np.append(X_var, var, axis=0) if X_var is not None else var
        y_survival = np.append(y_survival, int(survival[4]))
        y_labels = np.append(y_survival, labels)



y=y_survival
n_components = 2
vmin=np.percentile(y_survival,2)
vmax=np.percentile(y_survival,90)
fig = plt.figure(1, figsize=(20, 20))
plt.clf()
if n_components == 3:
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(X_z[:, 0], X_z[:, 1], X_z[:, 2], c=y, cmap='jet', vmin=vmin, vmax=vmax)
if n_components == 2:
    ax = fig.add_subplot(111)
    ax.scatter(X_z[:, 0], X_z[:, 1], c=y, cmap='jet', vmin=vmin, vmax=vmax)

ax.legend()

plt.savefig(
    os.path.join(constants.BASE_PROFILE, "output", "AE_by_samples_z.png"))

n_components = 2
fig = plt.figure(1, figsize=(20, 20))
plt.clf()
if n_components == 3:
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(X_mu[:, 0], X_mu[:, 1], X_mu[:, 2], c=y, cmap='jet', vmin=vmin, vmax=vmax)
if n_components == 2:
    ax = fig.add_subplot(111)
    ax.scatter(X_mu[:, 0], X_mu[:, 1], c=y, cmap='jet', vmin=vmin, vmax=vmax)

ax.legend()

plt.savefig(
    os.path.join(constants.BASE_PROFILE, "output", "AE_by_samples_mu.png"))

n_components = 2
fig = plt.figure(1, figsize=(20, 20))
plt.clf()
if n_components == 3:
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(X_var[:, 0], X_var[:, 1], X_var[:, 2], c=y, cmap='jet', vmin=vmin, vmax=vmax)
if n_components == 2:
    ax = fig.add_subplot(111)
    ax.scatter(X_var[:, 0], X_var[:, 1], c=y, cmap='jet', vmin=vmin, vmax=vmax)

ax.legend()

plt.savefig(
    os.path.join(constants.BASE_PROFILE, "output", "AE_by_samples_logvar.png"))







# X_mu = PCA(n_components=100).fit_transform(X_mu)

# label_ids_unique = np.unique(y_labels)
# label_ids = [trainset.get_labels_unique()[int(a)] for a in y_labels]
#
#
# samples_z_dict={}
# samples_mu_dict={}
# samples_var_dict={}
# for label, z, mu, var in zip(label_ids, X_z, X_mu, X_var):
#     if "unknown" not in label:
#         samples_z_dict[str(label)]=samples_z_dict[str(label)] + [z] if str(label) in samples_z_dict else []
#         samples_mu_dict[str(label)] = samples_mu_dict[str(label)] + [mu] if str(label) in samples_mu_dict else []
#         samples_var_dict[str(label)] = samples_var_dict[str(label)] + [var] if str(label) in samples_var_dict else []

# samples_z_mean_dict={}
# samples_mu_mean_dict={}
# samples_var_mean_dict={}
# for k,v in samples_z_dict.iteritems():
#     samples_z_mean_dict[k]=np.mean(v, axis=0)
#
# for k,v in samples_mu_dict.iteritems():
#     samples_mu_mean_dict[k]=np.mean(v, axis=0)
#
# for k,v in samples_var_dict.iteritems():
#     samples_var_mean_dict[k]=np.mean(v, axis=0)


# dataset = cancer_type_dataset.CANCER_TYPES[0]
# data_normalization = "fpkm"
# tested_gene_list_file_name = "protein_coding_long.txt"
#
# gene_expression_file_name, phenotype_file_name, survival_file_name, mutation_file_name, mirna_file_name, pval_preprocessing_file_name = \
#     param_builder.build_gdc_params(dataset, data_normalization)
#
#
# meta_groups = [json.load(file("groups/temp.json"))]
# filter_expression = None
# tested_gene_expression, tested_gene_expression_headers_rows, tested_gene_expression_headers_columns, labels_assignment, survival_dataset = \
#     infra.load_integrated_ge_data(tested_gene_list_file_name=tested_gene_list_file_name,
#                                   total_gene_list_file_name=tested_gene_list_file_name,
#                                   gene_expression_file_name=gene_expression_file_name,
#                                   phenotype_file_name=phenotype_file_name,
#                                   survival_file_name=survival_file_name,
#                                   var_th_index=None, meta_groups=meta_groups,
#                                   filter_expression=filter_expression)
#
# df_survival = pd.DataFrame(np.array(survival_dataset[1:, 1:]), index=survival_dataset[1:, 0],
#                            columns=survival_dataset[0, 1:])
#
# normal_vector_mean= samples_mu_mean_dict["{}, normal".format(dataset)]

# overall_survivals=[]
# distances=[]
# with torch.no_grad():
#     for i in range(trainset.__len__()):
#         try:
#             features, label, survival = trainset.get_full_item(i)
#         except:
#             continue
#         mu, var = encoder.encode(torch.stack([features]))
#
#         try:
#             distance=np.linalg.norm(mu.numpy() - normal_vector_mean)
#             overall_survival=df_survival.loc[tested_gene_expression_headers_rows[i]].iloc[4]
#             distances.append(distance)
#             overall_survivals.append(overall_survival)
#         except:
#             pass
#
#
#
#     print scipy.stats.pearsonr(np.array(distances, dtype=np.float), np.array(overall_survivals,dtype=np.float))
#     print scipy.stats.spearmanr(np.array(distances, dtype=np.float), np.array(overall_survivals, dtype=np.float))
