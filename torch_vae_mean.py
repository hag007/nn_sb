import torch
import torch.nn as nn
from torch.nn import functional as F
import torch.optim as optim
import os
import numpy as np
import constants
from cancer_type_dataset import CancerTypesDataset
import cancer_type_dataset
import simplejson as json
from utils.param_builder import build_gdc_params

from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as ml_colors

from matplotlib.lines import Line2D

import vae_model
import vae_bn_after_relu_model


from scipy.cluster.hierarchy import dendrogram, linkage


    # Reconstruction + KL divergence losses summed over all elements and batch


def loss_function(recon_x, x, mu, logvar):
    BCE = F.binary_cross_entropy(recon_x, x, reduction='sum')

    # see Appendix B from VAE paper:
    # Kingma and Welling. Auto-Encoding Variational Bayes. ICLR, 2014
    # https://arxiv.org/abs/1312.6114
    # 0.5 * sum(1 + log(sigma^2) - mu^2 - sigma^2)
    KLD = -0.5 * torch.sum(1 + logvar - mu.pow(2) - logvar.exp())

    return BCE + KLD

datasets=cancer_type_dataset.CANCER_TYPES
trainset = CancerTypesDataset(dataset_names=cancer_type_dataset.CANCER_TYPES, meta_groups_files=cancer_type_dataset.META_GROUPS, metagroups_names=cancer_type_dataset.CANCER_TYPES)
trainloader = torch.utils.data.DataLoader(trainset, batch_size=10,
                                          shuffle=True, num_workers=5, pin_memory=True)
testset = trainset
testloader = trainloader

criterion = nn.BCELoss()

net = vae_bn_after_relu_model.Net()

# create your optimizer
optimizer = optim.Adam(net.parameters(), lr=0.00001)


correct = 0
total = 0
X = None
X_z = None
X_mu = None
X_var = None
y = []



PATH="/home/hag007/Desktop/nn/VAE_model"
net.load_state_dict(torch.load(PATH))
net.eval()
epoch="inf"
with torch.no_grad():
    for i in range(trainset.__len__()):
        features, labels = trainset.__getitem__(i)
        features=torch.stack([features])
        labels=torch.stack([labels])
        _, labels = torch.max(labels, 1)
        outputs, z, mu, var = net(features)
        X_z = np.append(X_z, z, axis=0) if X_z is not None else z
        # mu=features.numpy()
        X_mu= np.append(X_mu, mu, axis=0) if X_mu is not None else mu
        X_var=np.append(X_var, var, axis=0) if X_var is not None else var
        y = np.append(y, labels)


# X_mu = PCA(n_components=2).fit_transform(X_mu)

label_ids_unique = np.unique(y)
label_ids = [trainset.get_labels_unique()[int(a)] for a in y]


samples_z_dict={}
samples_mu_dict={}
samples_var_dict={}
for label, z, mu, var in zip(label_ids, X_z, X_mu, X_var):
    if "unknown" not in label:
        samples_z_dict[str(label)]=samples_z_dict[str(label)] + [z] if str(label) in samples_z_dict else []
        samples_mu_dict[str(label)] = samples_mu_dict[str(label)] + [mu] if str(label) in samples_mu_dict else []
        samples_var_dict[str(label)] = samples_var_dict[str(label)] + [var] if str(label) in samples_var_dict else []


samples_z_mean_dict={}
samples_mu_mean_dict={}
samples_var_mean_dict={}

for k,v in samples_z_dict.iteritems():
    samples_z_mean_dict[k]=np.mean(v, axis=0)

for k,v in samples_mu_dict.iteritems():
    samples_mu_mean_dict[k]=np.mean(v, axis=0)

for k in sorted(samples_mu_dict.keys()):
    print "{}: {}".format(k,samples_mu_mean_dict[k])

for k,v in samples_var_dict.iteritems():
    samples_var_mean_dict[k]=np.mean(v, axis=0)

metagroups_names=[x.split(',')[0] for x in samples_z_mean_dict.keys()]


tumor_only_vectors_z={}
tumor_only_vectors_mu={}
tumor_only_vectors_var={}
for cur in metagroups_names:
    normal_key="{}, {}".format(cur, "normal")
    tumor_key = "{}, {}".format(cur, "tumor")
    tumor_only_vectors_z[cur]=np.subtract(samples_z_mean_dict[tumor_key], samples_z_mean_dict[normal_key])
    tumor_only_vectors_mu[cur] =np.subtract(samples_mu_mean_dict[tumor_key], samples_mu_mean_dict[normal_key])
    tumor_only_vectors_var[cur] =np.subtract(samples_var_mean_dict[tumor_key], samples_var_mean_dict[normal_key])


n_components = 2

if n_components == 2:
    fig = plt.figure(1, figsize=(10, 10))
    plt.clf()
    ax = fig.add_subplot(111)
    for k,v in tumor_only_vectors_mu.iteritems():
        ax.scatter(v[0], v[1])
        ax.annotate(k,
                    xy=(v[0], v[1]), xytext=(-20, 20), textcoords='offset points',
                    bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.5),
                    arrowprops=dict(facecolor='black', shrink=0.05, width=2, headwidth=3, headlength=2))
    plt.savefig(
        os.path.join(constants.BASE_PROFILE, "output", "mean_sub_by_samples_z.png"))

    fig = plt.figure(1, figsize=(10, 10))
    plt.clf()
    ax = fig.add_subplot(111)
    filtered_tumor_only_vectors_mu_values=[]
    filtered_tumor_only_vectors_mu_keys=[]
    for k,v in tumor_only_vectors_mu.iteritems():
        ax.scatter(v[0], v[1])
        ax.annotate(k,
                    xy=(v[0], v[1]), xytext=(-20, 20), textcoords='offset points',
                    bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.5),
                    arrowprops=dict(facecolor='black', shrink=0.05, width=2, headwidth=3, headlength=2))

        if True or (k!="KICH" and k!="PRAD") :
            filtered_tumor_only_vectors_mu_values.append(v)
            filtered_tumor_only_vectors_mu_keys.append(k)
    plt.savefig(
        os.path.join(constants.BASE_PROFILE, "output", "mean_sub_by_samples_mu.png"))

    # convert the redundant n*n square matrix form into a condensed nC2 array

    linked = linkage(filtered_tumor_only_vectors_mu_values, 'single')

    plt.figure(figsize=(10, 7))
    dendrogram(linked,
               orientation='top',
               labels=filtered_tumor_only_vectors_mu_keys,
               distance_sort='descending',
               show_leaf_counts=True)

    plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "mean_sub_by_samples_mu_hierarchical.png"))

    fig = plt.figure(1, figsize=(10, 10))
    plt.clf()
    ax = fig.add_subplot(111)
    for k,v in tumor_only_vectors_var.iteritems():
        ax.scatter(v[0], v[1])
        ax.annotate(k,
                    xy=(v[0], v[1]), xytext=(-20, 20), textcoords='offset points',
                    bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.5),
                    arrowprops=dict(facecolor='black', shrink=0.05, width=2, headwidth=3, headlength=2))
    plt.savefig(
        os.path.join(constants.BASE_PROFILE, "output", "mean_sub_by_samples_var.png"))


    fig = plt.figure(1, figsize=(10, 10))
    plt.clf()
    ax = fig.add_subplot(111)
    for k,v in samples_z_mean_dict.iteritems():
        if "normal" in k:
            ax.scatter(v[0], v[1])
            ax.annotate(k,
                    xy=(v[0], v[1]), xytext=(-20, 20), textcoords='offset points',
                    bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.5),
                    arrowprops=dict(facecolor='black', shrink=0.05, width=2, headwidth=3, headlength=2))
    plt.savefig(
        os.path.join(constants.BASE_PROFILE, "output", "mean_normal_by_samples_z.png"))

    fig = plt.figure(1, figsize=(10, 10))
    plt.clf()
    ax = fig.add_subplot(111)
    for k,v in samples_mu_mean_dict.iteritems():
        if "normal" in k:
            ax.scatter(v[0], v[1])
            ax.annotate(k,
                    xy=(v[0], v[1]), xytext=(-20, 20), textcoords='offset points',
                    bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.5),
                    arrowprops=dict(facecolor='black', shrink=0.05, width=2, headwidth=3, headlength=2))
    plt.savefig(
        os.path.join(constants.BASE_PROFILE, "output", "mean_normal_by_samples_mu.png"))

    fig = plt.figure(1, figsize=(10, 10))
    plt.clf()
    ax = fig.add_subplot(111)
    for k,v in samples_var_mean_dict.iteritems():
        if "normal" in k:
            ax.scatter(v[0], v[1])
            ax.annotate(k,
                    xy=(v[0], v[1]), xytext=(-20, 20), textcoords='offset points',
                    bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.5),
                    arrowprops=dict(facecolor='black', shrink=0.05, width=2, headwidth=3, headlength=2))
    plt.savefig(
        os.path.join(constants.BASE_PROFILE, "output", "mean_normal_by_samples_var.png"))


    fig = plt.figure(1, figsize=(10, 10))
    plt.clf()
    ax = fig.add_subplot(111)
    for k,v in samples_z_mean_dict.iteritems():
        if "tumor" in k:
            ax.scatter(v[0], v[1])
            ax.annotate(k,
                    xy=(v[0], v[1]), xytext=(-20, 20), textcoords='offset points',
                    bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.5),
                    arrowprops=dict(facecolor='black', shrink=0.05, width=2, headwidth=3, headlength=2))
    plt.savefig(
        os.path.join(constants.BASE_PROFILE, "output", "mean_tumor_by_samples_z.png"))

    fig = plt.figure(1, figsize=(10, 10))
    plt.clf()
    ax = fig.add_subplot(111)
    for k,v in samples_mu_mean_dict.iteritems():
        if "tumor" in k:
            ax.scatter(v[0], v[1])
            ax.annotate(k,
                    xy=(v[0], v[1]), xytext=(-20, 20), textcoords='offset points',
                    bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.5),
                    arrowprops=dict(facecolor='black', shrink=0.05, width=2, headwidth=3, headlength=2))

    COAD_new=(tumor_only_vectors_mu["LUAD"][0]+samples_mu_mean_dict["COAD, normal"][0], tumor_only_vectors_mu["LUAD"][1]+samples_mu_mean_dict["COAD, normal"][1])
    LUAD_new =(tumor_only_vectors_mu["COAD"][0] + samples_mu_mean_dict["LUAD, normal"][0], tumor_only_vectors_mu["COAD"][1]+samples_mu_mean_dict["LUAD, normal"][1])

    BRCA_new = (tumor_only_vectors_mu["KIRP"][0] + samples_mu_mean_dict["BRCA, normal"][0], tumor_only_vectors_mu["KIRP"][1] + samples_mu_mean_dict["BRCA, normal"][1])
    KIRP_new = (tumor_only_vectors_mu["BRCA"][0] + samples_mu_mean_dict["KIRP, normal"][0], tumor_only_vectors_mu["BRCA"][1] + samples_mu_mean_dict["KIRP, normal"][1])

    ax.scatter(*COAD_new)
    ax.scatter(*LUAD_new)

    ax.scatter(*BRCA_new)
    ax.scatter(*KIRP_new)

    ax.annotate("COAD_new",
                xy=COAD_new, xytext=(-20, 20), textcoords='offset points',
                bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.5),
                arrowprops=dict(facecolor='black', shrink=0.05, width=2, headwidth=3, headlength=2))
    ax.annotate("LUAD_new",
                xy=LUAD_new, xytext=(-20, 20), textcoords='offset points',
                bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.5),
                arrowprops=dict(facecolor='black', shrink=0.05, width=2, headwidth=3, headlength=2))

    ax.annotate("BRCA_new",
                xy=BRCA_new, xytext=(-20, 20), textcoords='offset points',
                bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.5),
                arrowprops=dict(facecolor='black', shrink=0.05, width=2, headwidth=3, headlength=2))
    ax.annotate("KIRP_new",
                xy=KIRP_new, xytext=(-20, 20), textcoords='offset points',
                bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.5),
                arrowprops=dict(facecolor='black', shrink=0.05, width=2, headwidth=3, headlength=2))


    plt.savefig(
        os.path.join(constants.BASE_PROFILE, "output", "mean_tumor_by_samples_mu.png"))

    fig = plt.figure(1, figsize=(10, 10))
    plt.clf()
    ax = fig.add_subplot(111)
    for k,v in samples_var_mean_dict.iteritems():
        if "tumor" in k:
            ax.scatter(v[0], v[1])
            ax.annotate(k,
                    xy=(v[0], v[1]), xytext=(-20, 20), textcoords='offset points',
                    bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.5),
                    arrowprops=dict(facecolor='black', shrink=0.05, width=2, headwidth=3, headlength=2))
    plt.savefig(
        os.path.join(constants.BASE_PROFILE, "output", "mean_tumor_by_samples_var.png"))
