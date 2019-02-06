import torch
import torch.nn as nn
from torch.nn import functional as F
import torch.optim as optim
import os
import numpy as np
import constants
from cancer_type_dataset import CancerTypesDataset
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


class Net(nn.Module):

    def __init__(self, factor=0.5, n_mito_input_layer=100, n_cancer_types=2, n_latent_vector=2):
        super(Net, self).__init__()
        # self.factor = factor
        # self.n_mito_input_layer=n_mito_input_layer

        self.fc1 = nn.Linear(n_mito_input_layer, int(n_mito_input_layer * factor))
        self.fc2 = nn.Linear(int(n_mito_input_layer * factor), int(n_mito_input_layer * factor ** 2))
        self.fc31 = nn.Linear(int(n_mito_input_layer * factor ** 2), n_latent_vector)
        self.fc32 = nn.Linear(int(n_mito_input_layer * factor ** 2), n_latent_vector)
        self.fc4 = nn.Linear(n_latent_vector, int(n_mito_input_layer * factor ** 2))
        self.fc5 = nn.Linear(int(n_mito_input_layer * factor ** 2), int(n_mito_input_layer * factor))
        self.fc6 = nn.Linear(int(n_mito_input_layer * factor), n_mito_input_layer)

    def encode(self, x):
        h1 = F.relu(self.fc1(x))
        h2 = F.relu(self.fc2(h1))
        h31 = self.fc31(h2)
        h32 = self.fc32(h2)
        return h31, h32

    def reparameterize(self, mu, logvar):
        std = torch.exp(0.5 * logvar)
        eps = torch.randn_like(std)
        return eps.mul(std).add_(mu)

    def decode(self, z):
        h4 = F.relu(self.fc4(z))
        h5 = F.relu(self.fc5(h4))
        h6 = F.sigmoid(self.fc6(h5))
        return h6

    def forward(self, x):
        mu, logvar = self.encode(x)
        z = self.reparameterize(mu, logvar)
        return self.decode(z), z, mu, logvar

    # Reconstruction + KL divergence losses summed over all elements and batch


def loss_function(recon_x, x, mu, logvar):
    BCE = F.binary_cross_entropy(recon_x, x, reduction='sum')

    # see Appendix B from VAE paper:
    # Kingma and Welling. Auto-Encoding Variational Bayes. ICLR, 2014
    # https://arxiv.org/abs/1312.6114
    # 0.5 * sum(1 + log(sigma^2) - mu^2 - sigma^2)
    KLD = -0.5 * torch.sum(1 + logvar - mu.pow(2) - logvar.exp())

    return BCE + KLD


csv_files = []
datasets = ['LUSC', 'LUAD', 'MESO', 'HNSC', 'BRCA', 'PRAD', 'SKCM', 'UVM', 'KIRP', 'KICH', 'KIRC', 'GBM', 'LGG', 'STAD',
            'PAAD']
for cur_ds in datasets:
    dataset = cur_ds
    constants.update_dirs(DATASET_NAME_u=dataset)
    data_normalizaton = "fpkm"
    gene_expression_file_name, phenotype_file_name, survival_file_name, mutation_file_name, mirna_file_name, pval_preprocessing_file_name = \
        build_gdc_params(dataset=dataset, data_normalizaton=data_normalizaton)
    csv_files.append(os.path.join(constants.DATA_DIR, gene_expression_file_name))

trainset = CancerTypesDataset(csv_files=csv_files, labels=datasets)
trainloader = torch.utils.data.DataLoader(trainset, batch_size=10,
                                          shuffle=True, num_workers=50, pin_memory=True)
testset = CancerTypesDataset(csv_files=csv_files, labels=datasets)
testloader = torch.utils.data.DataLoader(trainset, batch_size=10,
                                         shuffle=True, num_workers=50)

net = Net()
criterion = nn.BCELoss()

# create your optimizer
optimizer = optim.Adam(net.parameters(), lr=0.00001)

for epoch in range(100000):  # loop over the dataset multiple times

    running_loss = 0.0
    for i, data in enumerate(trainloader, 0):

        # get the inputs
        inputs, labels = data

        # zero the parameter gradients
        optimizer.zero_grad()

        # forward + backward + optimize
        outputs, z, mu, var = net(inputs)

        loss = loss_function(outputs, inputs, mu, var)
        loss.backward()
        running_loss += loss.item()
        optimizer.step()

        # print statistics
        if i % 100 == 99:  # print every 2000 mini-batches
            print('[%d, %5d] loss: %.3f' %
                  (epoch + 1, i + 1, running_loss / 100))
            running_loss = 0.0

    torch.save(net.state_dict(), os.path.join(constants.OUTPUT_GLOBAL_DIR, "VAE_model"))
    ###########################

    if epoch % 100 == 0:
        correct = 0
        total = 0
        X = None
        X_z = None
        X_mu = None
        X_var = None
        y = []

        with torch.no_grad():
            for data in testloader:
                features, labels = data
                _, labels = torch.max(labels, 1)
                outputs, z, mu, var = net(features)
                X_z = np.append(X_z, z, axis=0) if X_z is not None else z  
                X_mu=np.append(X_mu, mu, axis=0) if X_mu is not None else mu
                X_var=np.append(X_var, var, axis=0) if X_var is not None else var
                y = np.append(y, labels)

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
                   zip(datasets, colorlist_unique)]
        ax.legend(handles=patches)

        plt.savefig(
            os.path.join(constants.BASE_PROFILE, "output", "AE_by_samples_z_{}.png".format(epoch)))

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
                   zip(datasets, colorlist_unique)]
        ax.legend(handles=patches)

        plt.savefig(
            os.path.join(constants.BASE_PROFILE, "output", "AE_by_samples_mu_{}.png".format(epoch)))

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
                   zip(datasets, colorlist_unique)]
        ax.legend(handles=patches)

        plt.savefig(
            os.path.join(constants.BASE_PROFILE, "output", "AE_by_samples_logvar_{}.png".format(epoch)))

    ###########################

correct = 0
total = 0
X = None
X_r = None
y = []
with torch.no_grad():
    for data in testloader:
        features, labels = data
        X_r = np.append(X_r, features, axis=0) if X is not None else features
        _, labels = torch.max(labels, 1)
        outputs, z, mu, var = net(features)
        X = np.append(X, z, axis=0) if X is not None else z 
        y = np.append(y, labels)
        _, predicted = torch.max(outputs, 1)
        total += labels.size(0)
        correct += (predicted == labels).sum().item()

print('Accuracy of the network on  test batch: %d %%' % (
        100 * correct / total))


n_components=2
X = PCA(n_components=n_components).fit_transform(X_r)
fig = plt.figure(1, figsize=(20, 20))
plt.clf()
if n_components == 3:
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(X[:, 0], X[:, 1], X[:, 2], c=y, cmap='jet')
if n_components == 2:
    ax = fig.add_subplot(111)
    ax.scatter(X[:, 0], X[:, 1], c=y, cmap='jet')

colormap = cm.jet
label_ids_unique = np.unique(y)
label_ids = y

colorlist_unique = [ml_colors.rgb2hex(colormap(a)) for a in
                    label_ids_unique / float(max(label_ids))]
patches = [Line2D([0], [0], marker='o', color='gray', label=a,
                  markerfacecolor=c) for a, c in
           zip(datasets, colorlist_unique)]
ax.legend(handles=patches)

plt.savefig(
    os.path.join(constants.BASE_PROFILE, "output", "PCA_by_samples.png").format(constants.CANCER_TYPE))



n_components=2
X = TSNE(n_components=n_components, metric="correlation", perplexity=30.0).fit_transform(X_r)
fig = plt.figure(1, figsize=(20, 20))
plt.clf()
if n_components == 3:
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(X[:, 0], X[:, 1], X[:, 2], c=y, cmap='jet')
if n_components == 2:
    ax = fig.add_subplot(111)
    ax.scatter(X[:, 0], X[:, 1], c=y, cmap='jet')

colormap = cm.jet
label_ids_unique = np.unique(y)
label_ids = y

colorlist_unique = [ml_colors.rgb2hex(colormap(a)) for a in
                    label_ids_unique / float(max(label_ids))]
patches = [Line2D([0], [0], marker='o', color='gray', label=a,
                  markerfacecolor=c) for a, c in
           zip(datasets, colorlist_unique)]
ax.legend(handles=patches)

plt.savefig(
    os.path.join(constants.BASE_PROFILE, "output", "TSNE_by_samples.png").format(constants.CANCER_TYPE))
