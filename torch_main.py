import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
import os
import numpy as np
import constants
from cancer_type_dataset import CancerTypesDataset
import simplejson as json
from utils.param_builder import build_gdc_params
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

import matplotlib.cm as cm
import matplotlib.colors as ml_colors

from matplotlib.lines import Line2D

class Net(nn.Module):

    def __init__(self, factor = 0.5, n_mito_input_layer=100, n_cancer_types=2, n_latent_vector=2):
        super(Net, self).__init__()
        # self.factor = factor
        # self.n_mito_input_layer=n_mito_input_layer

        self.fc1 = nn.Linear(n_mito_input_layer, int(n_mito_input_layer*factor))
        self.fc2 = nn.Linear(int(n_mito_input_layer*factor), int(n_mito_input_layer*factor**2))
        self.fc3 = nn.Linear(int(n_mito_input_layer*factor**2), n_latent_vector)
        self.fc_out = nn.Linear(n_latent_vector, n_cancer_types)
    def forward(self, x):
        # Max pooling over a (2, 2) window
        # x = F.max_pool2d(F.relu(self.conv1(x)), (2, 2))
        # # If the size is a square you can only specify a single number
        # x = F.max_pool2d(F.relu(self.conv2(x)), 2)
        # x = x.view(-1, self.num_flat_features(x))
        h1 = F.relu(self.fc1(x))
        h2 = F.relu(self.fc2(h1))
        h3 = F.relu(self.fc3(h2))
        out = F.relu(self.fc_out(h3))
        return out, h3

csv_files=[]
datasets=['LUSC', 'LUAD']
for cur_ds in datasets:
    dataset=cur_ds
    constants.update_dirs(DATASET_NAME_u=dataset)
    data_normalizaton = "fpkm"
    gene_expression_file_name, phenotype_file_name, survival_file_name, mutation_file_name, mirna_file_name, pval_preprocessing_file_name = \
        build_gdc_params(dataset=dataset, data_normalizaton=data_normalizaton)
    csv_files.append(os.path.join(constants.DATA_DIR,gene_expression_file_name))


trainset = CancerTypesDataset(csv_files=csv_files, labels=datasets)
trainloader = torch.utils.data.DataLoader(trainset, batch_size=10,
                                          shuffle=True, num_workers=4, pin_memory=True)
testset = CancerTypesDataset(csv_files=csv_files, labels=datasets)
testloader = torch.utils.data.DataLoader(trainset, batch_size=10,
                                          shuffle=True, num_workers=4)



net = Net()
criterion = nn.CrossEntropyLoss()


# create your optimizer
optimizer = optim.SGD(net.parameters(), lr=0.0001)


for epoch in range(1000):  # loop over the dataset multiple times

    running_loss = 0.0
    for i, data in enumerate(trainloader, 0):

        # get the inputs
        inputs, labels = data

        # zero the parameter gradients
        optimizer.zero_grad()

        # forward + backward + optimize
        outputs, h1 = net(inputs)
        # labels = labels.squeeze_()
        labels = torch.max(labels, 1)[1]
        loss = criterion(outputs, labels)
        loss.backward()
        optimizer.step()

        # print statistics
        running_loss += loss.item()
        if i % 100 == 99:  # print every 2000 mini-batches
            print('[%d, %5d] loss: %.3f' %
                  (epoch + 1, i + 1, running_loss / 100))
            running_loss = 0.0

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
            outputs, h = net(features)
            X = np.append(X, h, axis=0) if X is not None else h
            y = np.append(y, labels)
            _, predicted = torch.max(outputs, 1)
            total += labels.size(0)
            correct += (predicted == labels).sum().item()

    print('Accuracy of the network on  test batch: %d %%' % (
            100 * correct / total))

    n_components=2
    fig = plt.figure(1, figsize=(20, 20))
    plt.clf()
    if n_components == 3:
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(X[:, 0], X[:, 1], X[:, 2], c=y, cmap='jet')
    if n_components == 2:
        ax = fig.add_subplot(111)
        ax.scatter(X[:, 0], X[:, 1], c=y, cmap='jet')

    colormap = cm.jet
    label_ids_unique=np.unique(y)
    label_ids=np.max(y)

    # colorlist_unique = [ml_colors.rgb2hex(colormap(i)) for i in
    #                     label_ids_unique / float(max(label_ids))]
    # patches = [Line2D([0], [0], marker='o', color='gray', label=a,
    #                   markerfacecolor=c) for a, c in
    #            zip(labels, colorlist_unique)]
    # ax.legend(handles=patches)

    plt.savefig(
        os.path.join(constants.BASE_PROFILE, "output", "AE_by_samples_{}.png".format(epoch)))

    ###########################



correct = 0
total = 0
X=None
X_r=None
y=[]
with torch.no_grad():
    for data in testloader:
        features, labels = data
        X_r=np.append(X_r,features,axis=0) if X is not None else features
        _, labels = torch.max(labels, 1)
        outputs, h = net(features)
        X=np.append(X,h,axis=0) if X is not None else h
        y=np.append(y, labels)
        _, predicted = torch.max(outputs,1)
        total += labels.size(0)
        correct += (predicted == labels).sum().item()

print('Accuracy of the network on  test batch: %d %%' % (
    100 * correct / total))

# n_components=2
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
# label_ids_unique=np.unique(y)
# label_ids=np.max(y)
#
# # colorlist_unique = [ml_colors.rgb2hex(colormap(i)) for i in
# #                     label_ids_unique / float(max(label_ids))]
# # patches = [Line2D([0], [0], marker='o', color='gray', label=a,
# #                   markerfacecolor=c) for a, c in
# #            zip(labels, colorlist_unique)]
# # ax.legend(handles=patches)
#
# plt.savefig(
#     os.path.join(constants.BASE_PROFILE, "output", "AE_by_samples.png").format(constants.CANCER_TYPE))
#
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
# label_ids_unique=np.unique(y)
# label_ids=np.max(y)
#
# # colorlist_unique = [ml_colors.rgb2hex(colormap(i)) for i in
# #                     label_ids_unique / float(max(label_ids))]
# # patches = [Line2D([0], [0], marker='o', color='gray', label=a,
# #                   markerfacecolor=c) for a, c in
# #            zip(labels, colorlist_unique)]
# # ax.legend(handles=patches)
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
# label_ids_unique=np.unique(y)
# label_ids=np.max(y)
#
# # colorlist_unique = [ml_colors.rgb2hex(colormap(i)) for i in
# #                     label_ids_unique / float(max(label_ids))]
# # patches = [Line2D([0], [0], marker='o', color='gray', label=a,
# #                   markerfacecolor=c) for a, c in
# #            zip(labels, colorlist_unique)]
# # ax.legend(handles=patches)
#
# plt.savefig(
#     os.path.join(constants.BASE_PROFILE, "output", "TSNE_by_samples.png").format(constants.CANCER_TYPE))
