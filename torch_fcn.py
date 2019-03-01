import torch
import torch.nn as nn
from torch.nn import functional as F
import torch.optim as optim
import os
import numpy as np
import constants
from torch_dataset_cancer import CancerTypesDataset
import torch_dataset_cancer
import simplejson as json
from utils.param_builder import build_gdc_params

from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as ml_colors
import fcn_bn_after_relu_flex_model

from matplotlib.lines import Line2D



num_workers=25
batch_size_train=100
batch_size_val=10


datasets=torch_dataset_cancer.CANCER_TYPES
torch_dataset=CancerTypesDataset(dataset_names=torch_dataset_cancer.CANCER_TYPES, meta_groups_files=torch_dataset_cancer.META_GROUPS, metagroups_names=["{}_{}".format(x.split("/")[1].split(".")[0], i_x) for i_x, x in enumerate(torch_dataset_cancer.META_GROUPS)])
train_dataset,test_dataset = torch.utils.data.random_split(torch_dataset, [torch_dataset.__len__()-torch_dataset.__len__()/100, torch_dataset.__len__()/100])

print "train: {}, test: {}".format(len(train_dataset), len(test_dataset))

trainloader = torch.utils.data.DataLoader(train_dataset, batch_size=batch_size_train,
                                          shuffle=True, num_workers=num_workers, pin_memory=True)

testloader = torch.utils.data.DataLoader(test_dataset, batch_size=batch_size_val,
                                          shuffle=True, num_workers=num_workers, pin_memory=True)

net = fcn_bn_after_relu_flex_model.Net(n_reduction_layers=3 ,factor=0.25)
load_model=True # False
PATH=os.path.join(constants.OUTPUT_GLOBAL_DIR, "FCN_model")
if load_model and os.path.exists(PATH):
   net.load_state_dict(torch.load(PATH))
   net.eval()

criterion = nn.MSELoss()

# create your optimizer
optimizer = optim.Adam(net.parameters(), lr=0.0001)

for epoch in range(0, 100000):  # loop over the dataset multiple times

    train_loss = 0.0
    val_loss = 0.0

    for i, data in enumerate(trainloader, 0):

        # get the inputs
        inputs, labels = data

        # zero the parameter gradients
        optimizer.zero_grad()

        # forward + backward + optimize
        outputs, h = net(inputs)

        loss = criterion(outputs, labels)
        loss.backward()
        train_loss += loss.item()
        optimizer.step()

        # print statistics
        if i % 10 == 9:  # print every 2000 mini-batches
            print('[%d, %5d] train loss: %.3f' %
                  (epoch + 1, i + 1, train_loss / 100))
            train_loss = 0.0

        torch.save(net.state_dict(), os.path.join(constants.OUTPUT_GLOBAL_DIR, "FCN_model"))

    for i, data in enumerate(testloader, 0):
        with torch.no_grad():
            # get the inputs
            inputs, labels = data
    
            # forward + backward + optimize
            outputs, h = net(inputs)
    
            loss = criterion(outputs, labels)
            val_loss += loss.item()
    
    # print statistics
    
    print('[%d, %5d] val loss: %.3f' %
          (epoch + 1, i + 1, val_loss / 100))
    val_loss = 0.0

    ###########################

    if epoch % 100==0: 
        correct = 0
        total = 0
        X = None
        X_z = None
        X_mu = None
        X_var = None
        y = []

        with torch.no_grad():
            for i, data in enumerate(testloader, 0):
                features, labels = data
                _, labels = torch.max(labels, 1)
                outputs, z = net(features)
                X_z = np.append(X_z, z, axis=0) if X_z is not None else z
                y = np.append(y, labels)

        print "len samples: {}".format(len(X_z))
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
            os.path.join(constants.BASE_PROFILE, "output", "FCN_by_samples_{}.png".format(epoch)))

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
