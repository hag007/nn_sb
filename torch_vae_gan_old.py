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
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import matplotlib.cm as cm
import matplotlib.colors as ml_colors

from matplotlib.lines import Line2D

from cancer_type_dataset import CANCER_TYPES

EPOCHS=10

class Encoder(nn.Module):

    def __init__(self, factor = 0.5, n_mito_input_layer=100, n_cancer_types=2, n_latent_vector=2):
        super(Encoder, self).__init__()
        # self.factor = factor
        # self.n_mito_input_layer=n_mito_input_layer

        self.fc1 = nn.Linear(n_mito_input_layer, int(n_mito_input_layer*factor))
        self.fc2 = nn.Linear(int(n_mito_input_layer*factor), int(n_mito_input_layer*factor**2))
        self.fc31 = nn.Linear(int(n_mito_input_layer*factor**2), n_latent_vector)
        self.fc32 = nn.Linear(int(n_mito_input_layer * factor ** 2), n_latent_vector)

    def encode(self, x):
        h_en1 = F.relu(self.fc1(x))
        h_en2 = F.relu(self.fc2(h_en1))
        h_en31 = self.fc31(h_en2)
        h_en32 = self.fc32(h_en2)

        return h_en31, h_en32

    def reparameterize(self, mu, logvar):
        std = torch.exp(0.5 * logvar)
        eps = torch.randn_like(std)

        return eps.mul(std).add_(mu)


    def forward(self, x):
        mu, logvar = self.encode(x)
        z = self.reparameterize(mu, logvar)
        return z, mu, logvar


class Decoder(nn.Module):

    def __init__(self, factor = 0.5, n_mito_input_layer=100, n_cancer_types=2, n_latent_vector=2):
        super(Decoder, self).__init__()

        self.fc_dec1 = nn.Linear(n_latent_vector, int(n_mito_input_layer * factor ** 2))
        self.fc_dec2 = nn.Linear(int(n_mito_input_layer*factor**2), int(n_mito_input_layer*factor))
        self.fc_dec3 = nn.Linear(int(n_mito_input_layer*factor), n_mito_input_layer)


    def decode(self, z):
        h_dec1 = F.relu(self.fc_dec1(z))
        h_dec2 = F.relu(self.fc_dec2(h_dec1))
        out_dec = F.sigmoid(self.fc_dec3(h_dec2))

        return out_dec


    def forward(self, input):
        z, mu, logvar = input
        decoded=self.decode(z)
        return decoded ,mu, logvar


class Discriminator(nn.Module):

    def __init__(self, factor = 0.5, n_mito_input_layer=100, n_cancer_types=2, n_latent_vector=2):
        super(Discriminator, self).__init__()
        # self.factor = factor
        # self.n_mito_input_layer=n_mito_input_layer

        self.fc_dis1 = nn.Linear(n_mito_input_layer, int(n_mito_input_layer*factor))
        self.fc_dis2 = nn.Linear(int(n_mito_input_layer*factor), int(n_mito_input_layer*factor**2))
        self.fc_dis3 = nn.Linear(int(n_mito_input_layer*factor**2), 1)


    def discriminator(self, x_hat):
        h_dis1 = F.relu(self.fc_dis1(x_hat))
        h_dis2 = F.relu(self.fc_dis2(h_dis1))
        out_dis = F.sigmoid(self.fc_dis3(h_dis2))

        return out_dis

    def forward(self, input):
        x_hat, _1, _2 = input
        dis_prediction=self.discriminator(x_hat)
        return dis_prediction


# Reconstruction + KL divergence losses summed over all elements and batch
def loss_function(recon_x, x, mu, logvar):
    BCE = F.binary_cross_entropy(recon_x, x, reduction='sum')

    # see Appendix B from VAE paper:
    # Kingma and Welling. Auto-Encoding Variational Bayes. ICLR, 2014
    # https://arxiv.org/abs/1312.6114
    # 0.5 * sum(1 + log(sigma^2) - mu^2 - sigma^2)
    KLD = -0.5 * torch.sum(1 + logvar - mu.pow(2) - logvar.exp())

    return BCE + KLD


def train_VAE(trainloader, VAE, optimizer):

    for k, v in VAE.state_dict().iteritems():
        if k.startswith("_"): continue

        v.requires_grad = True


    for epoch in range(EPOCHS):  # loop over the dataset multiple times
        print "EPOCH: {}".format(epoch)
        running_loss = 0.0
        for i, data in enumerate(trainloader, 0):

            # get the inputs
            inputs, labels = data

            # zero the parameter gradients
            optimizer.zero_grad()

            # forward + backward + optimize
            outputs, mu, logvar = VAE(inputs)

            loss = loss_function(outputs, inputs, mu, logvar)
            loss.backward()
            running_loss += loss.item()
            optimizer.step()


def train_GAN_dis(trainloader, GAN, optimizer, m_decoder, m_discriminator):

        for k, v in GAN.state_dict().iteritems():
            if k.startswith("_"): continue

            if "_dec" in k:
                v.requires_grad = False
            else:
                v.requires_grad = True

        batch_accuracies=[]
        for epoch in range(EPOCHS):  # loop over the dataset multiple times

            running_loss = 0.0
            for i, data in enumerate(trainloader, 0):

                # get the inputs
                inputs, labels = data
                inputs=torch.tensor(torch.stack([inputs[a] for a in np.arange(len(inputs)/2)]), dtype=torch.float)
                labels=torch.tensor(torch.tensor([1 for a in np.arange(len(labels)/2)]), dtype=torch.float)
                random_inputs=[]
                random_labels=[]
                with torch.no_grad():
                    for input in inputs:

                        random_labels.append(0)
                        outputs, _1, _2=m_decoder([torch.distributions.normal.Normal(loc=0,scale=1).sample(sample_shape=(2,)), None, None])
                        random_inputs.append(outputs)


                inputs=torch.cat((inputs, torch.stack(random_inputs)), dim=0)
                # inputs=torch.tensor(torch.stack(random_inputs), dtype=torch.float)
                labels=torch.cat((torch.tensor(labels), torch.tensor(random_labels, dtype=torch.float)), dim=0)
                # labels= torch.tensor(random_labels, dtype=torch.float)
                # zero the parameter gradients
                optimizer.zero_grad()

                # forward + backward + optimize
                outputs = m_discriminator([inputs,None,None])

                loss = F.binary_cross_entropy(outputs.view(-1,), torch.tensor(labels, dtype=torch.float), reduction='sum')
                loss.backward()
                running_loss += loss.item()
                optimizer.step()

                if epoch==EPOCHS-1:
                    correct = 0
                    total = 0
                    y = []
                    with torch.no_grad():
                        y_pred=m_discriminator([inputs,None,None])
                        y_pred=torch.round(y_pred).view(-1)
                        y = torch.tensor(labels,dtype=torch.float)
                        total += labels.size(0)
                        correct += (y_pred == y).sum().item()
                        batch_accuracies.append(100 * correct / total)

        print('Accuracy of the network on  test batch: {} %'.format(
            round(np.mean(np.array(batch_accuracies)),2)))


def train_GAN_gen(trainloader, GAN, optimizer):

    for k, v in GAN.state_dict().iteritems():
        if k.startswith("_"): continue
        if "_dis" in k:
            v.requires_grad = False
        else:
            v.requires_grad = True

    for epoch in range(EPOCHS):  # loop over the dataset multiple times

        running_loss = 0.0
        for i, data in enumerate(trainloader, 0):

            # get the inputs
            inputs, labels = data

            random_labels=[]
            random_inputs=[]
            for i in range(inputs.size()[0]*2):
                random_labels.append(1)
                random_inputs.append(torch.distributions.normal.Normal(loc=0, scale=1).sample(sample_shape=(2,)))

            # zero the parameter gradients
            optimizer.zero_grad()

            # forward + backward + optimize
            outputs = GAN([torch.stack(random_inputs), None, None])

            loss = F.binary_cross_entropy(outputs.view(-1, ), torch.tensor(random_labels, dtype=torch.float), reduction='sum')
            loss.backward()
            running_loss += loss.item()
            optimizer.step()



encoder=Encoder()
decoder=Decoder()
discriminator=Discriminator()

m_VAE = nn.Sequential(encoder,decoder)
m_GAN = nn.Sequential(decoder,discriminator)
m_FULL = nn.Sequential(encoder,decoder,discriminator)

csv_files = []
datasets = CANCER_TYPES

for cur_ds in datasets:
    dataset = cur_ds
    constants.update_dirs(DATASET_NAME_u=dataset)
    data_normalizaton = "fpkm"
    gene_expression_file_name, phenotype_file_name, survival_file_name, mutation_file_name, mirna_file_name, pval_preprocessing_file_name = \
        build_gdc_params(dataset=dataset, data_normalizaton=data_normalizaton)
    csv_files.append(os.path.join(constants.DATA_DIR, gene_expression_file_name))

trainset = CancerTypesDataset(csv_files=csv_files, labels=datasets)
trainloader = torch.utils.data.DataLoader(trainset, batch_size=100,
                                          shuffle=True, num_workers=40, pin_memory=True)
testset = trainset # CancerTypesDataset(csv_files=csv_files, labels=datasets)
testloader = trainloader # torch.utils.data.DataLoader(trainset, batch_size=10,
                                         # shuffle=True, num_workers=10)

criterion = nn.BCELoss()

# create your optimizer
vae_optimizer = optim.Adam(m_VAE.parameters(), lr=0.00001)
gan_optimizer = optim.Adam(m_GAN.parameters(), lr=0.00001)


for meta_epoch in range(3000):

    print "meta_epoch: {}".format(meta_epoch)

    ###########################
    if meta_epoch % 1 == 0:
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
                outputs, h, var = m_VAE(features)
                X = np.append(X, h, axis=0) if X is not None else h
                y = np.append(y, labels)

        n_components = 2
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
            os.path.join(constants.BASE_PROFILE, "output", "AE_by_samples_{}.png".format(meta_epoch)))

    ###########################

    # if meta_epoch%3==0:
    train_VAE(trainloader, m_VAE, vae_optimizer)
    # elif meta_epoch%3==1:
    train_GAN_dis(trainloader, m_GAN, gan_optimizer, decoder, discriminator)
    # else:
    train_GAN_gen(trainloader, m_GAN, gan_optimizer)

    torch.save(encoder.state_dict(), os.path.join(constants.OUTPUT_GLOBAL_DIR, "encoder_model"))
    torch.save(decoder.state_dict(), os.path.join(constants.OUTPUT_GLOBAL_DIR, "decoder_model"))
    torch.save(discriminator.state_dict(), os.path.join(constants.OUTPUT_GLOBAL_DIR, "discriminator_model"))




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
