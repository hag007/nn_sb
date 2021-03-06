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
import torch_dataset_sc
from torch.nn import functional as F
import constants
import torch_sc_vae_gan_copy_model

batch_size_train = 100
batch_size_val = 10
num_workers = 25

datasets = torch_dataset_sc.DATASETS
torch_dataset = torch_dataset_sc.SingleCellDataset(dataset_names=torch_dataset_sc.DATASETS,
                                                   meta_groups_files=torch_dataset_sc.META_GROUPS,
                                                   metagroups_names=["{}_{}".format(x, i_x) for i_x, x in
                                                                     enumerate(torch_dataset_sc.DATASETS)])
train_dataset, test_dataset = torch.utils.data.random_split(torch_dataset,
                                                            [torch_dataset.__len__() - torch_dataset.__len__() / 100,
                                                             torch_dataset.__len__() / 100])

print "n_train samples: {}, n_test samples: {}".format(len(train_dataset), len(test_dataset))

trainloader = torch.utils.data.DataLoader(train_dataset, batch_size=batch_size_train, shuffle=True,
                                          num_workers=num_workers, pin_memory=True)

testloader = torch.utils.data.DataLoader(test_dataset, batch_size=batch_size_val, shuffle=True, num_workers=num_workers,
                                         pin_memory=True)

# define constant
max_epochs = 100000
lr = 0.00005

beta = 5
alpha = 0.1
gamma = 5
delta = 5
n_latent_vector = 2
G = torch_sc_vae_gan_copy_model.VAE_GAN_Generator(n_latent_vector=n_latent_vector)
D = torch_sc_vae_gan_copy_model.Discriminator(n_latent_vector=n_latent_vector)

criterion = nn.BCELoss(reduction='sum')
criterion

opt_enc = optim.Adam(G.encoder.parameters(), lr=lr)
opt_dec = optim.Adam(G.decoder.parameters(), lr=lr)
opt_dis = optim.Adam(D.parameters(), lr=lr * alpha)
opt_vae = optim.Adam(G.parameters(), lr=lr)

fixed_noise = Variable(torch.randn(batch_size_train, n_latent_vector))
data, _ = next(iter(torch_dataset))
fixed_batch = Variable(data)
min_val = 10000000
min_val_epoch = -1

for epoch in range(max_epochs):
    train_loss = 0
    val_loss = 0
    print "cur epoch: {}".format(epoch)
    D_real_list, D_rec_enc_list, D_rec_noise_list, D_list = [], [], [], []
    g_loss_list, rec_loss_list, prior_loss_list = [], [], []
    for i, data_tuple in enumerate(trainloader, 0):
        data, _ = data_tuple
        batch_size = data.size()[0]
        ones_label = Variable(torch.ones(batch_size))
        zeros_label = Variable(torch.zeros(batch_size))
        # print("batch_size", batch_size)
        # print (data.shape)
        datav = Variable(data)
        rec_enc = G(datav)
        mean, logvar = G.encoder.encode(datav)

        similarity_rec_enc = rec_enc  # l_rec
        similarity_data = datav  # l_real

        rec_loss = nn.BCELoss(reduction='sum')(similarity_rec_enc, similarity_data.detach())

        # train encoder
        prior_loss = 1 + logvar - mean.pow(2) - logvar.exp()
        prior_loss = (-0.5 * torch.sum(prior_loss))

        err_enc = prior_loss + beta * rec_loss

        train_loss += err_enc.item()

        if i % 10 == 9:  # print every 2000 mini-batches
            print('[%d, %5d] train loss: %.3f' %
                  (epoch + 1, i + 1, train_loss / 100))
            train_loss = 0.0

        opt_vae.zero_grad()
        err_enc.backward()

        opt_vae.step()

    for i, data in enumerate(testloader, 0):
        with torch.no_grad():
            # get the inputs
            inputs, labels = data

            # forward + backward + optimize
            rec_enc = G(datav)
            mean, logvar = G.encoder.encode(datav)

            prior_loss = 1 + logvar - mean.pow(2) - logvar.exp()
            prior_loss = (-0.5 * torch.sum(prior_loss))  #

            similarity_rec_enc = rec_enc  # l_rec
            similarity_data = datav  # l_real

            rec_loss = nn.BCELoss(reduction='sum')(similarity_rec_enc, similarity_data.detach())
            loss = prior_loss + rec_loss
            val_loss += loss.item()

    print('[%d, %5d] val loss: %.3f' % (epoch + 1, i + 1, val_loss / 100))
    if min_val > val_loss / 100:
        min_val = val_loss / 100
        min_val_epoch = epoch
        torch.save(G.encoder.state_dict(), os.path.join(constants.OUTPUT_GLOBAL_DIR, "VAE_ENC_SC_mdl"))
        torch.save(G.decoder.state_dict(), os.path.join(constants.OUTPUT_GLOBAL_DIR, "VAE_DEC_SC_mdl"))
        torch.save(D.state_dict(), os.path.join(constants.OUTPUT_GLOBAL_DIR, "VAE_DIS_SC_mdl"))
    print "min_val epoch: {}, min_val: {}".format(min_val_epoch, min_val)



