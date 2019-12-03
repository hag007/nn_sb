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

class Encoder(nn.Module):
    def __init__(self, factor=0.5, n_mito_input_layer=torch_dataset_cancer.n_input_layer, n_cancer_types=2, n_latent_vector=100, n_reduction_layers=2):
        super(Encoder, self).__init__()
        self.n_reduction_layers = n_reduction_layers
        self.n_latent_vector = n_latent_vector

        for cur in np.arange(1, n_reduction_layers + 1):
            setattr(self, "fc_enc" + str(cur),
                    nn.Linear(int(n_mito_input_layer * factor ** (cur - 1)), int(n_mito_input_layer * factor ** cur)))

        self.fc_enc_l_mu = nn.Linear(int(n_mito_input_layer * factor ** n_reduction_layers), n_latent_vector)
        self.fc_enc_l_var = nn.Linear(int(n_mito_input_layer * factor ** n_reduction_layers), n_latent_vector)

    def encode(self, x):
        h=x
        for cur in np.arange(1, self.n_reduction_layers + 1):
            h=F.relu(getattr(self, "fc_enc" + str(cur))(h))

        l_mu = F.relu(getattr(self, "fc_enc_l_mu")(h))
        l_var = F.relu(getattr(self, "fc_enc_l_var")(h))

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

    def __init__(self, factor=0.5, n_mito_input_layer=torch_dataset_cancer.n_input_layer, n_cancer_types=2, n_latent_vector=100, n_reduction_layers=2):
        super(Decoder, self).__init__()
        self.n_reduction_layers = n_reduction_layers
        self.n_latent_vector= n_latent_vector

        self.fc_dec_l = nn.Linear(n_latent_vector, int(n_mito_input_layer * factor ** n_reduction_layers))
        self.fc_bn_dec_l = nn.BatchNorm1d(int(n_mito_input_layer * factor ** n_reduction_layers))

        for cur in np.arange(n_reduction_layers, 1, -1):
            setattr(self, "fc_dec"+str(cur), nn.Linear(int(n_mito_input_layer* factor ** cur), int(n_mito_input_layer * factor ** (cur-1))))

        setattr(self, "fc_dec1",
                nn.Linear(int(n_mito_input_layer * factor), int(n_mito_input_layer)))



    def decode(self, z):

        h = F.relu(getattr(self, "fc_dec_l")(z))
        for cur in np.arange(self.n_reduction_layers, 1, -1):
            h = F.relu(getattr(self,"fc_dec" + str(cur))(h))

        h = F.sigmoid(getattr(self,"fc_dec1")(h))

        return h



    def forward(self, input):
        z= input
        
        decoded=self.decode(z)
        return decoded


class VAE_GAN_Generator(nn.Module):
    def __init__(self, factor=0.5, n_mito_input_layer=torch_dataset_cancer.n_input_layer, n_cancer_types=2, n_latent_vector=100, n_reduction_layers=2):
        super(VAE_GAN_Generator, self).__init__()
        self.n_reduction_layers = n_reduction_layers
        self.n_latent_vector= n_latent_vector

        self.encoder = Encoder(factor=0.5, n_mito_input_layer=torch_dataset_cancer.n_input_layer, n_cancer_types=2, n_latent_vector=n_latent_vector, n_reduction_layers=2)
        self.decoder = Decoder(factor=0.5, n_mito_input_layer=torch_dataset_cancer.n_input_layer, n_cancer_types=2, n_latent_vector=n_latent_vector, n_reduction_layers=2)

    def forward(self, x):
        z = self.encoder(x)
        rec_images = self.decoder(z)

        return rec_images


class Discriminator(nn.Module):

    def __init__(self, factor=0.5, n_mito_input_layer=torch_dataset_cancer.n_input_layer, n_cancer_types=2, n_latent_vector=100, n_reduction_layers=2):
        super(Discriminator, self).__init__()
        # self.factor = factor
        # self.n_mito_input_layer=n_mito_input_layer
        self.n_reduction_layers = n_reduction_layers
        self.n_latent_vector = n_latent_vector

        for cur in np.arange(1, n_reduction_layers + 1):
            setattr(self, "fc_dis" + str(cur),
                    nn.Linear(int(n_mito_input_layer * factor ** (cur - 1)), int(n_mito_input_layer * factor ** cur)))

        self.fc_dis_l = nn.Linear(int(n_mito_input_layer * factor ** n_reduction_layers), n_latent_vector)
        self.fc_bn_dis_l = nn.BatchNorm1d(n_latent_vector)
        self.fc_out = nn.Linear(int(n_latent_vector), 1)

    def discriminate(self, x_hat):

        h = x_hat
        for cur in np.arange(1, self.n_reduction_layers + 1):
            h = F.relu(getattr(self, "fc_dis" + str(cur))(h))

        l=F.sigmoid(getattr(self, "fc_dis_l")(h))


        out_dis = F.sigmoid(self.fc_out(getattr(self, "fc_dis_l")(l)))

        return out_dis, l

    def forward(self, input):
        encoded = input
        dis_prediction, l =self.discriminate(encoded)
        return dis_prediction, l


    def similarity(self, x):
        batch_size = x.size()[0]
        outputs, features = self.discriminate(x)
        return features

