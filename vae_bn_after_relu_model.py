import torch
import torch.nn as nn
from torch.nn import functional as F


class Net(nn.Module):

    def __init__(self, factor=0.25, n_mito_input_layer=2000, n_cancer_types=2, n_latent_vector=2):
        super(Net, self).__init__()
        # self.factor = factor
        # self.n_mito_input_layer=n_mito_input_layer

        self.fc1 = nn.Linear(n_mito_input_layer, int(n_mito_input_layer * factor))
        self.fc1_bn = nn.BatchNorm1d(int(n_mito_input_layer * factor))
        self.fc2 = nn.Linear(int(n_mito_input_layer * factor), int(n_mito_input_layer * factor ** 2))
        self.fc2_bn = nn.BatchNorm1d(int(n_mito_input_layer * factor ** 2))
        self.fc31 = nn.Linear(int(n_mito_input_layer * factor ** 2), n_latent_vector)
        self.fc31_bn = nn.BatchNorm1d(n_latent_vector)
        self.fc32 = nn.Linear(int(n_mito_input_layer * factor ** 2), n_latent_vector)
        self.fc32_bn = nn.BatchNorm1d(n_latent_vector)
        self.fc4 = nn.Linear(n_latent_vector, int(n_mito_input_layer * factor ** 2))
        self.fc4_bn = nn.BatchNorm1d(int(n_mito_input_layer * factor ** 2))
        self.fc5 = nn.Linear(int(n_mito_input_layer * factor ** 2), int(n_mito_input_layer * factor))
        self.fc5_bn = nn.BatchNorm1d(int(n_mito_input_layer * factor))
        self.fc6 = nn.Linear(int(n_mito_input_layer * factor), n_mito_input_layer)

    def encode(self, x):
        h1 = self.fc1_bn(F.relu(self.fc1(x)))
        h2 = self.fc2_bn(F.relu(self.fc2(h1)))
        h31 = self.fc31_bn(self.fc31(h2))
        h32 = self.fc32_bn(self.fc32(h2))
        return h31, h32

    def reparameterize(self, mu, logvar):
        std = torch.exp(0.5 * logvar)
        eps = torch.randn_like(std)
        return eps.mul(std).add_(mu)

    def decode(self, z):
        h4 = self.fc4_bn(F.relu(self.fc4(z)))
        h5 = self.fc5_bn(F.relu(self.fc5(h4)))
        h6 = F.sigmoid(self.fc6(h5))
        return h6

    def forward(self, x):
        mu, logvar = self.encode(x)
        z = self.reparameterize(mu, logvar)
        return self.decode(z), z, mu, logvar
