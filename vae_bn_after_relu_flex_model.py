import torch
import torch.nn as nn
import cancer_type_dataset
from torch.nn import functional as F
import numpy as np

class Net(nn.Module):

    def __init__(self, factor=0.5, n_mito_input_layer=cancer_type_dataset.n_input_layer, n_cancer_types=2, n_latent_vector=2, n_reduction_layers=2):
        super(Net, self).__init__()
        # self.factor = factor
        # self.n_mito_input_layer=n_mito_input_layer
        self.n_reduction_layers=n_reduction_layers

        for cur in np.arange(1, n_reduction_layers+1):
            setattr(self, "fc_e"+str(cur), nn.Linear(int(n_mito_input_layer* factor ** (cur-1)), int(n_mito_input_layer * factor ** cur)))
            setattr(self,"fc_bn_e"+str(cur), nn.BatchNorm1d(int(n_mito_input_layer * factor ** cur)))

        self.fc_e_l_mu = nn.Linear(int(n_mito_input_layer * factor ** n_reduction_layers), n_latent_vector)
        self.fc_bn_e_l_mu = nn.BatchNorm1d(n_latent_vector)
        self.fc_e_l_var = nn.Linear(int(n_mito_input_layer * factor ** n_reduction_layers), n_latent_vector)
        self.fc_bn_e_l_var = nn.BatchNorm1d(n_latent_vector)

        self.fc_d_l = nn.Linear(n_latent_vector, int(n_mito_input_layer * factor ** n_reduction_layers))
        self.fc_bn_d_l = nn.BatchNorm1d(int(n_mito_input_layer * factor ** n_reduction_layers))

        for cur in np.arange(n_reduction_layers, 1, -1):
            setattr(self, "fc_d"+str(cur), nn.Linear(int(n_mito_input_layer* factor ** cur), int(n_mito_input_layer * factor ** (cur-1))))
            setattr(self, "fc_bn_d"+str(cur), nn.BatchNorm1d(int(n_mito_input_layer * factor ** (cur-1))))
        setattr(self, "fc_d1",
                nn.Linear(int(n_mito_input_layer * factor), int(n_mito_input_layer)))


    def encode(self, x):
        h=x
        for cur in np.arange(1, self.n_reduction_layers + 1):
            h=getattr(self, "fc_bn_e"+ str(cur))(F.relu(getattr(self, "fc_e" + str(cur))(h)))

        l_mu = getattr(self, "fc_bn_e_l_mu")(F.relu(getattr(self, "fc_e_l_mu")(h)))
        l_var = getattr(self, "fc_bn_e_l_var")(F.relu(getattr(self, "fc_e_l_var")(h)))

        return l_mu, l_var

    def reparameterize(self, mu, logvar):
        std = torch.exp(0.5 * logvar)
        eps = torch.randn_like(std)
        return eps.mul(std).add_(mu)

    def decode(self, z):

        h = getattr(self, "fc_bn_d_l")(F.relu(getattr(self, "fc_d_l")(z)))
        for cur in np.arange(self.n_reduction_layers, 1, -1):
            h = getattr(self,"fc_bn_d" + str(cur))(F.relu(getattr(self,"fc_d" + str(cur))(h)))

        h = F.sigmoid(getattr(self,"fc_d1")(h))

        return h

    def forward(self, x):
        mu, logvar = self.encode(x)
        z = self.reparameterize(mu, logvar)
        return self.decode(z), z, mu, logvar
