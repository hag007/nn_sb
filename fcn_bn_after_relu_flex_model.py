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

        for cur in np.arange(1, n_reduction_layers):
            setattr(self, "fc_e"+str(cur), nn.Linear(int(n_mito_input_layer* factor ** (cur-1)), int(n_mito_input_layer * factor ** cur)))
            setattr(self,"fc_bn_e"+str(cur), nn.BatchNorm1d(int(n_mito_input_layer * factor ** cur)))

        self.fc_l = nn.Linear(int(n_mito_input_layer * factor ** n_reduction_layers), n_latent_vector)
        self.fc_l_bn = nn.BatchNorm1d(n_latent_vector)
        self.fc_out = nn.Linear(n_latent_vector, len(cancer_type_dataset.CANCER_TYPES))


    def forward(self, x):

        h = x
        for cur in np.arange(1, self.n_reduction_layers):
            h = getattr(self, "fc_bn_e" + str(cur))(F.relu(getattr(self, "fc_e" + str(cur))(h)))

        h=self.fc_l_bn(F.relu(self.fc_l(h)))
        out = F.softmax(self.fc_out(h))

        return out, h
