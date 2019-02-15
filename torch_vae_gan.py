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
import vae_gan_bn_after_relu_flex_model, vae_bn_after_relu_flex_model

from matplotlib.lines import Line2D



num_workers=25
batch_size_train=100
batch_size_val=10

def loss_function(recon_x, x, mu, logvar, KLD_ratio):
    BCE = F.mse_loss(recon_x, x, reduction='sum')

    # see Appendix B from VAE paper:
    # Kingma and Welling. Auto-Encoding Variational Bayes. ICLR, 2014
    # https://arxiv.org/abs/1312.6114
    # 0.5 * sum(1 + log(sigma^2) - mu^2 - sigma^2)
    KLD = -0.5 * torch.sum(1 + logvar - mu.pow(2) - logvar.exp())

    return BCE + KLD*KLD_ratio


def backprop_vae(m_VAE, m_FULL, optimizer, n_epoches, trainloader, testloader, KLD_ratio):
    
    
     for k, v in m_FULL.state_dict().iteritems():
         if k.startswith("_"): continue

         if "_dis" in k:
             v.requires_grad = False
         elif ("_enc" in k or "_dec" in k) and "Float" in torch.typename(v):
             v.requires_grad = True


     for epoch in range(n_epoches):
        train_loss=0
        val_loss=0
        for i, data in enumerate(trainloader, 0):

            # get the inputs
            inputs, labels = data

            # zero the parameter gradients
            optimizer.zero_grad()

            # forward + backward + optimize
            outputs, z, mu, var = m_VAE(inputs)
            auth_decoded, _1, _2, _3, _4, dis_l_decoded = discriminator([outputs, None, None, None])
            auth_real, _1, _2, _3, _4, dis_l_real = discriminator([inputs, None, None, None])

            
            loss = loss_function(dis_l_decoded, dis_l_real.detach(), mu, var, KLD_ratio)
            loss.backward()
            train_loss += loss.item()
            optimizer.step()

            # print statistics
            if i % 10 == 9:  # print every 2000 mini-batches
                print('[%d, %5d] vae train loss: %.3f' %
                      (epoch + 1, i + 1, train_loss / 100))
                train_loss = 0.0

        for i, data in enumerate(testloader, 0):
            with torch.no_grad():
                # get the inputs
                inputs, labels = data

                # forward + backward + optimize
                outputs, z, mu, var = m_VAE(inputs)
                
                auth_decoded, _1, _2, _3, _4, dis_l_decoded = discriminator([outputs, None, None, None])
                auth_real, _1, _2, _3, _4, dis_l_real = discriminator([inputs, None, None, None])

                loss = loss_function(dis_l_decoded, dis_l_real.detach(), mu, var, KLD_ratio)
                val_loss += loss.item()

        # print statistics

        print('[%d, %5d] vae val loss: %.3f' %
              (epoch + 1, i + 1, val_loss / 100))
        val_loss = 0.0

    
        for k, v in m_FULL.state_dict().iteritems():
           if k.startswith("_") or "Float" not in torch.typename(v): continue
           v.requires_grad = True
     


def backprop_dis(m_VAE, m_discriminator, optimizer, n_epoches, trainloader):

    accuracy=0
    for epoch in range(n_epoches):  # loop over the dataset multiple times
        batch_accuracies=[]
        running_loss = 0.0
        for i, data in enumerate(trainloader, 0):

            inputs, labels = data
            real_inputs = torch.tensor(torch.stack([inputs[a] for a in np.arange(len(inputs)/2)]),
                                  dtype=torch.float)
            real_labels = torch.tensor(torch.tensor([0.1 for a in np.arange(len(labels)/2, (len(labels)/2)*2)]), dtype=torch.float)

            dummy_inputs=torch.tensor(torch.stack([inputs[a] for a in np.arange(len(inputs) / 2)]), dtype=torch.float) # + torch.tensor(torch.stack([inputs[a] for a in np.arange(len(inputs) / 2, (len(inputs)/2)*2)]),
              #                            dtype=torch.float))/2.0
            dummy_labels = torch.tensor(torch.tensor([0.9 for a in np.arange(len(labels) / 2)]), dtype=torch.float)

            # dummy_labels = torch.tensor(torch.tensor([0 for a in np.arange(len(labels) / 2)]), dtype=torch.float)
            random_inputs, _1, _2, _3 = m_VAE(dummy_inputs)

            inputs = torch.cat((real_inputs, random_inputs), dim=0)
            labels = torch.cat((torch.tensor(real_labels), torch.tensor(dummy_labels, dtype=torch.float)), dim=0)

            optimizer.zero_grad()

            # forward + backward + optimize
            outputs, _1, _2, _3, _4, _5 = m_discriminator([inputs, None, None, None])

            loss = F.mse_loss(outputs.view(-1, ), torch.tensor(labels, dtype=torch.float),
                                          reduction='sum')
            loss.backward()
            running_loss += loss.item()
            optimizer.step()

            correct = 0
            total = 0
            y = []
            with torch.no_grad():
                y_pred_raw, _1, _2, _3, _4, _5 = m_discriminator([inputs, None, None, None])
                y_pred = torch.round(y_pred_raw).view(-1)
                y = torch.tensor(labels, dtype=torch.float)
                total += labels.size(0)
                correct += (y_pred == torch.round(y)).sum().item()
                # print y_pred_raw
                # print y_pred
                batch_accuracies.append(100 * correct / total)
       
        accuracy=round(np.mean(np.array(batch_accuracies)), 2)
        print('Accuracy of the network on  test batch: {} %. loss: {}'.format(accuracy, running_loss))

    return running_loss, accuracy

def backprop_gen(m_VAE, m_GAN, m_discriminator, optimizer, n_epoches, trainloader):

    for k, v in m_GAN.state_dict().iteritems():
            if k.startswith("_"): continue

            if "_dis" in k:
                v.requires_grad = False
            elif ("_enc" in k or "_dec" in k) and "Float" in torch.typename(v):
                v.requires_grad = True

    for epoch in range(n_epoches):
        running_loss = 0.0
        print "cur epoch: {}".format(epoch)
        for i, data in enumerate(trainloader, 0):

            # get the inputs
            inputs, labels = data

            random_inputs, z,  mu, var = m_VAE(inputs)
            random_labels = torch.tensor([0.1 for x in range(labels.size()[0])], dtype=torch.float)

            # zero the parameter gradients
            optimizer.zero_grad()

            # forward + backward + optimize
            outputs, _1, _2, _3, _4, _5 = m_discriminator([random_inputs, None, None, None])

            # print outputs 
            loss = F.mse_loss(outputs.view(-1, ), torch.tensor(random_labels, dtype=torch.float),
                                          reduction='sum')
          
            loss.backward()
            running_loss += loss.item()
            optimizer.step()

        # print statistics
        print('[%d, %5d] gen train loss: %.3f' %
              (epoch + 1, i + 1, running_loss / 100))


    for k, v in m_GAN.state_dict().iteritems():
        if k.startswith("_") or "Float" not in torch.typename(v): continue
        v.requires_grad = True

    return running_loss/100


datasets=cancer_type_dataset.CANCER_TYPES
torch_dataset=CancerTypesDataset(dataset_names=cancer_type_dataset.CANCER_TYPES, meta_groups_files=cancer_type_dataset.META_GROUPS, metagroups_names=["{}".format(x) for i_x, x in enumerate(cancer_type_dataset.CANCER_TYPES)])
train_dataset,test_dataset = torch.utils.data.random_split(torch_dataset, [torch_dataset.__len__()-torch_dataset.__len__()/100, torch_dataset.__len__()/100])

print "train: {}, test: {}".format(len(train_dataset), len(test_dataset))

trainloader = torch.utils.data.DataLoader(train_dataset, batch_size=batch_size_train,
                                          shuffle=True, num_workers=num_workers, pin_memory=True)

testloader = torch.utils.data.DataLoader(test_dataset, batch_size=batch_size_val,
                                          shuffle=True, num_workers=num_workers, pin_memory=True)


encoder=vae_gan_bn_after_relu_flex_model.Encoder(n_latent_vector=100)
decoder=vae_gan_bn_after_relu_flex_model.Decoder(n_latent_vector=100)
discriminator=vae_gan_bn_after_relu_flex_model.Discriminator(n_latent_vector=100)
vae_old=vae_bn_after_relu_flex_model.Net(n_latent_vector=100)


model_base_folder=constants.OUTPUT_GLOBAL_DIR
PATH_DISCRIMINATOR= os.path.join(model_base_folder,"GAN_DIS_model") # os.path.join(constants.OUTPUT_GLOBAL_DIR, "VAE_model")
PATH_ENCODER= os.path.join(model_base_folder,"GAN_ENC_model")
PATH_DECODER= os.path.join(model_base_folder,"GAN_DEC_model")
load_model=True # False
if load_model and os.path.exists(PATH_ENCODER):
    encoder.load_state_dict(torch.load(PATH_ENCODER))
    encoder.eval()
    decoder.load_state_dict(torch.load(PATH_DECODER))
    decoder.eval()
    discriminator.load_state_dict(torch.load(PATH_DISCRIMINATOR))
    discriminator.eval()


m_VAE = nn.Sequential(encoder,decoder)
# vae_old.load_state_dict(torch.load(PATH_ENCODER))


m_GAN = nn.Sequential(decoder,discriminator)
m_FULL = nn.Sequential(encoder,decoder,discriminator)


# create your optimizer
lr=0.0005
optimizer_dis = optim.Adam(discriminator.parameters(), lr=lr)
optimizer_en = optim.Adam(encoder.parameters(), lr=lr)
optimizer_de = optim.Adam(decoder.parameters(), lr=lr)
optimizer_vae = optim.Adam(m_VAE.parameters(), lr=lr)
optimizer_vae_old = optim.Adam(vae_old.parameters(), lr=lr)
optimizer_gan = optim.Adam(m_GAN.parameters(), lr=lr)
optimizer_full = optim.Adam(m_FULL.parameters(), lr=lr)




n_epoches=1
for meta_epoch in range(0, 1000):  # loop over the dataset multiple times
        print "meta_epoch: {}".format(meta_epoch) 
        discrimination_loss=100
        gen_loss=100

        print "start backprop_vae.."
        backprop_vae(m_VAE, m_FULL, optimizer_full, n_epoches, trainloader, testloader, min(meta_epoch/100.0, 100.0))
        
        # backprop_vae(vae_old, optimizer_vae_old, n_epoches, trainloader, testloader)


        # if meta_epoch > 0:
        # while gen_loss > 0.1:
        print "start backprop_gen.."
        gen_loss=backprop_gen(m_VAE, m_GAN, discriminator, optimizer_gan, n_epoches , trainloader)
               
        # while discrimination_loss > 50: 
        print "start backprop_dis.."
        discrimination_loss, discrimination_accuracy=backprop_dis(m_VAE, discriminator, optimizer_dis, n_epoches, trainloader)
        


        torch.save(encoder.state_dict(), os.path.join(constants.OUTPUT_GLOBAL_DIR, "GAN_ENC_model"))
        torch.save(decoder.state_dict(), os.path.join(constants.OUTPUT_GLOBAL_DIR, "GAN_DEC_model"))
        torch.save(discriminator.state_dict(), os.path.join(constants.OUTPUT_GLOBAL_DIR, "GAN_DIS_model"))
