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
import torch_vae_gan_copy_model

batch_size_train=100
batch_size_val=10
num_workers=50

criterion = nn.BCELoss(reduction='sum')

datasets=torch_dataset_cancer.CANCER_TYPES
torch_dataset=torch_dataset_cancer.CancerTypesDataset(dataset_names=torch_dataset_cancer.CANCER_TYPES, meta_groups_files=torch_dataset_cancer.META_GROUPS, metagroups_names=["{}_{}".format(x, i_x) for i_x, x in enumerate(torch_dataset_cancer.CANCER_TYPES)])
train_dataset,test_dataset = torch.utils.data.random_split(torch_dataset, [torch_dataset.__len__()-torch_dataset.__len__()/10, torch_dataset.__len__()/10])

print "n_train samples: {}, n_test samples: {}".format(len(train_dataset), len(test_dataset))

trainloader = torch.utils.data.DataLoader(train_dataset, batch_size=batch_size_train, shuffle=True, num_workers=num_workers, pin_memory=True)

testloader = torch.utils.data.DataLoader(test_dataset, batch_size=batch_size_val, shuffle=True, num_workers=num_workers, pin_memory=True)


# define constant 
max_epochs = 100000
lr = 3e-4
beta = 5
alpha = 1
gamma = 5
delta = 1 
n_latent_vector=2
G = torch_vae_gan_copy_model.VAE_GAN_Generator(n_latent_vector=n_latent_vector)
D = torch_vae_gan_copy_model.Discriminator(n_latent_vector=n_latent_vector)

opt_enc = optim.Adam(G.encoder.parameters(), lr=lr)
opt_dec = optim.Adam(G.decoder.parameters(), lr=lr)
opt_dis = optim.Adam(D.parameters(), lr=lr * alpha)

fixed_noise = Variable(torch.randn(batch_size_train, n_latent_vector))
data, _ = next(iter(trainloader))
fixed_batch = Variable(data)
min_val=1000000
min_val_epoch=-1

for epoch in range(max_epochs):
    train_loss=0
    val_loss=0
    print "cur epoch: {}".format(epoch)
    D_real_list, D_rec_enc_list, D_rec_noise_list, D_list = [], [], [], []
    g_loss_list, rec_loss_list, prior_loss_list = [], [], []
    for i, data_tuple in enumerate(trainloader,0):
        data, _ = data_tuple
        batch_size = data.size()[0]
        ones_label = Variable(torch.ones(batch_size))
        zeros_label = Variable(torch.zeros(batch_size))

        datav = Variable(data)
        rec_enc = G(datav)
        mean, logvar = G.encoder.encode(datav)
        noisev = Variable(torch.randn(batch_size, n_latent_vector))
        rec_noise = G.decoder(noisev)

        # train discriminator
        output, l = D(datav)
        errD_real = criterion(output, ones_label)
        D_real_list.append(output.data.mean())
        output, l = D(rec_enc)
        errD_rec_enc = criterion(output, zeros_label)
        D_rec_enc_list.append(output.data.mean())
        output, l = D(rec_noise)
        errD_rec_noise = criterion(output, zeros_label)
        D_rec_noise_list.append(output.data.mean())

        dis_img_loss = errD_real + errD_rec_enc + errD_rec_noise
        D_list.append(dis_img_loss.data.mean())
        opt_dis.zero_grad()
        dis_img_loss.backward(retain_graph=True)
        opt_dis.step()

        # train decoder
        output, l_real = D(datav)
        errD_real = criterion(output, ones_label)
        output, l_rec = D(rec_enc)
        errD_rec_enc = criterion(output, zeros_label)
        output, l = D(rec_noise)
        errD_rec_noise = criterion(output, zeros_label)

        similarity_rec_enc = rec_enc
        similarity_data = datav

        dis_img_loss = errD_real + errD_rec_enc + errD_rec_noise
        gen_img_loss = - dis_img_loss

        g_loss_list.append(gen_img_loss.data.mean())
        rec_loss = criterion(similarity_rec_enc, similarity_data.detach())
        rec_loss_list.append(rec_loss.data.mean())
        err_dec = gamma * rec_loss + gen_img_loss
   
        opt_dec.zero_grad()
        err_dec.backward(retain_graph=True)
        opt_dec.step()
        
        # train encoder
        prior_loss = 1 + logvar - mean.pow(2) - logvar.exp()
        prior_loss = (-0.5 * torch.sum(prior_loss)) # / torch.numel(mean.data)
        # print (prior_loss, mean, std)
        prior_loss_list.append(prior_loss.data.mean())
        err_enc = prior_loss + beta * rec_loss

        train_loss+=prior_loss.item() + rec_loss.item()
        
        if i % 10 == 9:
            print('[%d, %5d] train loss: %.3f' %
                  (epoch + 1, i + 1, train_loss / 100)) 
            train_loss = 0.0

        opt_enc.zero_grad()
        err_enc.backward()
        opt_enc.step()

    for i, data in enumerate(testloader, 0):
        with torch.no_grad():
            inputs, labels = data

            rec_enc = G(datav)    
            mean, logvar = G.encoder.encode(datav)
    
            prior_loss = 1 + logvar - mean.pow(2) - logvar.exp()
            prior_loss = (-0.5 * torch.sum(prior_loss)) #  
        

            similarity_rec_enc = rec_enc
            similarity_data = datav

            
            rec_loss = criterion(similarity_rec_enc, similarity_data.detach())
            loss = prior_loss + rec_loss
            val_loss += loss.item()

    
    print('[%d, %5d] val loss: %.3f' %
          (epoch + 1, i + 1, val_loss / 100))

    samples = G.decoder(fixed_noise)

    localtime = time.asctime(time.localtime(time.time()))
    print (localtime)
    print ('[%d/%d]: D_real:%.4f, D_enc:%.4f, D_noise:%.4f, Loss_D:%.4f, Loss_G:%.4f, rec_loss:%.4f, prior_loss:%.4f'
           % (epoch,
              max_epochs,
              np.mean(D_real_list),
              np.mean(D_rec_enc_list),
              np.mean(D_rec_noise_list),
              np.mean(D_list),
              np.mean(g_loss_list),
              np.mean(rec_loss_list),
              np.mean(prior_loss_list)))

    model_base_folder=constants.OUTPUT_GLOBAL_DIR
    PATH_DISCRIMINATOR= os.path.join(model_base_folder,"GAN_DIS_mdl")
    PATH_ENCODER= os.path.join(model_base_folder,"GAN_ENC_mdl")
    PATH_DECODER= os.path.join(model_base_folder,"GAN_DEC_mdl")
    if min_val > val_loss/100:
        min_val=val_loss/100
        min_val_epoch=epoch
        torch.save(G.encoder.state_dict(), os.path.join(constants.OUTPUT_GLOBAL_DIR, "GAN_ENC_mdl"))
        torch.save(G.decoder.state_dict(), os.path.join(constants.OUTPUT_GLOBAL_DIR, "GAN_DEC_mdl"))
        torch.save(D.state_dict(), os.path.join(constants.OUTPUT_GLOBAL_DIR, "GAN_DIS_mdl"))

    torch.save(G.encoder.state_dict(), os.path.join(constants.OUTPUT_GLOBAL_DIR, "GAN_ENC_mdl_"))
    torch.save(G.decoder.state_dict(), os.path.join(constants.OUTPUT_GLOBAL_DIR, "GAN_DEC_mdl_"))
    torch.save(D.state_dict(), os.path.join(constants.OUTPUT_GLOBAL_DIR, "GAN_DIS_mdl_"))

    print "min_val epoch: {}, min_val: {}".format(min_val_epoch, min_val)

