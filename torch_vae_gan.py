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
import argparse
from matplotlib.lines import Line2D



num_workers=25
batch_size_train=100
batch_size_val=10
REDUCTION='sum'
GAMMA=15
ALPHA=0.1
BETA=5

def get_encoder_loss(m_VAEGAN, inputs, KLD_ratio):
    z, mu, var = m_VAEGAN.vae.encoder(inputs)
    ll_loss = get_ll_dis_loss(m_VAEGAN, inputs)
    prior_loss = get_prior_loss(mu, var, KLD_ratio)

    return ll_loss + prior_loss


def get_ll_dis_loss(m_VAEGAN, inputs):
	
    auth_decoded, l_decoded, _1, _2, _3, _4, = m_VAEGAN(inputs)

    # with torch.no_grad:
    auth_real, l_real, _1, _2, _3, _4 = m_VAEGAN(inputs)

    BCE = F.binary_cross_entropy(l_decoded, l_real.detach(),  reduction='sum')
    return BCE

def get_rec_loss(m_VAEGAN, inputs):

    decoded, z, mu, logvar = m_VAEGAN.vae(inputs)
    BCE = F.binary_cross_entropy(decoded, inputs, reduction='sum')
    # BCE = ((decoded - inputs) ** 2).sum() 
    return BCE


def get_prior_loss(mu, logvar, KLD_ratio):

    KLD = -0.5 * torch.sum(1 + logvar - mu.pow(2) - logvar.exp())

    return KLD # *KLD_ratio


def backprop_encoder(m_VAEGAN, optimizers, n_epoches, trainloader, testloader, KLD_ratio):

     for epoch in range(n_epoches):
        train_loss=0
        val_loss=0
        for i, data in enumerate(trainloader, 0):

            inputs, labels = data
            [optimizer.zero_grad() for optimizer in optimizers]
     
            z, mu, var = m_VAEGAN.vae.encoder(inputs)
            ll_loss = get_ll_dis_loss(m_VAEGAN, inputs)
            prior_loss = get_prior_loss(mu, var, KLD_ratio)
            loss=BETA*ll_loss + prior_loss
          
            loss.backward()
            train_loss += loss.item()
            [optimizer.step() for optimizer in optimizers]

            if i % 10 == 9:
                print('[%d, %5d] vae train loss: %.3f' %
                      (epoch + 1, i + 1, train_loss / 100))
                train_loss = 0.0

        for i, data in enumerate(testloader, 0):
            with torch.no_grad():
                inputs, labels = data
                loss = get_encoder_loss(m_VAEGAN, inputs, KLD_ratio)
                val_loss += loss.item()


        print('[%d, %5d] vae val loss: %.3f' %
              (epoch + 1, i + 1, val_loss / 100))



def backprop_dis(m_VAEGAN, optimizer, n_epoches, trainloader):

    accuracy=0
    for epoch in range(n_epoches):
        batch_accuracies=[]
        running_loss = 0.0

        for i, data in enumerate(trainloader, 0):

            inputs, labels = data

            true_labels = torch.ones(len(data[0])) * 0.05
            fake_labels = torch.ones(len(data[0])) * 0.95


            real_inputs = inputs
            real_labels = true_labels

            decoded_inputs, _1, _2, _3 = m_VAEGAN.vae(inputs)
            decoded_labels = fake_labels

            random_inputs = m_VAEGAN.vae.decoder.decode(torch.randn(len(labels), m_VAEGAN.vae.n_latent_vector))
            random_labels = fake_labels


            inputs = torch.cat((real_inputs, decoded_inputs,random_inputs), dim=0)
            labels = torch.cat((real_labels, decoded_labels, random_labels), dim=0)

            auth, l = m_VAEGAN.discriminator.discriminate(inputs)


            optimizer.zero_grad()
            loss = torch.nn.BCELoss(reduction=REDUCTION)(auth.view(-1, ), labels)
            loss.backward()
            running_loss += loss.item()
            optimizer.step()

            correct = 0
            total = 0
            y = []
            with torch.no_grad():
                y_pred_raw, l = m_VAEGAN.discriminator.discriminate(inputs)
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


def backprop_gen(m_VAEGAN, optimizer, n_epoches, trainloader):

    for epoch in range(n_epoches):
        running_loss = 0.0
        print "cur epoch: {}".format(epoch)
        for i, data in enumerate(trainloader, 0):

            # get the inputs
            inputs, labels = data


            true_labels = torch.ones(len(data[0])) * 0.95
            fake_labels = torch.ones(len(data[0])) * 0.05

            real_inputs = inputs
            real_labels = true_labels

            decoded_inputs, z,  mu, var = m_VAEGAN.vae(inputs)
            decoded_labels = fake_labels

            random_inputs = m_VAEGAN.vae.decoder.decode(torch.randn(len(labels), m_VAEGAN.n_latent_vector))
            random_labels = fake_labels

            inputs = torch.cat((real_inputs, decoded_inputs,random_inputs), dim=0)
            labels = torch.cat((real_labels, decoded_labels, random_labels), dim=0)

            auth, l = m_VAEGAN.discriminator.discriminate(inputs)

            optimizer.zero_grad()

            loss = torch.nn.BCELoss(reduction=REDUCTION)(auth.view(-1, ), labels) + GAMMA*get_ll_dis_loss(m_VAEGAN, inputs)
          
            loss.backward()
            running_loss += loss.item()
            optimizer.step()

        # print statistics
        print('[%d, %5d] gen train loss: %.3f' %
              (epoch + 1, i + 1, running_loss / 100))


    return running_loss/100



def main():
    parser = argparse.ArgumentParser(description='args')
    parser.add_argument('--n_latent_vector', dest='n_latent_vector', default='100')
    parser.add_argument('--load_model', dest='load_model', default="false")

    args = parser.parse_args()

    load_model=args.load_model=='true'
    n_latent_vector=int(args.n_latent_vector)

    torch_dataset=CancerTypesDataset(dataset_names=cancer_type_dataset.CANCER_TYPES, meta_groups_files=cancer_type_dataset.META_GROUPS, metagroups_names=["{}".format(x) for i_x, x in enumerate(cancer_type_dataset.CANCER_TYPES)])
    train_dataset,test_dataset = torch.utils.data.random_split(torch_dataset, [torch_dataset.__len__()-torch_dataset.__len__()/100, torch_dataset.__len__()/100])

    print "n_train samples: {}, n_test samples: {}".format(len(train_dataset), len(test_dataset))

    trainloader = torch.utils.data.DataLoader(train_dataset, batch_size=batch_size_train,
                                              shuffle=True, num_workers=num_workers, pin_memory=True)

    testloader = torch.utils.data.DataLoader(test_dataset, batch_size=batch_size_val,
                                              shuffle=True, num_workers=num_workers, pin_memory=True)

    m_encoder=vae_gan_bn_after_relu_flex_model.Encoder(n_latent_vector=n_latent_vector)
    m_decoder=vae_gan_bn_after_relu_flex_model.Decoder(n_latent_vector=n_latent_vector)
    m_discriminator=vae_gan_bn_after_relu_flex_model.Discriminator(n_latent_vector=n_latent_vector)


    model_base_folder=constants.OUTPUT_GLOBAL_DIR
    PATH_DISCRIMINATOR= os.path.join(model_base_folder,"GAN_DIS_model") # os.path.join(constants.OUTPUT_GLOBAL_DIR, "VAE_model")
    PATH_ENCODER= os.path.join(model_base_folder,"GAN_ENC_model")
    PATH_DECODER= os.path.join(model_base_folder,"GAN_DEC_model")

    if load_model and os.path.exists(PATH_ENCODER):
        m_encoder.load_state_dict(torch.load(PATH_ENCODER))
        m_encoder.eval()
        m_decoder.load_state_dict(torch.load(PATH_DECODER))
        m_decoder.eval()
        m_discriminator.load_state_dict(torch.load(PATH_DISCRIMINATOR))
        m_discriminator.eval()


    m_VAE = vae_gan_bn_after_relu_flex_model.VAE(m_encoder,m_decoder)
    m_GAN = vae_gan_bn_after_relu_flex_model.GAN(m_decoder,m_discriminator)
    m_VAEGAN=vae_gan_bn_after_relu_flex_model.VAEGAN(m_VAE,m_discriminator)


    # create your optimizer
    lr=0.001
    o_discriminator = optim.Adam(m_discriminator.parameters(), lr=lr*ALPHA)
    o_encoder = optim.Adam(m_encoder.parameters(), lr=lr)
    o_decoder = optim.Adam(m_decoder.parameters(), lr=lr)
    # o_vae = optim.Adam(m_VAE.parameters(), lr=lr)
    # o_gan = optim.Adam(m_GAN.parameters(), lr=lr)
    # o_vaegan = optim.Adam(m_VAEGAN.parameters(), lr=lr)




    n_epoches=1
    for meta_epoch in range(0, 1000):  # loop over the dataset multiple times
            print "meta_epoch: {}".format(meta_epoch)
            discrimination_loss=100
            gen_loss=100

            print "start backprop_vae.."
            backprop_encoder(m_VAEGAN, [o_encoder], n_epoches, trainloader, testloader, min(meta_epoch / 100.0, 1.0))

            # while gen_loss > 0.1:
            print "start backprop_gen.."
            gen_loss=backprop_gen(m_VAEGAN, o_decoder, n_epoches, trainloader)

            # while discrimination_loss > 50:
            print "start backprop_dis.."
            discrimination_loss, discrimination_accuracy=backprop_dis(m_VAEGAN, o_discriminator, n_epoches, trainloader)


            torch.save(m_encoder.state_dict(), os.path.join(constants.OUTPUT_GLOBAL_DIR, "GAN_ENC_model"))
            torch.save(m_decoder.state_dict(), os.path.join(constants.OUTPUT_GLOBAL_DIR, "GAN_DEC_model"))
            torch.save(m_discriminator.state_dict(), os.path.join(constants.OUTPUT_GLOBAL_DIR, "GAN_DIS_model"))


if __name__ == "__main__":
    main()
