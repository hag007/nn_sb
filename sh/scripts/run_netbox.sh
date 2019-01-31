# base_folder=/media/hag007/Data/bnet
# network_folder=$base_folder/networks
# datasets_folder=$base_folder/datasets
# data_folder=$dataset_folder/data
# cache_folder=$dataset_folder/cache
# output_folder=$dataset_folder/output
# dataset_folder=/media/hag007/Data/bnet/datasets/GE_random_TNFa_2
export NETBOX_HOME=/home/hag007/repos/bnetworks_alg/netbox # ../repos/netbox
export PATH=$PATH:$NETBOX_HOME/bin


python $NETBOX_HOME/bin/netAnalyze.py ${CONFIG_FILE_NAME}

