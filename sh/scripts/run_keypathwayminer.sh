base_folder=/media/hag007/Data/bnet
network_folder=$base_folder/networks
datasets_folder=$base_folder/datasets
data_folder=$dataset_folder/data
cache_folder=$dataset_folder/cache
output_folder=$dataset_folder/output
dataset_folder=/media/hag007/Data/bnet/datasets/GE_random_TNFa_2

pwd

java -jar -Xmx2G KPM-4.0.jar -strategy=INES -algo=GREEDY -K=1  -L1=420
