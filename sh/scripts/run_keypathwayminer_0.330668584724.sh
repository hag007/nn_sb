base_folder=/specific/netapp5/gaga/hagailevi/evaluation/bnet
network_folder=$base_folder/networks
datasets_folder=$base_folder/datasets
data_folder=$dataset_folder/data
cache_folder=$dataset_folder/cache
output_folder=$dataset_folder/output
dataset_folder=/specific/netapp5/gaga/hagailevi/evaluation/bnet/datasets/GE_random_SHEZH_1_keypathwayminer_INES_GREEDY_168

pwd

java -jar -Xmx2G KPM-4.0.jar -strategy=INES -algo=GREEDY -K=1  -L1=420
