base_folder=/specific/netapp5/gaga/hagailevi/evaluation/bnet
networks_folder=$base_folder/networks
datasets_folder=$base_folder/datasets
dataset_folder=/specific/netapp5/gaga/hagailevi/evaluation/bnet/datasets/GE_random_HC12_jactivemodules_sa_2159
data_folder=$dataset_folder/data
cache_folder=$dataset_folder/cache
output_folder=$dataset_folder/output
is_greedy=False
algo_dir=/specific/netapp5/gaga/hagailevi/evaluation/bnetworks_alg/jactivemodules
num_of_modules=10
overlap_threshold=0

echo /specific/netapp5/gaga/hagailevi/evaluation/bnet/datasets/GE_random_HC12_jactivemodules_sa_2159/cache/deg_edger.tsv


java -jar $algo_dir/jactivemodules.jar \
                                 /specific/netapp5/gaga/hagailevi/evaluation/bnet/networks/dip.sif \
                                 /specific/netapp5/gaga/hagailevi/evaluation/bnet/datasets/GE_random_HC12_jactivemodules_sa_2159/cache/deg_edger.tsv \
                                 $is_greedy \
                                 $num_of_modules \
                                 $overlap_threshold \
                                 /specific/netapp5/gaga/hagailevi/evaluation/bnet/datasets/GE_random_HC12_jactivemodules_sa_2159/output/jactivemodules_sa_results.txt

