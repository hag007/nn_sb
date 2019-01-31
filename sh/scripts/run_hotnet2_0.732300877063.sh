hotnet2=/specific/netapp5/gaga/hagailevi/evaluation/bnetworks_alg/hotnet2 # /home/hag007/networks_algo/hotnet2
cache_folder=/specific/netapp5/gaga/hagailevi/evaluation/bnet/datasets/GE_random_TNFa_2_hotnet2_1551/cache # /home/hag007/bnet/datasets/TNFa_2/cache
output_folder=/specific/netapp5/gaga/hagailevi/evaluation/bnet/datasets/GE_random_TNFa_2_hotnet2_1551/output # /home/hag007/bnet/datasets/TNFa_2/output
network_name=dip
num_cores=5
num_network_permutations=10
num_heat_permutations=100

source $hotnet2/venv/bin/activate

# Run HotNet2.
$hotnet2/venv/bin/python $hotnet2/HotNet2.py \
    -nf  $cache_folder/${network_name}/${network_name}_ppr_0.5.h5 \
    -pnp $cache_folder/${network_name}/permuted/${network_name}_ppr_0.5_1.h5 \
    -hf  $cache_folder/heatfile.json \
    -np  $num_network_permutations \
    -hp  $num_heat_permutations \
    -o   $output_folder/results \
    -c   $num_cores
