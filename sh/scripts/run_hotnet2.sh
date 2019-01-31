hotnet2=/home/hag007/repos/bnetworks_alg/hotnet2 # /home/hag007/networks_algo/hotnet2
cache_folder=/media/hag007/Data/bnet/datasets/GE_random_TNFa_2/cache # /home/hag007/bnet/datasets/TNFa_2/cache
output_folder=/media/hag007/Data/bnet/datasets/GE_random_TNFa_2/output # /home/hag007/bnet/datasets/TNFa_2/output
network_name=dip
num_cores=-1
num_network_permutations=1 # 100
num_heat_permutations=1 # 1000

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
