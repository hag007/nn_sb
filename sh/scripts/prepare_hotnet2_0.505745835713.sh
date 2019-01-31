hotnet2=/specific/netapp5/gaga/hagailevi/evaluation/bnetworks_alg/hotnet2 # /home/hag007/networks_algo/hotnet2
cache_folder=/specific/netapp5/gaga/hagailevi/evaluation/bnet/datasets/GE_random_ROR_2_hotnet2_1388/cache # /home/hag007/bnet/datasets/TNFa_2/cache
num_cores=5
num_network_permutations=10
num_heat_permutations=100

source $hotnet2/venv/bin/activate

# Create network data.
$hotnet2/venv/bin/python $hotnet2/makeNetworkFiles.py \
    -e  $cache_folder/hotnet2_edges.txt \
    -i  $cache_folder/hotnet2_vertices.txt \
    -nn dip \
    -p  dip \
    -b  0.5 \
    -q 3 \
    -o  $cache_folder/dip \
    -np $num_network_permutations \
    -c  $num_cores

$hotnet2/venv/bin/python $hotnet2/makeHeatFile.py \
    scores \
    -hf $cache_folder/heatfile.txt \
    -o  $cache_folder/heatfile.json \
    -n  heatfile
