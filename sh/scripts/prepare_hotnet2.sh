hotnet2=/home/hag007/repos/bnetworks_alg/hotnet2 # /home/hag007/networks_algo/hotnet2
cache_folder=/media/hag007/Data/bnet/datasets/GE_random_TNFa_2/cache # /home/hag007/bnet/datasets/TNFa_2/cache
num_cores=-1
num_network_permutations=1 # 100
num_heat_permutations=1 # 1000

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