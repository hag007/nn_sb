from infra import *
import random

def generate_random_set(random_size, meta_gene_set):
    gene_list = load_gene_list(meta_gene_set)
    random.shuffle(gene_list)
    gene_list = gene_list[:random_size]
    output_file_name = "random_list_{}.txt".format(random_size)
    f = open(os.path.join(constants.LIST_DIR,output_file_name),"w+")
    f.writelines("\n".join(gene_list))
    f.close()
    return output_file_name