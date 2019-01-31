from ensembl2gene_symbol import get_g2e_dictionary
from matplotlib import style
import os
style.use("ggplot")
from infra import *
import logging
import constants
from utils.ensembl2gene_symbol import e2g_convertor
sh = logging.StreamHandler()
logger = logging.getLogger("log")
logger.addHandler(sh)


def load_gene_list(gene_list_file_name, gene_list_path=None):
    if gene_list_path == None:
        gene_list_path = os.path.join(constants.LIST_DIR,gene_list_file_name)
    f = open(gene_list_path,'r')
    lines = [l.strip() for l in f]
    f.close()
    return lines

def load_mirna_clusters(gene_list_file_name, gene_list_path=None):
    if gene_list_path == None:
        gene_list_path = os.path.join(constants.DICTIONARIES_DIR,gene_list_file_name)
    f = open(gene_list_path,'r')
    lines = [l.strip().split() for l in f]
    f.close()
    return lines

def output_mirna_list_file(mrna_id, mirnas, origin_list, gene_list_path=None, source="GDC-TCGA",dataset="melanoma"):
    if gene_list_path == None:
        gene_list_path = os.path.join(constants.LIST_DIR,"mir-{}-{}.txt".format(origin_list, mrna_id))
    f = open(gene_list_path,'w')
    f.write('\n'.join(mirnas))
    f.close()
    return gene_list_path.split("\\")[-1]

def main():
    mirna_clusters = load_mirna_clusters("mir_clusters_by_targets.txt")
    mir_dict = {}
    counter = 0
    while len(mirna_clusters) !=0:
        print counter
        cur_mrna = mirna_clusters[0]
        for cur_mirna in cur_mrna[1:]:
            if mir_dict.has_key(cur_mirna):
                mir_dict[cur_mirna].append(cur_mrna[0])
            else:
                mir_dict[cur_mirna] = [cur_mrna[0]]
        del mirna_clusters[0]
        counter+=1

    with file(os.path.join(constants.DICTIONARIES_DIR, "mir_to_mrna.txt"), "w+") as f:
        for k,v in mir_dict.iteritems():
            f.write("{}\t{}\r\n".format(k, "\t".join(v)))

main()