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

def main(mrna_list_file_names, mir_list_file_names):
    output_files = []
    mirna_clusters = load_mirna_clusters("mir_clusters_by_targets.txt")
    associated_mirna = []
    for cur_mrna_list in mrna_list_file_names:
        mrna_list = load_gene_list(cur_mrna_list)
        for cur_mir_list in mir_list_file_names:
            mir_list = load_gene_list(cur_mir_list)
            for cur in mirna_clusters:
                if cur[0].split(".")[0] in mrna_list  and len(set(cur[1:]).intersection(mir_list)) != 0:
                    associated_mirna = associated_mirna + list(set(cur[1:]).intersection(mir_list))

    associated_mirna = list(set(associated_mirna))
    associated_mirna = e2g_convertor(associated_mirna)
    f = file(os.path.join(constants.LIST_DIR, "mir_{}.txt".format("_".join([x.split(".")[0] for x in mrna_list_file_names]))), "w+")
    f.write("\r\n".join(associated_mirna))
    f.close()
    print associated_mirna

    return associated_mirna

main()