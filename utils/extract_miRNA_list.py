from ensembl2gene_symbol import get_g2e_dictionary
from matplotlib import style
import os
style.use("ggplot")
import logging
sh = logging.StreamHandler()
import constants
logger = logging.getLogger("log")
logger.addHandler(sh)

LABEL_ID = "sample_type.samples"
PRIMARY_TUMOR = "Primary Tumor"
METASTATIC = "Metastatic"

def load_gene_list(gene_list_file_name, gene_list_path=None, source="GDC-TCGA",dataset="melanoma"): #  ="TCGA-SKCM.htseq_counts.tsv"
    if gene_list_path == None:
        gene_list_path = os.path.join(constants.BASE_PROFILE,source,dataset,gene_list_file_name)
    f = open(gene_list_path,'r')
    lines = [l.strip() for l in f]
    f.close()
    return lines

def load_mirna_clusters(gene_list_file_name, gene_list_path=None, source="GDC-TCGA",dataset="melanoma"):
    if gene_list_path == None:
        gene_list_path = os.path.join(constants.BASE_PROFILE,source,dataset,gene_list_file_name)
    f = open(gene_list_path,'r')
    lines = [l.strip() for l in f]
    f.close()
    return lines

def output_mirna_list_file(mrna_id, mirnas, origin_list, gene_list_path=None, source="GDC-TCGA",dataset="melanoma"):
    if gene_list_path == None:
        gene_list_path = os.path.join(BASE_PROFILE,source,dataset,"list","mir-{}-{}.txt".format(origin_list, mrna_id))
    f = open(gene_list_path,'w')
    f.write('\n'.join(mirnas))
    f.close()
    return gene_list_path.split("\\")[-1]


output_files = []
origin_lists = ["tami_metabolism.txt", "tami_mito.txt", "tami_extra.txt"]
for origin_list in origin_lists:
    print "fetch mrna keys from {}".format(origin_list)
    mRNA_list = load_gene_list(origin_list)
    gene_symbols2ensembl = get_g2e_dictionary()
    mirna_clusters = load_mirna_clusters("mir_cluster_out.txt")
    ensembl_ids = []

    counter=0
    for cur in mRNA_list:
        if gene_symbols2ensembl.has_key(cur):
            ensembl_ids.append(gene_symbols2ensembl[cur].split('.')[0])
            print "{} retrieved from mirna cluster file".format(cur)
            counter+=1
        else:
            print "{} not found in mirna cluster file".format(cur)
    print "total mrna found in dictionary: {}/{}".format(counter, len(mRNA_list))

    counter=0

    for cur_line in mirna_clusters:
        cur_gene = cur_line.split()[0]
        if cur_gene in ensembl_ids:
            print "{} is found in mirna cluster file".format(cur_gene)
            counter+=1
            output_files.append(output_mirna_list_file(mRNA_list[ensembl_ids.index(cur_gene)], cur_line.split()[1:], origin_list.split('.')[0]))
            #ensembl_ids.remove(cur_gene)
    print "total mrna found in mirna cluser: {}/{}".format(counter, len(ensembl_ids))

print "output files: \n{}".format(output_files)



