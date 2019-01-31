from infra import load_phenotype_data
from infra import load_gene_list
import numpy as np
import os

SEPARATOR = "@%@"

LABEL_ID = "sample_type.samples"
PRIMARY_TUMOR = "Primary Tumor"
METASTATIC = "Metastatic"

# return gene expression table filtered according an external list with proper orientation (0 degree angle according genes, 90 degree angle according patients)
def filter_gene_expression_profile(gene_list_file_name, gene_expression_file_name, gene_filter_file_name=None, gene_list_path=None, gene_expression_path=None, gene_filter_path=None, source="GDC-TCGA",dataset="melanoma",by_gene=False):
    pd = load_phenotype_data("TCGA-SKCM.GDC_phenotype.tsv", phenotype_list_path=None, source="GDC-TCGA",
                             dataset="melanoma")
    pd_headers = pd[0]
    label_index = [i for i, v in enumerate(pd_headers) if v == LABEL_ID][0]
    invalid_labels = [cur[0] for i, cur in enumerate(pd) if
                      cur[label_index] == METASTATIC or cur[label_index] == PRIMARY_TUMOR or i == 0]

    gene_list = load_gene_list(gene_list_file_name=gene_list_file_name, gene_list_path=gene_list_path, source=source,dataset=dataset)
    if gene_filter_file_name:
        filter_gene_list = load_gene_list(gene_list_file_name=gene_filter_file_name, gene_list_path=gene_filter_path,
                                          source=source, dataset=dataset)
        gene_list = [cur for cur in gene_list if cur in filter_gene_list or cur[:cur.find('.')] in filter_gene_list]

    if gene_expression_path == None:
        gene_expression_path = os.path.join(BASE_PROFILE, source, dataset, gene_expression_file_name)
    f = open(gene_expression_path,'r')
    expression_profiles_filtered = [l for i, l in enumerate(f) if i==0 or any([l.strip()[0:l.strip().find('\t')] in gene_list or l.strip()[0:l.strip().find('\t')].split(".")[0] in gene_list])]
    f.close()
    expression_profiles_filtered_out = []
    expression_headers = expression_profiles_filtered[0].split()
    for cur in expression_profiles_filtered:
        splited = cur.split("\t")
        splited = [cur for i, cur in enumerate(splited) if expression_headers[i] not in invalid_labels ]
        expression_profiles_filtered_out.append("\t".join(splited))

    f = open(gene_expression_path+"_filtered", 'w+')
    f.writelines(expression_profiles_filtered)
    f.close()

filter_gene_expression_profile("protein_coding.txt", "TCGA-SKCM.htseq_counts.tsv", gene_filter_file_name=None, gene_list_path=None, gene_expression_path=None, gene_filter_path=None, source="GDC-TCGA",dataset="melanoma", by_gene = True)
