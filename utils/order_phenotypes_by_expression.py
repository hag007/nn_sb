from infra import load_phenotype_data
from infra import load_gene_expression_profile_by_genes

SEPARATOR = "@%@"

LABEL_ID = "sample_type.samples"
PRIMARY_TUMOR = "Primary Tumor"
METASTATIC = "Metastatic"
pd = load_phenotype_data("TCGA-SKCM.GDC_phenotype.tsv", phenotype_list_path=None, source="GDC-TCGA",dataset="melanoma")
# pd = np.flip(np.rot90(pd, k=1, axes=(1, 0)), 1)
epd = load_gene_expression_profile_by_genes("protein_coding.txt", "TCGA-SKCM.htseq_counts.tsv", gene_filter_file_name=None, gene_list_path=None, gene_expression_path=None, gene_filter_path=None, source="GDC-TCGA",dataset="melanoma")
pd_headers = pd[0]
label_index = [i for i, v in enumerate(pd_headers) if v == LABEL_ID][0]
pd = sorted(filter(lambda i: i[0] in epd[0], pd),key = lambda i: epd[0].index(i[0]))

# for i, cur in enumerate(pd):
#     print pd[i]
#     print epd[0][i]

map = {'Primary Tumor' : '0',
       'Metastatic' : '1',
       'Additional Metastatic' : '3',
       'Solid Tissue Normal' : '2'}
for cur in pd:
    print map[cur[label_index]] + "\t",