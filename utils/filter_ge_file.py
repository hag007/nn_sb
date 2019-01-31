import constants
import infra
import os

# constants.update_dirs(DATASET_NAME="TNFa")
#
# ge_data = infra.load_gene_expression_profile_by_genes(gene_list_file_name="dip_bg.txt")
#
# lns = ["\t".join(cur) for  cur in ge_data]
# file(os.path.join(constants.OUTPUT_DIR, "ge.tsv"), "w+").write("\n".join(lns))

constants.update_dirs(DATASET_NAME_u="TNFa")

ge_data = infra.load_gene_expression_profile_by_genes(gene_expression_path="/media/hag007/Data/omics/GDC-TCGA/UVM/tcga_data/TCGA-UVM.htseq_fpkm.tsv")
lns = ["\t".join([cur[0][:cur[0].index(".")]] +cur[1:]) if i > 0 else "\t".join(cur) for  i, cur in enumerate(ge_data)]
file("/media/hag007/Data/omics/GDC-TCGA/UVM/tcga_data/TCGA-UVM.htseq_fpkm_pz.tsv", "w+").write("\n".join(lns))