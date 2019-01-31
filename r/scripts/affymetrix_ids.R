library("biomaRt")
ensembl = useMart(biomart= "ensembl",dataset="hsapiens_gene_ensembl")
affy_ensembl= c("affy_hg_u133_plus_2", "ensembl_gene_id")
getBM(attributes= affy_ensembl, mart= ensembl, values = "*", uniqueRows=T)