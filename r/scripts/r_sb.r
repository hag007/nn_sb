## try http:// if https:// URLs are not supported

# source("https://bioconductor.org/biocLite.R")
# biocLite("edgeR")
# biocLite("locfit")
# biocLite("STRINGdb")

# biocLite("DLBCL")
# biocLite('biomaRt')
# biocLite("Rgraphviz")

# library("edgeR")
# library("DESeq2")
# library('biomaRt')
# # library("STRINGdb")
# library("Rgraphviz")
# library("BioNet")



library("Rgraphviz")
library("BioNet")

library(edgeR)
  
group <- read.delim("/home/hag007/bnet/datasets/BRCA/data/classes.tsv", header=FALSE)
data <- read.delim("/home/hag007/bnet/datasets/BRCA/data/ge.tsv", row.names = 1)
genes <- rownames(data)
y <- DGEList(counts=data,group=group)
y <- calcNormFactors(y)
design <- model.matrix(~group)
y <- estimateDisp(y,design)
if (is.na(y$common.dispersion)){
  y$common.dispersion=0.05
}
et <- exactTest(y)
edgeR_results = topTags(et, n = nrow(data))$table
# rownames(edgeR_results) <- genes
result <- data.frame(edgeR_results[order(edgeR_results$FDR),])




##load dictionary
# p.dict <- read.delim("/home/hag007/bnet/dictionaries/esnp_dictionary.txt")
# p.dict <- p.dict[!(duplicated(p.dict[["ENSP"]])),]
# rownames(p.dict) <- p.dict[["ENSP"]]
# p.dict <- p.dict["ENSG"]

##load STRING ppi
# ppi <- read.delim("/home/hag007/bnet/dictionaries/string_ppi_small.txt")
# ppi <- ppi[c(1,2)] 
# ppi <- data.frame(lapply(ppi, function(x) {
#   }))
# 
# ppi <- data.frame(lapply(ppi, function(x) {
#   x <-  p.dict[x,1]
# }))
# # flatten.df <- c(t(ppi))
# ig <- graph_from_edgelist((as.matrix(ppi)), directed = FALSE)

# pvals <- cbind(t = dataLym$t.pval, s = dataLym$s.pval)
# rownames(pvals) <- dataLym$label
# pval <- aggrPvals(pvals, order = 2, plot = FALSE)
# subnet <- subNetwork(dataLym$label, interactome)
# subnet <- rmSelfLoops(subnet)
# fb <- fitBumModel(pval, plot = FALSE)
# scores <- scoreNodes(subnet, fb, fdr = 0.001)
# module <- runFastHeinz(subnet, scores)
# logFC <- dataLym$diff
# names(logFC) <- dataLym$label
# plotModule(module, scores = scores, diff.expr = logFC)


# library('biomaRt')
# mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
# genes <- rownames(data)
# G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)

# e2g = data.frame(data=G_list[2])
# rownames(e2g) <-G_list[[1]]
# combined.dataset <- merge(data,e2g, by=0)
# combined.dataset <- combined.dataset[!duplicated(combined.dataset['hgnc_symbol']), ]
# ensemblids <- combined.dataset["Row.names"]
# hgnc.symbols <- combined.dataset[["hgnc_symbol"]]
# rownames(combined.dataset) <- hgnc.symbols
# combined.dataset['ensembl_id'] = ensemblids
# combined.dataset = combined.dataset[ , !(names(combined.dataset) %in% c("Row.names", "hgnc_symbol"))]


##load DIP ppi
# ig <- loadNetwork.sif("/home/hag007/bnet/output/dip_out.sif")
# data <- read.delim("/home/hag007/bnet/output/DEG_deseq.tsv", sep="\t", stringsAsFactors=FALSE)
# rownames(data)  <- data[["ensembl_id"]]
# data  <- data[names(data) !="ensembl_id"]
# pval <- data[["deseq_p"]]
# names(pval) <- rownames(data)
# combined.dataset <- data
# subnet <- subNetwork(rownames(combined.dataset), ig)
# subnet <- rmSelfLoops(subnet)
# 
# pval <-na.omit(pval)
# # pval <- pval[!pval==1]
# fb <- fitBumModel(pval, plot = TRUE, starts = 10)
# scores <- scoreNodes(subnet, fb, fdr = 1) #0.9999999999 
# scores = -log10(pval)
# module <- runFastHeinz(subnet, scores)
# # logFC <- dataLym$diff
# # names(logFC) <- dataLym$label
# plotModule(module, scores = scores) # , diff.expr = logFC
# write(nodes(module), file = "data.txt", sep="\n")


# x <- read.delim("/home/hag007/bnet/datasets/TNFa_2/ge.tsv",row.names="id")
# y <- DGEList(counts=x,group=group)
# y <- calcNormFactors(y)
# design <- model.matrix(~group)
# y <- estimateDisp(y,design)
# #plotBCV(y)
# et <- exactTest(y)
# results_edgeR <- topTags(et, n = nrow(x))
#
# coldata <- data.frame(row.names=colnames(x), group)
# dds <- DESeqDataSetFromMatrix(countData = x, colData = coldata, design= ~ group)
# dds <- DESeq(dds)
# resultsNames(dds) # lists the coefficients
# res <- results(dds)
# bb <- res$padj < .05
# bb


# b <- sum(results_edgeR$table$FDR < .05)
## plotSmear(et, de.tags = rownames(results_edgeR)[results_edgeR$table$FDR < .05])
## abline(h = c(-2, 2), col = "blue")