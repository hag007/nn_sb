.libPaths("/specific/netapp5/gaga/hagailevi/evaluation/Renv")
library(edgeR, lib.loc  = "/specific/netapp5/gaga/hagailevi/evaluation/Renv")
group=unlist(group)
print(dim(group))
print(dim(data))
y <- DGEList(counts=data,group=group)
y <- calcNormFactors(y)
design <- model.matrix(~group)
y <- estimateDisp(y,design)
if (is.na(y$common.dispersion)){
  y$common.dispersion=0.1
}
et <- exactTest(y)
edgeR_results = topTags(et, n = nrow(data))$table
# rownames(edgeR_results) <- genes
result <- data.frame(edgeR_results[order(edgeR_results$FDR),])