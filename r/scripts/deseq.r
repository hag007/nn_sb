library(DESeq2)

coldata <- data.frame(row.names=colnames(data), group)
dds <- DESeqDataSetFromMatrix(countData = data, colData = coldata, design= ~ group)
dds <- DESeq(dds)
DESeq_results <- data.frame(results(dds))
result <- data.frame(DESeq_results[order(DESeq_results$padj),])