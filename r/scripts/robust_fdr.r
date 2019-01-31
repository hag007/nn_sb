# source("https://bioconductor.org/biocLite.R")
# biocLite("prot2D")
.libPaths("/specific/netapp5/gaga/hagailevi/evaluation/Renv")
library(prot2D, lib.loc  = "/specific/netapp5/gaga/hagailevi/evaluation/Renv")

# emp_file_name="/home/hag007/Desktop/aggregate_report/oob/emp_diff_ERS_1_bionet_passed_oob.tsv"
# discrete_interval=0.001

data<-read.delim(emp_file_name, row.names = 1)
pval<-data["emp_pval"][!is.na(data["emp_pval"])]
pval[pval==0]=discrete_interval

res<-robust.fdr(pval, sides = 1, discrete = T, use8 = T)

result<-res$q
