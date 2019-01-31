import re
import os
import math
from numpy import *
import numpy.random
from sklearn.datasets import fetch_mldata
import sklearn.preprocessing
import numpy as np
import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import style
style.use("ggplot")
from sklearn import svm
from sklearn import svm
from sklearn.model_selection import GridSearchCV, cross_val_score
from sklearn.model_selection import PredefinedSplit
from sklearn.metrics import accuracy_score
import scipy
from statsmodels.sandbox.stats.multicomp import fdrcorrection0
import time
from matplotlib.ticker import FormatStrFormatter
import math
from ensembl2gene_symbol import e2g_convertor
import constants
import logging
sh = logging.StreamHandler()
logger = logging.getLogger("log")
logger.addHandler(sh)

LABEL_ID = "sample_type.samples"
PRIMARY_TUMOR = "Primary Tumor"
METASTATIC = "Metastatic"


def load_gene_list(gene_list_file_name, gene_list_path=None): #  ="TCGA-SKCM.htseq_counts.tsv"
    if gene_list_path == None:
        gene_list_path = os.path.join(constants.LIST_DIR,gene_list_file_name)
    f = open(gene_list_path,'r')
    lines = [l.strip() for l in f]
    f.close()
    return lines


lines_1 = ["hsa-mir-513a-1", "hsa-mir-513a-2", "hsa-mir-4781", "hsa-mir-4797", "hsa-mir-548b", "hsa-mir-2116", "hsa-mir-510", "hsa-mir-5588", "hsa-mir-548s", "hsa-mir-1468", "hsa-mir-561", "hsa-mir-513b", "hsa-mir-125b-2", "hsa-mir-513c", "hsa-mir-1262", "hsa-mir-6509", "hsa-mir-6507", "hsa-mir-181b-2", "hsa-mir-181b-1", "hsa-mir-873", "hsa-mir-221", "hsa-mir-181a-2", "hsa-mir-508", "hsa-mir-192", "hsa-mir-1910", "hsa-mir-548y", "hsa-mir-506", "hsa-mir-197", "hsa-mir-5683", "hsa-mir-363", "hsa-mir-514a-3", "hsa-mir-876", "hsa-mir-514a-1", "hsa-mir-378g", "hsa-mir-514a-2", "hsa-mir-514b", "hsa-mir-548v", "hsa-mir-509-1", "hsa-mir-30a", "hsa-mir-509-2", "hsa-mir-509-3", "hsa-mir-935", "hsa-mir-507", "hsa-mir-942", "hsa-mir-211", "hsa-let-7c", "hsa-mir-30c-2", "hsa-mir-99a", "hsa-mir-29c", "hsa-mir-194-2", "hsa-mir-651", "hsa-mir-335", "hsa-mir-219a-1", "hsa-mir-222", "hsa-mir-1249", "hsa-mir-1271", "hsa-mir-203a", "hsa-mir-4733"]
lines_2 = ["hsa-mir-4781", "hsa-mir-513a-1", "hsa-mir-513a-2", "hsa-mir-548b", "hsa-mir-510", "hsa-mir-2116", "hsa-mir-5588", "hsa-mir-4797", "hsa-mir-1468", "hsa-mir-548s", "hsa-mir-508", "hsa-mir-513b", "hsa-mir-181a-2", "hsa-mir-125b-2", "hsa-mir-561", "hsa-mir-506", "hsa-mir-181b-2", "hsa-mir-181b-1", "hsa-let-7c", "hsa-mir-6509", "hsa-mir-935", "hsa-mir-1262", "hsa-mir-6507", "hsa-mir-30a", "hsa-mir-192", "hsa-mir-873", "hsa-mir-514a-3", "hsa-mir-514a-1", "hsa-mir-514a-2", "hsa-mir-1910", "hsa-mir-548y", "hsa-mir-211", "hsa-mir-548v", "hsa-mir-514b", "hsa-mir-197", "hsa-mir-378g", "hsa-mir-363", "hsa-mir-507", "hsa-mir-509-1", "hsa-mir-509-2", "hsa-mir-509-3", "hsa-mir-99a", "hsa-mir-876", "hsa-mir-1249", "hsa-mir-5683", "hsa-mir-221", "hsa-mir-219a-1", "hsa-mir-30c-2", "hsa-let-7g", "hsa-mir-942", "hsa-mir-29c", "hsa-mir-222", "hsa-mir-194-2", "hsa-mir-651", "hsa-mir-335", "hsa-mir-26a-2", "hsa-mir-203a", "hsa-mir-26a-1", "hsa-mir-5093", "hsa-mir-6802", "hsa-mir-4733", "hsa-mir-1296"]
lines_uppered = []
for i, cur in enumerate(lines_2):
    lines_2[i] = lines_2[i].upper()

print "genes_1 list size: {}".format(len(lines_1))
print "genes_2 list size: {}".format(len(lines_2))

included_genes = []
for cur in lines_1:
    if cur in lines_2 or cur.upper().split('.')[0] in lines_2:
        included_genes.append(cur)

for cur in included_genes:
    print cur

print "gene count group 1: {} gene count groups two: {}, mutual:{}".format(len(lines_1), len(lines_2),len(included_genes))