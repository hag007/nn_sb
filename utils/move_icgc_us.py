import sys
sys.path.insert(0, "../")
import constants
import os
import shutil

datasets=["KIRC", "KIRP", "LUSC", "LUAD", "COAD", "BRCA", "STAD", "LIHC", "READ", "PRAD", "BLCA", "HNSC", "THCA", "UCEC", "OV", "PAAD"] 
for cur_ds in datasets:
    shutil.move(os.path.join(constants.DATASETS_DIR, "ICGC_{}".format(cur_ds)), os.path.join(constants.DATASETS_DIR, "ICGC_{}_US".format(cur_ds)))
   
   
