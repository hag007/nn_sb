import pandas as pd
import os
import constants

dataset="SOC"
algo="jactivemodules_greedy"

path_1="/media/hag007/Data/bnet/output/emp_fdr_5000/fdr/5000_5000/intersection_names_{}_{}.txt".format(dataset, algo)
path_2="/media/hag007/Data/bnet/output/emp_fdr_5000/fdr/1000_200/intersection_names_{}_{}.txt".format(dataset, algo)
path_3="/media/hag007/Data/bnet/output/emp_fdr_5000/fdr/5000_200/intersection_names_{}_{}.txt".format(dataset, algo)



set_1=file(path_1).readlines()
set_2=file(path_2).readlines()
set_3=file(path_3).readlines()

print len(set(set_1))
print len(set(set_2))
print len(set(set_3))

print "==="
print len(set(set_1).intersection(set(set_2)))
print len(set(set_1).intersection(set(set_3)))
print len(set(set_3).intersection(set(set_2)))
print len(set(set_3).intersection(set(set_2)).intersection(set(set_1)))
