import sys
sys.path.insert(0, '../')

import pandas as pd

from fastsemsim.SemSim import *
from fastsemsim.Ontology import ontologies
from fastsemsim.Ontology import AnnotationCorpus

import matplotlib
matplotlib.use("Agg")

import constants

import multiprocessing

from utils.daemon_multiprocessing import func_star

import argparse

import math
from scipy.cluster.hierarchy import dendrogram, linkage
from matplotlib import pyplot as plt

import scipy.spatial.distance as ssd

ontology_type = 'GeneOntology'
ignore_parameters = {'ignore': {}}
source_type = 'obo'
source = os.path.join(os.path.join(constants.GO_DIR, constants.GO_FILE_NAME))

print "\n######################"
print "# Loading ontology... #"
print "######################\n"

ontology = ontologies.load(source=source, source_type=source_type, ontology_type=ontology_type,
                           parameters=ignore_parameters)

print "\n######################"
print "# Loading Annotation Corpus... #"
print "######################\n"
ac = AnnotationCorpus.AnnotationCorpus(ontology)
ac.parse(os.path.join(constants.GO_DIR, "goa_human.gaf"), "gaf-2.0")
ac.isConsistent()

print "\n#################################"
print "# Annotation corpus successfully loaded."
print "#################################\n"

semsim = GSESAMESemSim(ontology, ac)  # maxSemSim(ontology, ac) #

def calc_similarity(mat_adj, i_x, i_y, x, y, x_score, y_score, norm):
    key="{}_{}".format(i_x,i_y)
    key_inv="{}_{}".format(i_y,i_x)
    if mat_adj[key] != -2: return
    mat_adj[key] = semsim.SemSim(x, y) # * (norm-np.abs(x_score-y_score))/norm# , ResnikSemSim(ontology,ac))
    # print mat_adj[key]
    if np.isnan(mat_adj[key]):
        mat_adj[key] = -1
    mat_adj[key_inv] = mat_adj[key]

def main(base_folder="/home/hag007/Downloads/top_2000_14", pf=5):

    go_terms_tables={}
    go_terms = {}
    norm=0
    for cur in [x for x in os.listdir(base_folder) if not os.path.isdir(x) and x.startswith("GO_") and x.endswith(".xls")]:
        go_terms_tables[cur[3:]]=pd.read_csv(os.path.join(base_folder, cur),sep='\t', index_col=0)
        norm=max(-np.log10(go_terms_tables[cur[3:]]["P-value"].min()),norm)

    for k,v in go_terms_tables.iteritems():
        go_terms[k]=v[v["B"]<=500][v["B"]>=10] #

    df_summary=pd.DataFrame()
    for cur_x in go_terms.keys():
        for cur_y in go_terms.keys():
            if go_terms[cur_y].index.size == 0 or go_terms[cur_x].index.size == 0: continue
            if cur_y==cur_x:
                df_summary.loc[cur_x, cur_y] = 0

            if go_terms.keys().index(cur_y) > go_terms.keys().index(cur_x):

                manager=multiprocessing.Manager()
                adj=manager.dict()
                for x in range(len(go_terms[cur_x].index)):
                    for y in range(len(go_terms[cur_y].index)):
                        adj["{}_{}".format(x,y)]=-2

                params=[]
                for i_x, x in enumerate(go_terms[cur_x].index):
                    for i_y, y in enumerate(go_terms[cur_y].index):
                        params.append([calc_similarity, [adj, i_x, i_y, x, y, -np.log10(go_terms[cur_x].loc[x,"P-value"]), -np.log10(go_terms[cur_y].loc[y,"P-value"]), norm]])

                p = multiprocessing.Pool(pf)
                p.map(func_star,params)

                # adj_sum = sum(
                #     [adj["{}_{}".format(x, y)] for x in range(len(go_terms[cur_x])) for y in range(len(go_terms[cur_y]))
                #      if adj["{}_{}".format(x, y)] != -1])
                # adj_count = float(len(
                #     [adj["{}_{}".format(x, y)] for x in range(len(go_terms[cur_x])) for y in range(len(go_terms[cur_y]))
                #      if adj["{}_{}".format(x, y)] != -1]))

                adj_sum_x_max=[[adj["{}_{}".format(x,y)] for y in range(len(go_terms[cur_y])) if adj["{}_{}".format(x,y)]!=-1 ] for x in range(len(go_terms[cur_x]))]
                adj_sum_x_max=[max(x) for x in adj_sum_x_max]
                adj_sum_y_max = [[adj["{}_{}".format(x, y)] for x in range(len(go_terms[cur_x])) if adj["{}_{}".format(x, y)] != -1] for y in range(len(go_terms[cur_y]))]
                adj_sum_y_max = [max(x) for x in adj_sum_y_max]

                adj_sum_max=adj_sum_x_max+adj_sum_y_max
                adj_sum=sum(adj_sum_max) # - len(set(go_terms[cur_x]).intersection(go_terms[cur_y]))
                adj_count= len(list(np.append(go_terms[cur_x], go_terms[cur_y])))

                print "adj_sum: ",adj_sum
                print "adj_count: ",adj_count

                if adj_sum >0:
                    df_summary.loc[cur_x,cur_y]=1 - adj_sum/adj_count
                    df_summary.loc[cur_y, cur_x]=df_summary.loc[cur_x, cur_y]
                print cur_x, cur_y, go_terms.keys().index(cur_x), go_terms.keys().index(cur_y)
                p.close()
                p.join()

    print df_summary
    df_summary.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "cancer_type_go_distance.tsv"), sep='\t')

    distArray = ssd.squareform(df_summary[df_summary.index.values].values)

    linked = linkage(distArray, method='single', metric='euclidean')

    plt.figure(figsize=(10, 7))
    dendrogram(linked,
               orientation='top',
               labels=[x.split(".")[0] for x in df_summary.index.values],
               distance_sort='descending',
               show_leaf_counts=True)
    plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "hierarchical_go_clustering.png"))

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='args')
    parser.add_argument('--datasets', dest='datasets', default="SOC")
    parser.add_argument('--prefix', dest='prefix', default="GE")
    parser.add_argument('--algos', dest='algos', default="jactivemodules_greedy")
    parser.add_argument('--pf', dest='pf', default=5)
    args = parser.parse_args()

    prefix = args.prefix
    datasets=["{}_{}".format(prefix,x) for x in args.datasets.split(",")]
    algos = args.algos.split(",")
    pf=int(args.pf)
    print "test" 
    ds_summary=pd.DataFrame()
    for cur_ds in datasets:
        print "current dataset: {}".format(cur_ds)
        main(pf=pf)
