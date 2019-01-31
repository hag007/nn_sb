import sys
sys.path.insert(0, '../')

import pandas as pd

from fastsemsim.SemSim import *
from fastsemsim.Ontology import ontologies
from fastsemsim.Ontology import AnnotationCorpus

import matplotlib
matplotlib.use("Agg")

from rpy2.robjects import pandas2ri
pandas2ri.activate()

import constants

import multiprocessing

from utils.daemon_multiprocessing import func_star

import argparse

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

def calc_similarity(mat_adj, i_x, i_y, x, y):
    key="{}_{}".format(i_x,i_y)
    key_inv="{}_{}".format(i_y,i_x)
    if mat_adj[key] != -2: return
    mat_adj[key] = semsim.SemSim(x, y)  # , ResnikSemSim(ontology,ac))
    print mat_adj[key]
    if np.isnan(mat_adj[key]):
        mat_adj[key] = -1
    mat_adj[key_inv] = mat_adj[key]

def main(datasets, algos, pf=10):


    if not os.path.exists(os.path.join(constants.OUTPUT_GLOBAL_DIR,"emp_fdr", "ds_2_alg_scores")):
        try:
            os.makedirs(os.path.join(constants.OUTPUT_GLOBAL_DIR,"emp_fdr", "ds_2_alg_scores"))
        except Exception, e:
            print "error while creating ds_2_alg_scores folder: {}".format(e)

    for cur_ds in datasets:
        constants.update_dirs(DATASET_NAME_u=cur_ds)
        algo_go_sim_score = []
        total_num_genes = []
        algos_signals = []

        for i_algo, cur_algo in enumerate(algos):
            print "current cur_algo: {}".format(cur_algo)
            try:
                emp_results = pd.read_csv(
                    os.path.join(constants.OUTPUT_GLOBAL_DIR, "emp_fdr", "MAX", 
                                 "emp_diff_{}_{}_passed_oob.tsv".format(cur_ds[cur_ds.index("_") + 1:], cur_algo)), sep='\t', index_col=0)
                # emp_results = pd.read_csv(
                #     os.path.join("/home/hag007/Desktop/fdr_terms/fdr_005_i_1000/terms/emp_diff_{}_{}_passed_oob.tsv"
                #                  .format(cur_ds[cur_ds.index("_") + 1:], cur_algo)),
                #     sep='\t', index_col=0)
            except:
                total_num_genes.append(0)
                algos_signals.append(0)
                algo_go_sim_score.append(1)
                continue

            emp_results=emp_results.sort_values(by='emp_rank')
            emp_results_fdr=emp_results.loc[emp_results["passed_oob_permutation_test"].values,:]["GO name"]

            algos_signals.append(len(emp_results_fdr.index))
            all_go_terms = emp_results_fdr.index.values

            try:
                total_num_genes.append(pd.read_csv(
                    os.path.join(constants.OUTPUT_GLOBAL_DIR, constants.DATASET_NAME, cur_algo, "all_modules_general.tsv"),
                    sep="\t")["total_num_genes"][0])
            except:
                total_num_genes.append(0)


            manager=multiprocessing.Manager()
            adj=manager.dict()
            for x in range(len(all_go_terms)):
                for y in range(len(all_go_terms)):
                    adj["{}_{}".format(x,y)]=-2

            params=[]
            for i_x, x in enumerate(all_go_terms):
                for i_y, y in enumerate(all_go_terms):
                    params.append([calc_similarity, [adj, i_x, i_y, x, y]])

            p = multiprocessing.Pool(pf)
            p.map(func_star,params)

            adj_sum=sum([adj["{}_{}".format(x,y)] for x in range(len(all_go_terms)) for y in range(len(all_go_terms)) if adj["{}_{}".format(x,y)]!=-1])
            adj_count=float(len([adj["{}_{}".format(x,y)] for x in range(len(all_go_terms)) for y in range(len(all_go_terms)) if adj["{}_{}".format(x,y)]!=-1]))
            print "adj_sum: ",adj_sum
            print "adj_count: ",adj_count

            p.close()
            p.join()
            file(os.path.join(constants.OUTPUT_GLOBAL_DIR,"emp_fdr", "ds_2_alg_scores", "{}_{}_{}".format(cur_ds,cur_algo, "n_sig.txt")), 'w+').write(str(len(all_go_terms)))
            file(os.path.join(constants.OUTPUT_GLOBAL_DIR, "emp_fdr", "ds_2_alg_scores", "{}_{}_{}".format(cur_ds, cur_algo, "var.txt")), 'w+').write(str(1 - adj_sum / adj_count) if adj_count>0 else str(0) )


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='args')
    parser.add_argument('--datasets', dest='datasets', default="SOC")
    parser.add_argument('--prefix', dest='prefix', default="GE")
    parser.add_argument('--algos', dest='algos', default="jactivemodules_greedy")
    parser.add_argument('--pf', dest='pf', default=10)
    args = parser.parse_args()

    prefix = args.prefix
    datasets=["{}_{}".format(prefix,x) for x in args.datasets.split(",")]
    algos = args.algos.split(",")
    pf=int(args.pf)
    print "test" 
    ds_summary=pd.DataFrame()
    for cur_ds in datasets:
        print "current dataset: {}".format(cur_ds)
        main(datasets=datasets, algos=algos, pf=pf)
