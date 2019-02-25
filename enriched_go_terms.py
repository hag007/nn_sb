import sys
sys.path.insert(0, '../')

import pandas as pd

from fastsemsim.SemSim import *
from fastsemsim.Ontology import ontologies
from fastsemsim.Ontology import AnnotationCorpus
from fastsemsim.SemSim.SetSemSim import SetSemSim
import matplotlib
matplotlib.use("Agg")

import constants
import utils.go as go



ontology_type = 'GeneOntology'
ignore_parameters = {'ignore': {}}
source_type = 'obo'
source = os.path.join(os.path.join(constants.GO_DIR, constants.GO_FILE_NAME))


QVAL_FIELD="qval"
PVAL_FIELD="pval"
ENRICHED_FILE_NAME="enriched_genes_temp.txt"
TOTAL_FILE_NAME="total_genes_temp.txt"
FDR_TH=0.05
def main(base_folder=os.path.join(constants.OUTPUT_GLOBAL_DIR,"deg")):
    sgd_tables={}
    go_terms = {}
    norm=0
    for cur in [x for x in os.listdir(base_folder)]:
        cancer_type=cur.split("_")[1].split('.')[0]
        sgd_tables[cur]=pd.read_csv(os.path.join(base_folder, cur),sep='\t', index_col=0)
        file(os.path.join(constants.LIST_DIR,ENRICHED_FILE_NAME), 'w+').write("\n".join(list(sgd_tables[cur][sgd_tables[cur][QVAL_FIELD]<FDR_TH].index)))
        file(os.path.join(constants.LIST_DIR, TOTAL_FILE_NAME), 'w+').write("\n".join(list(sgd_tables[cur].index)))
        results=go.check_group_enrichment(ENRICHED_FILE_NAME, TOTAL_FILE_NAME, method="tango", module=cancer_type)
        pd.DataFrame(results).set_index("GO id").to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "go_terms", cancer_type+".tsv"),sep='\t')

main()