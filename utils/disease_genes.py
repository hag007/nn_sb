import pandas as pd
import constants
import os
import numpy as np
from utils.ensembl2gene_symbol import g2e_convertor


def get_genes_by_disease(disease_name=None, disease_id=None, score=None):
    disease_gene_association = pd.read_csv(os.path.join(constants.DICTIONARIES_DIR, "curated_gene_disease_associations.tsv"), sep="\t")
    if disease_name is not None:
        disease_gene_association=disease_gene_association[disease_gene_association["diseaseName"] == disease_name]
    if disease_id is not None:
        disease_gene_association=disease_gene_association[disease_gene_association["diseaseId"] == disease_id]
    if score is not None:
        disease_gene_association=disease_gene_association[disease_gene_association["score"].astype(np.float) >= score]

    return g2e_convertor(list(disease_gene_association["geneSymbol"]))

def get_all_disease_names():
    disease_gene_association = pd.read_csv(os.path.join(constants.DICTIONARIES_DIR, "curated_gene_disease_associations.tsv"), sep="\t")
    return list(disease_gene_association["diseaseName"].unique())

def get_all_disease_ids():
    disease_gene_association = pd.read_csv(os.path.join(constants.DICTIONARIES_DIR, "curated_gene_disease_associations.tsv"), sep="\t")
    return list(disease_gene_association["diseaseId"])
