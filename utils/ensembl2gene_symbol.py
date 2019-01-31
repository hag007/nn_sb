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
import logging
import constants
sh = logging.StreamHandler()
logger = logging.getLogger("log")
logger.addHandler(sh)
g2e_dict = None
e2g_dict = None

def load_gene_dictionary(gene_list_file_name, gene_list_path=None, source="GDC-TCGA",dataset="melanoma"): #  ="TCGA-SKCM.htseq_counts.tsv"
    if gene_list_path == None:
        gene_list_path = os.path.join(constants.DICTIONARIES_DIR,gene_list_file_name)
    f = open(gene_list_path,'r')
    lines = [l.strip() for l in f]
    f.close()
    return lines


def get_g2e_dictionary():
    lines_dict = load_gene_dictionary(constants.ENSEMBL_TO_GENE_SYMBOLS)

    gene_symbols2ensembl = {}
    for cur in lines_dict:
        splited_line = cur.split()
        if splited_line[0].find('.') > 0:
            limit = splited_line[0].find('.')
        else:
            limit = len(splited_line[0])
        gene_symbols2ensembl[splited_line[1]] = splited_line[0][:limit]
    return gene_symbols2ensembl

def get_e2g_dictionary():
    lines_dict = load_gene_dictionary(constants.ENSEMBL_TO_GENE_SYMBOLS)

    ensembl2gene_symbols = {}
    for cur in lines_dict:
        splited_line = cur.split()
        if splited_line[0].find('.') > 0:
            limit = splited_line[0].find('.')
        else:
            limit = len(splited_line[0])
        ensembl2gene_symbols[splited_line[0][:limit]] = splited_line[1]
    return ensembl2gene_symbols


def e2g_convertor(e_ids):
    if type(e_ids) is str:
        e_ids=[e_ids]

    global e2g_dict
    if e2g_dict is None:
        e2g_dict = get_e2g_dictionary()
    results = []
    for cur in e_ids:
        if e2g_dict.has_key(cur.split(".")[0]):
            results.append(e2g_dict[cur.split(".")[0]])
        else:
            results.append(cur.split(".")[0])
    return results


def g2e_convertor(g_ids):
    if type(g_ids) is str:
        g_ids=[g_ids]

    global g2e_dict
    if g2e_dict is None:
        g2e_dict = get_g2e_dictionary()
    results = []
    for cur in g_ids:
        if g2e_dict.has_key(cur.split(".")[0]):
            results.append(g2e_dict[cur.split(".")[0]])
    return results