import simplejson as json
import os

dir_path = os.path.dirname(os.path.realpath(__file__))
PATH_TO_CONF = "env/config/conf.json"
config_json = json.load(file(os.path.join(dir_path, PATH_TO_CONF)))

REPO_DIR = os.path.dirname(os.path.realpath(__file__))
SH_DIR = os.path.join(REPO_DIR, "sh","scripts")

ALGO_BASE_DIR = config_json['ALGO_BASE_PROFILE']

SERVER_MODE = False
REPORTS = True
DISEASE_MODE = False
HG_MODE = True
ALGO_HG_MODE=False
EMB_MODE = True
USE_CACHE = False
PHENOTYPE_FORMAT = "GDC"
DATASET_TYPE = "GDC"
CANCER_TYPE = "SKCM"
BASE_PROFILE= config_json['BASE_PROFILE']
DATASET_NAME = "TNFa_2"
DATASETS_DIR= os.path.join(BASE_PROFILE, "datasets")
NETWORKS_DIR = os.path.join(BASE_PROFILE, "networks")
TEMPLATES_DIR = os.path.join(BASE_PROFILE, "templates")
DATASET_DIR = os.path.join(DATASETS_DIR, DATASET_NAME)
DATA_DIR = os.path.join(DATASET_DIR, "data")
CACHE_DIR = os.path.join(DATASET_DIR, "cache")
DICTIONARIES_DIR = os.path.join(BASE_PROFILE, "dictionaries")
OUTPUT_DIR = os.path.join(DATASET_DIR, "output")
OUTPUT_GLOBAL_DIR = os.path.join(BASE_PROFILE, "output")
TCGA_DATA_DIR = os.path.join(DATASET_DIR, "data")
GO_DIR = os.path.join(BASE_PROFILE, "GO")
CACHE_GLOBAL_DIR = os.path.join(BASE_PROFILE, "cache_global")
LIST_DIR = os.path.join(BASE_PROFILE, "list")
REPOS_DIR = os.path.join(BASE_PROFILE, "repos")
RAW_DIR = os.path.join(BASE_PROFILE, "raw")

DEG_DESEQ = "deseq"
DEG_EDGER = "edger"
DEG_T = "t"
PREDEFINED_SCORE = "predefined_scored"
IS_PVAL_SCORES = True

SEPARATOR = "@%@"

LABEL_ID = "sample_type.samples"
PRIMARY_TUMOR = "Primary Tumor"
METASTATIC = "Metastatic"

LABELS_NORMAL = "labels_normal"
LABELS_SHUFFLE = "labels_shuffle"
LABELS_RANDOM = "labels_random"
LABELS_ALTERNATED = "labels_alternated"
LABELS_INVERTED = "labels_inverted"

ENSEMBL_TO_GENE_SYMBOLS = "ensembl2gene_symbol.txt"
ENSEMBL_TO_ENTREZ = "ensembl2entrez.txt"

GO_OBO_URL = 'http://purl.obolibrary.org/obo/go/go-basic.obo'
GO_ASSOCIATION_GENE2GEO_URL = 'https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz'
GO_FILE_NAME = 'go_bp.obo' #'go-basic.obo'
GO_ASSOCIATION_FILE_NAME = "gene2go"

NUM_GTE = "gte"
NUM_GT = "gt"
NUM_EQ = "eq"
NUM_LTE = "lte"
NUM_LT = "lt"
NUM_NE = "ne"
NUM_ALL_OPS = [NUM_EQ, NUM_GTE, NUM_GT, NUM_LTE, NUM_LT, NUM_NE]


FROM_DISK = "FROM_DISK"
ON_THE_FLY = "ON_THE_FLY"
FILTER_KEYWORDS = ["_label", "_name"]
ALL_CANCER_TYPES = ["ESCA", "LAML", "ACC", "CHOL", "BLCA", "BRCA", "CESC", "COAD", "UCEC", "GBM", "HNSC", "KICH", "KIRC", "KIRP", "DLBC", "LIHC", "LGG", "LUAD", "LUSC", "SKCM", "MESO", "UVM", "PANCAN", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "STAD", "TGCT", "THYM", "THCA", "UCS"]
ALL_TUMOR_TYPES = ["Primary Tumor", "Metastatic", "Additional - New Parimary", "Additional Metatatic", "Primary Blood Derived Cancer - Peripheral Blood", "Blood Derived Cancer - Bone Marrow, Post-treatment", "Primary Blood Derived Cancer - Bone Marrow", "Recurrent Blood Derived Cancer - Peripheral Blood", "Recurrent Tumor"]
def update_dirs(BASE_DIR=config_json["BASE_PROFILE"], DATASET_NAME_u=None, DATASET_TYPE_u ="GDC-TCGA", CANCER_TYPE_u ="SKCM"):

    global BASE_PROFILE
    global CACHE_DIR
    global OUTPUT_DIR
    global LIST_DIR
    global TCGA_DATA_DIR
    global DICTIONARIES_DIR
    global CANCER_TYPE
    global DATASET_TYPE
    global DATASETS_DIR
    global GO_DIR
    global CACHE_GLOBAL_DIR
    global OUTPUT_GLOBAL_DIR
    global DATA_DIR
    global DATASET_DIR
    global RAW_DIR
    global DATASET_NAME

    BASE_PROFILE=BASE_DIR
    DATASET_TYPE = DATASET_TYPE_u
    CANCER_TYPE = CANCER_TYPE_u
    if DATASET_NAME_u is None:
        DATASET_NAME_u = CANCER_TYPE_u
    DATASET_NAME = DATASET_NAME_u


    DATASETS_DIR= os.path.join(BASE_PROFILE, "datasets")
    DATASET_DIR = os.path.join(DATASETS_DIR, DATASET_NAME)  # ""{}/{}/".format(DATASET_TYPE,CANCER_TYPE)
    CACHE_DIR = os.path.join(DATASET_DIR, "cache")
    DICTIONARIES_DIR = os.path.join(BASE_PROFILE, "dictionaries")
    OUTPUT_DIR = os.path.join(DATASET_DIR, "output")
    OUTPUT_GLOBAL_DIR = os.path.join(BASE_PROFILE, "output")
    TCGA_DATA_DIR = os.path.join(DATASET_DIR, "data")
    DATA_DIR = os.path.join(DATASET_DIR, "data")
    GO_DIR = os.path.join(BASE_PROFILE, "GO")
    LIST_DIR = os.path.join(BASE_PROFILE, "list")
    CACHE_GLOBAL_DIR = os.path.join(BASE_PROFILE, "cache_global")
    RAW_DIR = os.path.join(BASE_PROFILE, "raw")
update_dirs()