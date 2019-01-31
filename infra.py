
import re
import os
from fractions import Fraction
import numpy as np
from utils.stopwatch import Stopwatch
import logging
sh = logging.StreamHandler()
logger = logging.getLogger("log")
logger.addHandler(sh)
import constants
import sys
############################################ infra ########################################

def load_gene_list(gene_list_file_name, gene_list_path=None): #  ="TCGA-SKCM.htseq_counts.tsv"
    if gene_list_path == None:
        gene_list_path = os.path.join(constants.LIST_DIR,gene_list_file_name)
    f = open(gene_list_path,'r')
    lines = [l.strip() for l in f]
    f.close()
    return lines

def load_dictionary(gene_list_file_name, gene_list_path=None): #  ="TCGA-SKCM.htseq_counts.tsv"
    if gene_list_path == None:
        gene_list_path = os.path.join(constants.DICTIONARIES_DIR,gene_list_file_name)
    f = open(gene_list_path,'r')
    lines = [l.strip().split() for l in f]
    f.close()
    return lines

def load_classes(groups_file_name="classes.tsv", gene_list_path=None): #  ="TCGA-SKCM.htseq_counts.tsv"
    if gene_list_path == None:
        gene_list_path = os.path.join(constants.DATA_DIR,groups_file_name)
    f = open(gene_list_path,'r')
    lines = [l.strip().split() for l in f][0]
    f.close()
    return lines

# return gene expression table filtered according an external list with proper orientation (0 degree angle according genes, 90 degree angle according patients)
def load_gene_expression_profile(gene_list_file_name=None, gene_expression_file_name="ge.tsv", gene_filter_file_name=None, gene_list_path=None, gene_expression_path=None, gene_filter_path=None ,by_gene=False, list_mode = "FROM_DISK"):
    stopwatch = Stopwatch()
    stopwatch.start()
    if gene_list_file_name is None:
        gene_list = None
    elif list_mode == "ON_THE_FLY":
        gene_list = gene_list_file_name
    else:
        gene_list = load_gene_list(gene_list_file_name=gene_list_file_name, gene_list_path=gene_list_path)

    if gene_expression_path == None:
        gene_expression_path = os.path.join(constants.DATA_DIR, gene_expression_file_name)





    expression_profiles_filtered=None
    if gene_list == None:
        stopwatch.start()
        f = open(gene_expression_path, 'r')
        expression_profiles_filtered = [l.strip().split() for i, l in enumerate(f)]

    else:
        gene_list = [l.split(".")[0] for i, l in enumerate(gene_list)]
        print stopwatch.stop("done loading gene list")
        # random.shuffle(gene_list)
        # gene_list = gene_list[:400]
        if gene_filter_file_name:
            stopwatch.start()
            filter_gene_list = load_gene_list(gene_list_file_name=gene_filter_file_name, gene_list_path=gene_filter_path)
            gene_list = [cur for cur in gene_list if cur in filter_gene_list]
            print stopwatch.stop("done filter gene list")

        if gene_expression_path == None:
            gene_expression_path = os.path.join(constants.DATA_DIR, gene_expression_file_name)
            stopwatch.start()
        f = open(gene_expression_path,'r')
        expression_profiles_filtered = [l.strip().split() for i, l in enumerate(f) if i==0 or l[:l.strip().find('\t')].split(".")[0] in gene_list]

    expression_profiles_filtered = [x for x in expression_profiles_filtered if len(x) == len(expression_profiles_filtered[0])]
    # or l.strip()[0:l.strip().find('\t')] in gene_list or l.strip()[0:l.strip().find('\t')].split(".")[0] in gene_list
    f.close()
    print stopwatch.stop("done filter gene expression")
    if not by_gene:
        stopwatch.start()
        expression_profiles_filtered = np.flip(np.rot90(expression_profiles_filtered, k=1, axes=(1,0)),1)
        print stopwatch.stop("done rotate gene expression")

    return expression_profiles_filtered


def load_gene_expression_profile_by_genes(gene_list_file_name=None, gene_expression_file_name="ge.tsv", gene_filter_file_name=None, gene_list_path=None, gene_expression_path=None, gene_filter_path=None, list_mode="FROM_DISK"):
    return load_gene_expression_profile(gene_list_file_name=gene_list_file_name, gene_expression_file_name=gene_expression_file_name, gene_filter_file_name=gene_filter_file_name, gene_list_path=gene_list_path, gene_expression_path=gene_expression_path, gene_filter_path=gene_filter_path, by_gene=True, list_mode=list_mode)

def load_gene_expression_profile_by_patients(gene_list_file_name=None, gene_expression_file_name="ge.tsv", gene_filter_file_name=None, gene_list_path=None, gene_expression_path=None, gene_filter_path=None, list_mode="FROM_DISK"):
    return load_gene_expression_profile(gene_list_file_name=gene_list_file_name, gene_expression_file_name=gene_expression_file_name, gene_filter_file_name=gene_filter_file_name, gene_list_path=gene_list_path, gene_expression_path=gene_expression_path, gene_filter_path=gene_filter_path,by_gene=False, list_mode=list_mode)


def load_phenotype_data(phenotype_file_name, phenotype_list_path=None, source="GDC-TCGA",dataset="melanoma"):
    if not phenotype_list_path:
        phenotype_list_path = constants.TCGA_DATA_DIR
    f = open(os.path.join(phenotype_list_path,phenotype_file_name), 'r')
    phenotype_profile = [re.split("[\t]", l.rstrip('\n')) for l in f]
    f.close()
    return phenotype_profile

def load_survival_data(survival_file_name, survival_list_path=None, source="GDC-TCGA",dataset="melanoma"):
    if not survival_list_path:
        survival_list_path = constants.TCGA_DATA_DIR
    f = open(os.path.join(survival_list_path,survival_file_name), 'r')
    survival_profile = [l.strip().split('\t') for l in f]
    if constants.PHENOTYPE_FORMAT=="TCGA":
        survival_profile = [survival_profile[0]] + [[x[0][:-1]] + x[1:] for x in survival_profile[1:]]

    f.close()
    return survival_profile

def load_mutation_data(mutation_file_name, mutation_list_path=None, source="GDC-TCGA", dataset="melanoma"):
    if not mutation_list_path:
        mutation_list_path = os.path.join(constants.TCGA_DATA_DIR,mutation_file_name)
    f = open(mutation_list_path, 'r')

    mutation_profile = []
    mutation_profile.append(f.readline().strip().split('\t'))
    default_list = ["" for x in mutation_profile[0]]
    mutation_profile = mutation_profile + [(l.strip().split('\t') + default_list)[:len(default_list)] for l in f]
    f.close()
    return mutation_profile


def divided_patient_ids_by_label(phenotype_list_file_name, phenotype_list_path=None, labels=None, label_values=None, groups=None):
    thismodule = sys.modules[__name__]
    if not groups and not labels:
        groups = [{"sample_type.samples" :{"type": "string", "value": ["Primary Tumor"]}},
                  {"sample_type.samples": {"type": "string", "value": ["Metastatic"]}}]
    elif not groups and any(labels):
        divided_patient_ids_by_label_old(phenotype_list_file_name, phenotype_list_path, labels, label_values)
        return
    phenotype_data_formatted = load_phenotype_data(phenotype_list_file_name, phenotype_list_path)
    headers = phenotype_data_formatted[0]
    phenotype_profiles = phenotype_data_formatted[1:]
    patients_by_labeling = []
    for i in range(len(groups)):
        patients_by_labeling.append([])
    label_indices = []
    duplicates_counter = 0
    try:
        for i, pp in enumerate(phenotype_profiles):
            for j, cur_group in enumerate(groups):
                dup=0
                is_hold_constraits = True
                for k,v in cur_group.iteritems():
                    if k in constants.FILTER_KEYWORDS: continue
                    if len(pp) <= headers.index(k):
                        is_hold_constraits = False
                        continue
                    if v['type'] == "string":
                        if not any([pp[headers.index(k)] == cur for cur in v["value"]]):
                            is_hold_constraits = False
                    elif v['type'] == "number":
                        if not getattr(thismodule, "num_op_{}".format(v['op']))(pp[headers.index(k)], v["value"]):
                                is_hold_constraits = False
                if is_hold_constraits:
                    patients_by_labeling[j].append(pp[0])
                    dup+=1
            if dup > 1:
                duplicates_counter+=1
        print "number of duplicated patients: {}".format(duplicates_counter)
    except ValueError, e:
        print ValueError, e
        pass
    return (patients_by_labeling)


def divided_patient_ids_by_label_old(phenotype_list_file_name, phenotype_list_path=None, labels=None, label_values=None,
                                 group_0=None, group_1=None):
    if not labels:
        labels = [constants.LABEL_ID]
    if type(labels) is not list:
        labels = [labels]
    if not label_values:
        label_values = [[constants.PRIMARY_TUMOR], [constants.METASTATIC]]
    if type(label_values[0]) is not list:
        label_values = [[label_values[0]], [label_values[1]]]
    phenotype_data_formatted = load_phenotype_data(phenotype_list_file_name, phenotype_list_path)
    headers = phenotype_data_formatted[0]
    phenotype_profiles = phenotype_data_formatted[1:]
    label_0_patients = []
    label_1_patients = []
    label_indices = []
    for label in labels:
        label_indices.append([i for i, v in enumerate(headers) if v == label][0])
    for pp in phenotype_profiles:
        if all([pp[cur_idx] == label_values[0][i] for i, cur_idx in enumerate(label_indices)]):
            label_0_patients.append(pp[0])
        elif all([pp[cur_idx] == label_values[1][i] for i, cur_idx in enumerate(label_indices)]):
            label_1_patients.append(pp[0])
    return (label_0_patients, label_1_patients)

# def transform_to_map(tbl):
#     map = {}
#     for cur in tbl[1:]:
#         map[cur[0]] = cur[1:]
#     return map

# load expression profile filtered by an external genes list and divided according tumor type label
def load_expression_profile_by_labelling(gene_list_file_name, gene_expression_file_name, phenotype_file_name, label=None, label_values=None, gene_filter_file_name=None, tested_gene_path=None, gene_expression_path=None, phenotype_path=None, gene_filter_path=None, groups=None):

    expression_profiles_formatted = load_gene_expression_profile_by_patients(gene_list_file_name, gene_expression_file_name, gene_filter_file_name=gene_filter_file_name, gene_list_path=tested_gene_path, gene_expression_path=gene_expression_path, gene_filter_path=gene_filter_path)
    patients_by_labeling = divided_patient_ids_by_label(phenotype_file_name, phenotype_path, label, label_values, groups)

    expression_profiles_by_labeling = []
    for i in patients_by_labeling:
        expression_profiles_by_labeling.append([])
    logger.info("about to split expression by primary tumor and metastatic")

    logger.info("expression_profile size: {},{}".format(*np.shape(expression_profiles_formatted)))
    for i,cur in enumerate(expression_profiles_formatted):
        if i==0: # that is, vertical headers
            for cur_group in expression_profiles_by_labeling:
                cur_group.append(cur)
        cur_id = expression_profiles_formatted[i][0]
        if not str.isdigit(cur_id[-1]) and constants.PHENOTYPE_FORMAT == "TCGA":
            cur_id = cur_id[:-1]
        patient_found = False
        for k, cur_group in enumerate(expression_profiles_by_labeling):
            if cur_id in patients_by_labeling[k]:
                expression_profiles_by_labeling[k].append(cur)
                patient_found = True
        if not patient_found:
            logger.info("no labeling option were found for {}".format(expression_profiles_formatted[i][0]))

    print "done split expression"

    return expression_profiles_by_labeling

# def nCk(n,k):
#   return long( reduce(mul, (Fraction(n-i, i+1) for i in range(k)), 1) )

def save_sets(st,fl_name):
    fl = open(fl_name,'w+')
    for cur in st:
        for n in cur:
            fl.write(str(n)+constants.SEPARATOR)
        fl.write('\n')

def load_sets(fl_name):
    lst = []
    fl = open(fl_name,'r')
    for cur in fl.readlines():
        lst.append(cur.split(constants.SEPARATOR)[:-1])
        lst[-1][-1] = float(lst[-1][-1])
    return lst


def labels_assignments(meta_groups, phenotype_file_name, patients_list):
    labels_assignment = []
    for i, cur_groups in enumerate(meta_groups):
        labeled_patients = divided_patient_ids_by_label(phenotype_file_name, groups=cur_groups)
        cur_labeled = np.array(patients_list)
        for j, cur_patients_group in enumerate(labeled_patients):
            cur_labeled[np.in1d(cur_labeled, cur_patients_group)] = j + 1
        cur_labeled[~np.core.defchararray.isdigit(cur_labeled)] = 0
        labels_assignment.append(cur_labeled.astype(np.int32))
    return labels_assignment


def separate_headers(dataset, is_numbers=True):
    dataset = np.array(dataset)
    dataset_headers_columns = dataset[0][1:]
    dataset_headers_rows = dataset[1:, 0]
    dataset = dataset[1:]
    dataset = dataset[:, 1:]
    if is_numbers:
        dataset = dataset.astype(np.float64)
    return dataset_headers_rows, dataset_headers_columns, dataset

def load_integrated_mutation_data(mutation_file_name,
                            survival_file_name, phenotype_file_name, gene_filter_file_name=None, filter_expression=None,
                            meta_groups=None, phenotype_labels_heatmap=None, var_th_index=None):


    cache_path = os.path.join(constants.CACHE_DIR, "datasets",
                              "datasets_{}".format(mutation_file_name[:mutation_file_name.rindex(".")]))

    if constants.USE_CACHE and os.path.exists(cache_path):
        print "loading datasets from cache"
        mutations_headers_rows =  np.load(os.path.join(cache_path,"header_rows.npy"))
        mutations_headers_columns =  np.load(os.path.join(cache_path,"header_columns.npy"))
        mutation_dataset =  np.load(os.path.join(cache_path,"data.npy"))
    else:
        print "loading datasets from files"
        mutation_dataset = np.array(load_mutation_data(mutation_file_name))
        print "separating dataset headers"
        mutations_headers_rows, mutations_headers_columns, mutation_dataset = separate_headers(
            mutation_dataset, is_numbers=False)

        myfunc_vec = np.vectorize(lambda x: x if str.isdigit(x[-1]) else x[:-1])
        if constants.PHENOTYPE_FORMAT == "TCGA":
            # col_original_len = len(mutations_headers_rows)
            mutations_headers_rows = myfunc_vec(mutations_headers_rows)
            # duplicated_i = [i for i, x in enumerate(tested_gene_expression_headers_columns) if
            #                 x in tested_gene_expression_headers_columns[i + 1:]]
            # tested_gene_expression_headers_columns = tested_gene_expression_headers_columns[
            #     [False if a in duplicated_i else True for a in range(col_original_len)]]
            # mutation_dataset = mutation_dataset[[False if a in duplicated_i else True for a in range(col_original_len)],:]





    if constants.USE_CACHE:
        if not os.path.exists(cache_path):
            os.makedirs(cache_path)
        print "saving data to cahce"
        np.save(os.path.join(cache_path,"header_rows.npy"), mutations_headers_rows)
        np.save(os.path.join(cache_path,"header_columns.npy"), mutations_headers_columns)
        np.save(os.path.join(cache_path,"data.npy"), mutation_dataset)

    survival_dataset = np.array(load_survival_data(survival_file_name, survival_list_path=None))
    phenotype_dataset = np.array(load_phenotype_data(phenotype_file_name, phenotype_list_path=None))

    if np.shape(mutation_dataset)[0] < 2:
        print "no mutations were found for the specific gene list {}. skipping...".format(
            mutation_file_name.split(".")[0])
        return None

    if filter_expression is not None:
        filtered_patients =[ y for x in divided_patient_ids_by_label(phenotype_file_name, groups=filter_expression) for y in x]
        print "number of filtered patients from phenotypes: {}".format(len(filtered_patients))
    else:
        filtered_patients = np.append(mutations_headers_rows, survival_dataset[1:, 0])

    mutation_dataset, mutations_headers_rows = filter_patients_dataset_by_patients(filtered_patients, mutations_headers_rows, mutation_dataset)
    if np.shape(mutation_dataset)[0] == 1:
        print "no expressions were found after filtering by labels {}. skipping...".format(filter_expression)
        return None

    survival_dataset = filter_survival_by_patients(filtered_patients, survival_dataset)
    if np.shape(survival_dataset)[0] == 1:
        print "no survival were found after filtering by labels {}. skipping...".format(filter_expression)
        return None

    print "total patients taken into account: {}".format(
        len([x for x in survival_dataset[:, 0] if x in mutations_headers_rows]))

    unique_mutations_headers_rows = np.unique(mutations_headers_rows)
    labels_assignment = None
    if meta_groups is not None:
        labels_assignment = labels_assignments(meta_groups, phenotype_file_name, unique_mutations_headers_rows)

    phenotype_heatmap = None
    if phenotype_labels_heatmap is not None:
        phenotype_heatmap = phenotype_dataset[1:, np.in1d(phenotype_dataset[0], phenotype_labels_heatmap)]
        phenotype_heatmap = phenotype_heatmap[[np.where(phenotype_dataset[1:,0]==x)[0][0] for x in [y for y in  unique_mutations_headers_rows if y in filtered_patients]]]
        phenotype_heatmap[phenotype_heatmap==""] = "-1"
        phenotype_heatmap = phenotype_heatmap.astype(np.float32)
        phenotype_heatmap[phenotype_heatmap > 5] = 5
    return (mutation_dataset, mutations_headers_rows, mutations_headers_columns, labels_assignment, survival_dataset, phenotype_heatmap)


def load_integrated_ge_data(tested_gene_list_file_name=None, total_gene_list_file_name=None, gene_expression_file_name="ge.tsv",
                            survival_file_name=None, phenotype_file_name=None, gene_filter_file_name=None, filter_expression=None,
                            meta_groups=None, var_th_index=None):

    cache_path = os.path.join(constants.CACHE_DIR, "datasets",
                              "datasets_{}_{}".format(gene_expression_file_name.split(".")[0],
                                                      tested_gene_list_file_name[:tested_gene_list_file_name.rindex(".")]))

    tested_gene_expression_headers_rows = None
    tested_gene_expression_headers_columns = None
    tested_gene_expression = None
    if constants.USE_CACHE and os.path.exists(cache_path):
        print "loading datasets from cache"
        tested_gene_expression_headers_rows =  np.load(os.path.join(cache_path,"header_rows.npy"))
        tested_gene_expression_headers_columns =  np.load(os.path.join(cache_path,"header_columns.npy"))
        tested_gene_expression =  np.load(os.path.join(cache_path,"data.npy"))
    else:
        print "loading datasets from files"
        tested_gene_expression = np.array(
            load_gene_expression_profile_by_genes(tested_gene_list_file_name, gene_expression_file_name,
                                                  gene_filter_file_name))
        print "separating dataset headers"
        tested_gene_expression_headers_rows, tested_gene_expression_headers_columns, tested_gene_expression = separate_headers(
            tested_gene_expression)
        myfunc_vec = np.vectorize(lambda x: x if str.isdigit(x[-1]) else x[:-1])
        if constants.PHENOTYPE_FORMAT=="TCGA":
            col_original_len = len(tested_gene_expression_headers_columns)
            tested_gene_expression_headers_columns = myfunc_vec(tested_gene_expression_headers_columns)
            duplicated_i = [i for i, x in enumerate(tested_gene_expression_headers_columns) if
             x in tested_gene_expression_headers_columns[i + 1:]]
            tested_gene_expression_headers_columns=tested_gene_expression_headers_columns[
                [False if a in duplicated_i else True for a in range(col_original_len)]]
            tested_gene_expression=tested_gene_expression[:,[False if a in duplicated_i else True for a in range(col_original_len)]]



    if not os.path.exists(cache_path):
        os.makedirs(cache_path)
    print "saving data to cahce"
    np.save(os.path.join(cache_path,"header_rows.npy"), tested_gene_expression_headers_rows)
    np.save(os.path.join(cache_path,"header_columns.npy"), tested_gene_expression_headers_columns)
    np.save(os.path.join(cache_path,"data.npy"), tested_gene_expression)

    survival_dataset=None
    if os.path.exists(os.path.join(constants.DATA_DIR,survival_file_name)):
        print "loading data survival data"
        survival_dataset = np.array(load_survival_data(survival_file_name, survival_list_path=None))
    else:
        "no survival data. continue..."

    if np.shape(tested_gene_expression)[0] < 1:
        print "no expressions were found for the specific gene list {}. skipping...".format(
            tested_gene_list_file_name.split(".")[0])
        return None


    if filter_expression is not None:
        filtered_patients = [y for x in divided_patient_ids_by_label(phenotype_file_name, groups=filter_expression) for y in x]
        print "number of filtered patients from phenotypes: {}".format(len(filtered_patients))
    else:
        print "no filter applied"
        filtered_patients = tested_gene_expression_headers_columns

    tested_gene_expression, tested_gene_expression_headers_columns = filter_genes_dataset_by_patients(filtered_patients, tested_gene_expression_headers_columns, tested_gene_expression)
    if np.shape(tested_gene_expression)[1] == 1:
        print "no expressions were found after filtering by labels {}. skipping...".format(filter_expression)
        return None

    num_filtered_patients = len(tested_gene_expression_headers_columns)
    if survival_dataset is not None:
        filtered_patients = np.append(tested_gene_expression[0, 1:], survival_dataset[1:, 0])
        survival_dataset = filter_survival_by_patients(filtered_patients, survival_dataset)
        if np.shape(survival_dataset)[0] == 1:
            print "no survivors were found after filtering by labels {}. skipping...".format(filter_expression)
            return None

        num_filtered_patients = len([x for x in survival_dataset[:, 0] if x in tested_gene_expression_headers_columns])

    print "total patients taken into account: {}".format(num_filtered_patients)



    labels_assignment = None
    if meta_groups is not None:
        print "clustering patients by groups"
        labels_assignment = labels_assignments(meta_groups, phenotype_file_name,
                                               tested_gene_expression_headers_columns)

    if var_th_index is not None:
        print "filtering top vars"
        tested_gene_expression, tested_gene_expression_headers_rows, tested_gene_expression_headers_columns = filter_top_var_genes(
        tested_gene_expression, tested_gene_expression_headers_columns,
        tested_gene_expression_headers_rows, var_th_index)
    else:
        print "skipping filter top vars"

    tmp = tested_gene_expression_headers_rows
    tested_gene_expression_headers_rows = tested_gene_expression_headers_columns
    tested_gene_expression_headers_columns = tmp
    tested_gene_expression = np.rot90(np.flip(tested_gene_expression, 1), k=-1, axes=(1, 0))

    return (tested_gene_expression, tested_gene_expression_headers_rows, tested_gene_expression_headers_columns, labels_assignment, survival_dataset)


def filter_top_var_genes(tested_gene_expression, tested_gene_expression_headers_columns,
                         tested_gene_expression_headers_rows, var_th_index):
    row_var = np.var(tested_gene_expression, axis=1)
    row_var_sorted = np.sort(row_var)[::-1]
    if var_th_index is None:
        var_th_index = len(row_var_sorted) - 1
    row_var_th = row_var_sorted[var_th_index]
    row_var_masked_indices = np.where(row_var_th > row_var)[0]
    gene_expression_top_var = np.delete(tested_gene_expression, row_var_masked_indices, axis=0)

    gene_expression_top_var_headers_rows = np.delete(tested_gene_expression_headers_rows, row_var_masked_indices,
                                                     axis=0)
    gene_expression_top_var_headers_columns = tested_gene_expression_headers_columns

    return gene_expression_top_var, gene_expression_top_var_headers_rows, gene_expression_top_var_headers_columns



def filter_survival_by_patients(filtered_patients, survival_dataset):
    filtered_survival_bool = np.in1d(survival_dataset[:, 0], filtered_patients)
    # filtered_survival_bool[0] = True
    print "Total n patients in survival before filtering: {}".format(np.shape(survival_dataset)[0] - 1)
    survival_dataset = survival_dataset[filtered_survival_bool, :]
    print "Total n patients in survival after filtering: {}".format(np.shape(survival_dataset)[0] - 1)
    return survival_dataset


def filter_patients_dataset_by_patients(filtered_patients, total_patients, dataset_by_patients):
    filtered_gene_expression_bool = np.in1d(total_patients, filtered_patients)
    # filtered_gene_expression_bool[0] = True
    print "Total n patients in expression before filtering: {}".format(np.shape(dataset_by_patients)[0] - 1)
    dataset_by_patients = dataset_by_patients[filtered_gene_expression_bool]
    total_patients = total_patients[filtered_gene_expression_bool]
    print "Total n patients in expression after filtering: {}".format(np.shape(dataset_by_patients)[0] - 1)
    return dataset_by_patients, total_patients

def filter_genes_dataset_by_patients(filtered_patients, total_patients, dataset_by_genes):
    filtered_gene_expression_bool = np.in1d(total_patients, filtered_patients)
    # filtered_gene_expression_bool[0] = True
    print "Total n patients in expression before filtering: {}".format(np.shape(dataset_by_genes)[1] - 1)
    dataset_by_genes = dataset_by_genes[:, filtered_gene_expression_bool]
    total_patients = total_patients[filtered_gene_expression_bool]
    print "Total n patients in expression after filtering: {}".format(np.shape(dataset_by_genes)[1] - 1)
    return dataset_by_genes, total_patients

def num_op_gt(ds_value, q_value):
    try:
        return str.isdigit(ds_value) and float(ds_value) > float(q_value)
    except ValueError:
        return False

def num_op_lt(ds_value, q_value):
    try:
        return str.isdigit(ds_value) and float(ds_value) < float(q_value)
    except ValueError:
        return False
