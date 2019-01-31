import scipy.special

import matplotlib.pyplot as plt
from matplotlib import style
style.use("ggplot")
import scipy
from scipy.stats import hypergeom
from statsmodels.sandbox.stats.multicomp import fdrcorrection0
import logging
sh = logging.StreamHandler()
logger = logging.getLogger("log")
logger.addHandler(sh)
import constants
from infra import *
import math
MIN_FC_VAL = 1
import time
############################ (2) significance expression and proportion differntiations #############################


def plot_pvalues(y_axis, x_axis, th,output_file_name, is_bigger_better=False):
    n, bins, patches = plt.hist(y_axis, x_axis)
    for c, p in zip(bins, patches):
        if is_bigger_better:
            th_condition = c < th
        else:
            th_condition = c > th

        if th_condition:
            color = 'blue'
        else:
            color = 'red'
        plt.setp(p, 'facecolor', color)
    plt.savefig(os.path.join(constants.OUTPUT_DIR, output_file_name))
    plt.cla()


def plot_pvalues_log_scaled(y_axis, x_axis, th,output_file_name):
    plot_pvalues([-math.log(cur, 10) for cur in y_axis], x_axis, -math.log(th, 10), output_file_name, is_bigger_better=True)


def summarize_genes_proportion(tested_gene_list, total_gene_list, gene_pval_pair, true_counter):
    significant_tested_gene_list = [((i + 1), cur[0], cur[1]) for i, cur in enumerate(gene_pval_pair) if
                                    cur[0] in tested_gene_list and i < true_counter]
    included_tested_genes = [cur for cur in tested_gene_list if
                             cur in total_gene_list or cur[:cur.find('.')] in total_gene_list]
    included_tested_genes_size = len(included_tested_genes)
    significant_tested_gene_list_size = len(significant_tested_gene_list)
    print "total tested genes in true hypothsis: {} out of possible {}".format(significant_tested_gene_list_size, included_tested_genes_size)
    tested_gene_list_size = len(tested_gene_list)
    total_gene_list_size = len(total_gene_list)
    results_table = []
    expected_actual_difference_list = []
    expected_actual_ratio_difference_list = []
    rank_in_total_list_list = []
    z_test_proportion_test_list = []
    for i, cur in enumerate(significant_tested_gene_list):
        rank_in_tesed_list = (i + 1)
        rank_in_total_list = cur[0]
        ensembel_id = cur[1]
        p_val = cur[2]
        expected_quantity = rank_in_total_list * (included_tested_genes_size / (total_gene_list_size * 1.0))
        expected_proportion = included_tested_genes_size / (total_gene_list_size * 1.0)
        actual_quantity = (i + 1)
        actual_proportion = (i + 1) / (cur[0] * 1.0)
        expected_actual_difference = actual_quantity - expected_quantity
        expected_actual_ratio_difference = expected_actual_difference / (expected_quantity * 1.0)
        z_test_proportion_test = (actual_proportion - expected_proportion) / math.sqrt(
            (expected_proportion * (1 - expected_proportion)) / rank_in_total_list)
        results_table.append([rank_in_tesed_list, rank_in_total_list, ensembel_id, p_val,
                              expected_quantity, expected_proportion,
                              actual_quantity, actual_proportion,
                              expected_actual_difference, expected_actual_ratio_difference,
                              z_test_proportion_test])

        expected_actual_difference_list.append(expected_actual_difference)
        expected_actual_ratio_difference_list.append(expected_actual_ratio_difference)
        z_test_proportion_test_list.append(z_test_proportion_test)
        rank_in_total_list_list.append(rank_in_total_list)
    return expected_actual_difference_list, expected_actual_ratio_difference_list, rank_in_total_list_list, z_test_proportion_test_list


def plot_genes_proportion(expected_actual_difference_list, expected_actual_ratio_difference_list, z_test_proportion_test_list, rank_in_total_list_list, total_significant_hypotheses_size, expected_tested_genes_ratio, tested_gene_list_file_name):
    z_score_threshold_two_way = 1.96
    tested_genes_size = len(rank_in_total_list_list)
    y_counter = [min(i, tested_genes_size) for i in range(1,tested_genes_size+1)]
    plt.plot(rank_in_total_list_list[1:], y_counter[1:], label="number of significant values (n)")
    plt.plot(rank_in_total_list_list[1:], expected_actual_difference_list[1:], label="actual-expected significant hypo. difference (n)")
    plt.plot([total_significant_hypotheses_size, total_significant_hypotheses_size], [-20, tested_genes_size + 5], label="True hypotheses threshold")
    plt.plot([0, total_significant_hypotheses_size], [0, tested_genes_size],
             color="gray")
    plt.plot([0, total_significant_hypotheses_size], [0, expected_tested_genes_ratio],
             color="black")
    plt.legend()
    plt.savefig(os.path.join(constants.OUTPUT_DIR, "{}_sum_n".format(tested_gene_list_file_name[:tested_gene_list_file_name.find('.')])))
    plt.cla()
    plt.plot(rank_in_total_list_list[1:], expected_actual_ratio_difference_list[1:], label="actual/expected proportion ratio")
    plt.plot(rank_in_total_list_list[1:], z_test_proportion_test_list[1:], label="z_score")
    plt.plot(rank_in_total_list_list[1:], [z_score_threshold_two_way for i in range(1,tested_genes_size)], label="z_score threshold (two-way)")
    plt.plot([total_significant_hypotheses_size, total_significant_hypotheses_size], [-0.5, 3.5], label="True hypotheses threshold")
    plt.legend()
    plt.savefig(os.path.join(constants.OUTPUT_DIR, "{}_sum_p".format(tested_gene_list_file_name[:tested_gene_list_file_name.find('.')])))

# mHGT DP
def calc_num_of_non_extremer_paths(non_extremer_paths_DP_table, HGTs, mHGT, n, b):
    if n==0 and b==0:
        # non_extremer_paths_DP_table[b][n] = 1
        return 1
    elif b==-1 or b>n or (HGTs[b][n] < mHGT and (b<len(HGTs)-1 or n<len(HGTs[0])-1)):
        # non_extremer_paths_DP_table[b][n] = 0
        return 0
    elif non_extremer_paths_DP_table[b][n] == -1:
        non_extremer_paths_DP_table[b][n] = long(calc_num_of_non_extremer_paths(non_extremer_paths_DP_table, HGTs, mHGT, n-1, b)) + long(calc_num_of_non_extremer_paths(non_extremer_paths_DP_table, HGTs, mHGT, n-1, b-1))


    return non_extremer_paths_DP_table[b][n]

# (2) main
def deg(tested_gene_file_name, total_gene_file_name, gene_expression_file_name, phenotype_file_name, gene_filter_file_name=None, tested_gene_list_path=None, total_gene_list_path=None, gene_expression_path=None, phenotype_path=None, gene_filter_path=None, groups=None, groups_name=None):
    print "about ot analyse: {}".format(tested_gene_file_name)
    # fetch gene expression by gene_id, divided by tumor type11111
    groups_results = load_expression_profile_by_labelling(gene_list_file_name=total_gene_file_name, gene_expression_file_name=gene_expression_file_name, phenotype_file_name=phenotype_file_name, gene_filter_file_name=gene_filter_file_name, tested_gene_path=total_gene_list_path, gene_expression_path=gene_expression_path, phenotype_path=phenotype_path, gene_filter_path=gene_filter_path, groups=groups)
    group_0_expression = groups_results[0]
    group_1_expression = groups_results[1]
    group_0_expression = np.rot90(np.flip(group_0_expression, 1), k=-1, axes=(1,0))
    group_1_expression = np.rot90(np.flip(group_1_expression, 1), k=-1, axes=(1, 0))

    # test pval for significance differentiation between label values (primar vs metastatic)

    pvals = []
    gene_symbols = []
    for  i in range(1,len(group_0_expression)):
        mean_differences = np.average([float(c) for c in group_0_expression[i][1:]]) - np.average([float(c) for c in group_1_expression[i][1:]])

        mean_foldchange = max(np.average([float(c) for c in group_0_expression[i][1:]]),1)/ max(np.average(
            [float(c) for c in group_1_expression[i][1:]]),1)


        cur_pval = scipy.stats.ttest_ind([float(c) for c in group_0_expression[i][1:]], [float(c) for c in group_1_expression[i][1: ]])[1]
        direction = None
        if not math.isnan(cur_pval):
            if mean_differences > 0:
                direction = "downregulated"
            if mean_differences < 0:
                direction = "upregulated"
            pvals.append((group_0_expression[i][0], direction, mean_differences, cur_pval, mean_foldchange))

    pvals.sort(key=lambda x: (x[1], x[3]), reverse=False)
    fdr_results = fdrcorrection0([x[3] for x in pvals], alpha=0.05, method='indep', is_sorted=False)
    pvals = [(cur_pval[0],cur_pval[1],cur_pval[2],cur_pval[3], fdr_results[1][i], cur_pval[4]) for i, cur_pval in enumerate(pvals)]
    true_counter = len([cur for cur in fdr_results[0] if cur == True])
    print "true hypothesis: {}/{}".format(true_counter, np.size(fdr_results[0]))
    # sort gene_id-pval pairs by pval
    with file(os.path.join(constants.OUTPUT_DIR, "deg_{}_{}_{}.txt".format(constants.CANCER_TYPE, groups_name, time.time())), "w+") as f:
        output = ""
        for cur_pval in pvals:
            output+="{}\t{}\t{}\t{}\t{}\t{}\n".format(*cur_pval)
        f.write(output)
        print "pval saved to file"
