from scipy.stats import hypergeom

def calc_HG_test(total_gene_list_N, tests_gene_list_B, total_gene_list_n):
    b = len(set(total_gene_list_n).intersection(set(tests_gene_list_B)))
    B = len(set(tests_gene_list_B)) # .intersection(set(total_gene_list_N)))
    N = len(total_gene_list_N)
    n = len(total_gene_list_n)
    print "run HG test with {},{},{},{}".format(b, N, B, n)
    return "{}\t({} {} {} {})".format(hypergeom.sf(b - 1, N, B, n), b, N, B, n)