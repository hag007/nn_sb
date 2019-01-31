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
from constants import *
from infra import *

pvals=[[0, 0, 0, 0, 0, 0, 0.001, 0.001, 0.001, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.003, 0.004, 0.004, 0.006, 0.013, 0.018, 0.018, 0.021, 0.023, 0.028, 0.028, 0.029, 0.03, 0.031, 0.035, 0.036, 0.037, 0.04, 0.041, 0.042, 0.043, 0.043, 0.043, 0.044, 0.044, 0.045, 0.046, 0.046, 0.049, 0.054, 0.055, 0.056, 0.062, 0.067, 0.069, 0.069, 0.07, 0.071, 0.081, 0.084, 0.085, 0.085, 0.091, 0.093, 0.098, 0.098, 0.098, 0.098, 0.098, 0.098, 0.098, 0.098, 0.105, 0.105, 0.108, 0.111, 0.128, 0.134, 0.138, 0.139, 0.143, 0.144, 0.144, 0.152, 0.172, 0.184, 0.184, 0.187, 0.19, 0.19, 0.196, 0.198, 0.203, 0.207, 0.215, 0.217, 0.226, 0.239, 0.255, 0.281, 0.281, 0.281, 0.285, 0.321, 0.321, 0.327, 0.328, 0.33, 0.33, 0.333, 0.338, 0.348, 0.365, 0.419, 0.419, 0.476, 0.48, 0.544, 0.547, 0.557, 0.592, 0.604, 0.676, 0.706, 0.714, 0.729, 0.747, 0.749, 0.785, 0.785, 0.785, 0.785, 0.785, 0.785, 0.785, 0.785, 0.787, 0.803, 0.807, 0.832, 0.832, 0.842, 0.846, 0.846, 0.847, 0.847, 0.847, 0.85, 0.902, 0.902, 0.942, 0.944, 0.986, 0.986, 0.986, 0.986, 0.986, 0.986, 0.986, 0.986, 0.986, 0.986, 0.986, 0.986, 0.987, 0.995, 0.995, 0.995, 0.996, 0.996, 1, 1, 1, 1, 1, 1, 1]]

for pval in pvals:
    # pval=np.append(pval,np.zeros(5990-np.size(pval)))
    # pval = [10**(-x) for x in pval]
    pval.sort()
    fdr_results = fdrcorrection0(pval, alpha=0.05, method='indep', is_sorted=False)
    true_counter = len([cur for cur in fdr_results[0] if cur == True])
    print "cutoff: {}".format(-np.log10(pval[true_counter]))
    # print list(fdr_results[1])
    print "true hypothesis: {}".format(true_counter)
    print "total hypothesis: {}".format(np.size(fdr_results[0]))