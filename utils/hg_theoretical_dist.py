import json
from matplotlib import style
from pandas._libs.parsers import k

style.use("ggplot")
import seaborn as sns
sns.set(color_codes=True)
import logging
sh = logging.StreamHandler()
logger = logging.getLogger("log")
logger.addHandler(sh)
from constants import *
from infra import *
from param_builder import build_gdc_params
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import networkx as nx
import shutil
from pandas.errors import EmptyDataError
from utils.permute_network import EdgeSwapGraph
import scipy
from scipy.stats import hypergeom
import numpy as np

if __name__ == "__main__":

        n_draws=100
        n_set=50
        n_bins=10
        fig, ax = plt.subplots(figsize=(12, 10))
        # ax.set_yscale("log")
        x_axis=np.linspace(0,1, 10000)
        samples = np.random.hypergeometric(n_set, 1900, n_draws, 10000)

        [M, n, N] = [2000, n_set, n_draws]
        rv = hypergeom(M, n, N)
        x = np.arange(0, n + 1)
        pmf_balls = rv.pmf(x)
        plt.plot(x, pmf_balls)


        plt.hist(samples, normed=True, bins=n_bins)
        # sns.distplot(samples, norm_hist=True, kde=False)
        plt.title("GO RV dist.")
        plt.legend()
        plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "hg_rv_dist_uniform_pvals.png"))
        plt.clf()
        pvals=[]
        for i, cur_sample in enumerate(samples):
                print "current iteration: {}".format(i)
                pvals.append(1-hypergeom.cdf(cur_sample, 2000, n_set, n_draws))

        plt.hist(pvals, normed=True, bins=n_bins)


        beta_params = scipy.stats.beta.fit(pvals)
        chi2_params = scipy.stats.chi2.fit(pvals)
        x = np.linspace(0, np.max(pvals) + 12, 1000)

        y, x = np.histogram(pvals, bins=n_bins, normed=True)
        x = (x + np.roll(x, -1))[:-1] / 2.0
        # x=x[1:]
        pdf_beta = scipy.stats.beta(*beta_params).pdf(x)
        sse_beta = np.sum(np.power(y - pdf_beta, 2.0))
        pdf_chi2 = scipy.stats.chi2(*chi2_params).pdf(x)
        sse_chi2 = np.sum(np.power(y - pdf_chi2, 2.0))
        plt.plot(x, scipy.stats.beta(*beta_params).pdf(x), color='red',
                 label="beta (error: {})".format(round(sse_beta, 2)))
        plt.plot(x, scipy.stats.chi2(*chi2_params).pdf(x), color='blue',
                 label="chi2 (error: {})".format(round(sse_chi2, 2)))




        plt.title("GO pval dist.")
        plt.legend()
        plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "hg_pval_dist_uniform_pvals.png"))

        plt.clf()

        log_pvals=-np.log10(pvals)
        plt.hist(log_pvals, normed=True, bins=n_bins)

        beta_params = scipy.stats.beta.fit(log_pvals)
        chi2_params = scipy.stats.chi2.fit(log_pvals)
        x = np.linspace(0, np.max(log_pvals) + 12, 1000)

        y, x = np.histogram(log_pvals, bins=n_bins, normed=True)
        x = (x + np.roll(x, -1))[:-1] / 2.0
        # x=x[1:]
        pdf_beta = scipy.stats.beta(*beta_params).pdf(x)
        sse_beta = np.sum(np.power(y - pdf_beta, 2.0))
        pdf_chi2 = scipy.stats.chi2(*chi2_params).pdf(x)
        sse_chi2 = np.sum(np.power(y - pdf_chi2, 2.0))
        plt.plot(x, scipy.stats.beta(*beta_params).pdf(x), color='red',
                 label="beta (error: {})".format(round(sse_beta, 2)))
        plt.plot(x, scipy.stats.chi2(*chi2_params).pdf(x), color='blue',
                 label="chi2 (error: {})".format(round(sse_chi2, 2)))


        plt.title("GO pval dist.")
        plt.legend()
        plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "hg_log_pval_dist_uniform_pvals.png"))

        plt.clf()
        log_log_pvals=np.log10(-np.log10(pvals))
        log_log_pvals=log_log_pvals*(-1)
        log_log_pvals=log_log_pvals+2
        plt.hist(log_log_pvals, normed=True, bins=n_bins)

        beta_params = scipy.stats.beta.fit(log_log_pvals)
        chi2_params = scipy.stats.chi2.fit(log_log_pvals)
        x = np.linspace(0, np.max(log_log_pvals) + 12, 1000)

        y, x = np.histogram(log_log_pvals, bins=n_bins, normed=True)
        x = (x + np.roll(x, -1))[:-1] / 2.0
        # x=x[1:]
        pdf_beta = scipy.stats.beta(*beta_params).pdf(x)
        sse_beta = np.sum(np.power(y - pdf_beta, 2.0))
        pdf_chi2 = scipy.stats.chi2(*chi2_params).pdf(x)
        sse_chi2 = np.sum(np.power(y - pdf_chi2, 2.0))
        plt.plot(x, scipy.stats.beta(*beta_params).pdf(x), color='red',
                 label="beta (error: {})".format(round(sse_beta, 2)))
        plt.plot(x, scipy.stats.chi2(*chi2_params).pdf(x), color='blue',
                 label="chi2 (error: {})".format(round(sse_chi2, 2)))




        plt.title("GO pval dist.")
        plt.legend()
        plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "hg_log_log_pval_dist_uniform_pvals.png"))


        x=1

