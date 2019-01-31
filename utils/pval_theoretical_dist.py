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

if __name__ == "__main__":


        fig, ax = plt.subplots(figsize=(12, 10))
        # ax.set_yscale("log")
        x_axis=np.linspace(0,1, 10000)

        samples_u=np.random.uniform(0,1,10000)
        log_samples_u=-np.log10(samples_u)
        pval=np.array(log_samples_u)
        beta_params = scipy.stats.beta.fit(pval)
        chi2_params = scipy.stats.chi2.fit(pval)
        print beta_params

        fig, ax = plt.subplots(figsize=(12, 10))
        x=np.linspace(0,np.max(pval)+12, 1000)

        sse_chi2=0
        sse_beta=0
        y, x = np.histogram(pval, bins=50, normed=True)
        x = (x + np.roll(x, -1))[:-1] / 2.0
        # x=x[1:]
        pdf_beta = scipy.stats.beta(*beta_params).pdf(x)
        sse_beta = np.sum(np.power(y - pdf_beta, 2.0))
        pdf_chi2 = scipy.stats.chi2(*chi2_params).pdf(x)
        sse_chi2 = np.sum(np.power(y - pdf_chi2, 2.0))
        print sse_chi2
        print sse_beta

        x_axis= x # [0.0001, 0.1, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2]
        plt.plot(x_axis, scipy.stats.beta(*beta_params).pdf(x_axis), color='red', label="beta (error: {})".format(round(sse_beta,2)))
        plt.plot(x_axis, scipy.stats.chi2(*chi2_params).pdf(x_axis), color='blue', label="chi2 (error: {})".format(round(sse_chi2,2)))
        # plt.hist(pval,bins=50,normed=True)
        sns.distplot(pval, norm_hist=True, kde=False, bins=50)
        plt.title("GO pval dist.")
        real_n_term=np.size(pval)
        plt.legend()
        plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "go_pval_dist_uniform_pvals.png"))
        plt.clf()

