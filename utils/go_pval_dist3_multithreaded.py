
import sys
sys.path.append('../')

import seaborn as sns
sns.set(color_codes=True)
import pandas as pd
import numpy as np
import scipy
import os
import multiprocessing
from multiprocessing import Process
import utils.mixture_dist

import matplotlib.pyplot as plt


if __name__ == "__main__":
    output=pd.read_csv("/home/hag007/Desktop/emp_fdr/TNFa_2_2/emp_diff_TNFa_2_jactivemodules_greedy.tsv", sep='\t')
    output=output.loc[output["filtered_pval"].sort_values(ascending=False).index, :]
    output = output.dropna()
    output = output.loc[output['dist_n_samples'].values != '0',:]
    counter=0
    for index, cur in output.iterrows():
        if counter>300: #  or not "chemotaxis" == cur["GO name"] : # or "neuro" not in cur["GO name"]: #
            counter += 1
            continue

        pval=np.array([float(x) for x in cur["dist_n_samples"][1:-1].split(", ")])
        pval=pval[pval!=0]
        # pval = pval-np.min(pval)

        if np.size(pval) == 0:
            pval=np.array([0])
        n_bins='sqrt' #int(np.max(pval))

        # beta_params = scipy.stats.beta.fit(pval)
        chi2_params = scipy.stats.chi2.fit(pval)
        # print beta_params

        fig, ax = plt.subplots(figsize=(17, 10))

        y, x = np.histogram(pval, bins=n_bins, normed=True)
        x = (x + np.roll(x, -1))[:-1] / 2.0

        expected = (1, 1, 1, 1, 1, np.max(pval),1)
        # params, cov = curve_fit(bimodal, x, y, expected , bounds=([0,0,0.5,-300,-300,0,0],[60,60,500,300,300,60,60]), ftol=7)
        # from scipy.special import gamma
        # Draw 100 samples from an exponential distribution with beta=5.5
        # data = np.append(np.random.normal(20, 1, 1000), np.random.exponential(1, 1000))
        data=pval
        # Define the model for an exponential distribution

        param_counter = 0
        opt_objective_values = np.inf
        optional_mixture_params = [{"v": 0.5}, {"v": 0.75}, {"v": 1.0}, {"v": 0.25}, {"v": 0.01}]
        mixture_result = None
        opt_mixture_params = None
        manager = multiprocessing.Manager()
        return_list = manager.list()
        prcs = []
        for i, cur_params in enumerate(optional_mixture_params):
            prcs.append(Process(target=utils.mixture_dist.fit, args=[pval, cur_params['v'], return_list]))
            prcs[-1].start()
            prcs[-1].join()
        # for prc in prcs:
        #     prc.join()


        print(return_list)
        for i, mixture_result in enumerate(return_list):
            if mixture_result is None:
                print("no values for a process.. continue to next one")
                continue
            model, mixture_params, objective_value = mixture_result
            print("cur params: {}, cur objective: {}".format(mixture_params, objective_value))
            if opt_objective_values > np.abs(objective_value):
                opt_objective_values = np.abs(objective_value)
                opt_mixture_params = mixture_params


            ls=np.linspace(0, 100, 1000)
            plt.plot(ls, scipy.stats.chi2(*chi2_params).pdf(ls))

            # plt.hist(data, normed=True, bins=50)

            plt.plot(ls, model(ls,**mixture_params) )

            # plt.title("args: {}".format(opt_mixture_params))
            # plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "test.png"))
            # plt.legend()

            # print "a={}, b={}, k={}, l1={}, s1={}, l2={}, s2={}".format(*params)
            # print chi2_params
            # exit(0)
            # pdf_beta = scipy.stats.beta(*beta_params).pdf(x)
            # sse_beta = np.sum(np.power(y - pdf_beta, 2.0))
            # pdf_chi2 = scipy.stats.chi2(*chi2_params).pdf(x)
            # sse_chi2 = np.sum(np.power(y - pdf_chi2, 2.0))
            # print sse_chi2
            # print sse_beta


            # plt.plot(x, scipy.stats.beta(*beta_params).pdf(x) , color='red', label="beta (error: {})".format(round(sse_beta,2)))
            # plt.plot(x, scipy.stats.chi2(*chi2_params).pdf(x), color='blue', label="chi2 (error: {})".format(round(0,2))) # sse_chi2

            # plt.plot(x, scipy.stats.chi2(params[2],params[3],params[4]).pdf(x), color='blue',
            #          label="chi2 (error: {})".format(round(0, 2)))  # sse_chi2

            # plot_title="terms: {}\npval: {} (real: {})\nk: {},mu: {}, sigma2: {}, A: {}, B: {}"\
            #     .format(cur["GO name"], opt_mixture_params['A'] * scipy.stats.norm(opt_mixture_params['mu'], opt_mixture_params['sigma2']).sf(cur["filtered_pval"]) + \
            # opt_mixture_params['B'] * scipy.stats.expon(opt_mixture_params['k']).sf(cur["filtered_pval"]),
            #              cur["filtered_pval"], opt_mixture_params['k'], opt_mixture_params['mu'], opt_mixture_params['sigma2'], opt_mixture_params['A'], opt_mixture_params['B'])
            plot_title = "terms: {}\n".format(cur["GO name"])
            # print plot_title
            sns.distplot(pval, norm_hist=True, kde=False, bins=n_bins)
            # plt.title("GO pval dist.")
            real_n_term=np.size(pval)
            plt.legend()
            print(cur["GO name"])
            print("k={}, loc={}, scale={}".format(*chi2_params))
            print("k={}, loc={}, scale={}, mu={}, sigma2={}, V={}".format(mixture_params['k'], mixture_params['B'], mixture_params['A'], mixture_params['mu'], mixture_params['sigma2'], mixture_params['V']))
            # plot_title+="k={}, loc={}, scale={}, mu={}, sigma2={}, V={}".format(opt_mixture_params['k'], opt_mixture_params['B'], opt_mixture_params['A'], opt_mixture_params['mu'], opt_mixture_params['sigma2'], opt_mixture_params['V'])
            plot_title+="mixture_params :{}\nobjective:{}".format(mixture_params, objective_value)
            plt.title(plot_title)
            plt.savefig(os.path.join("/media/hag007/Data/bnet/output", "go_pval_dist_TNFa_2_{}_{}.png".format(cur["GO name"], i)))
            plt.clf()
            counter+=1
            # if counter >50:
            #     break

