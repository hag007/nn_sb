
import sys
sys.path.append('../')

import constants
import seaborn as sns
sns.set(color_codes=True)
import pandas as pd
import numpy as np
import scipy
import os
import multiprocessing
from multiprocessing import Process
import utils.mixture_dist
import subprocess
import matplotlib.pyplot as plt


from symfit import Parameter, Variable, Fit, sqrt, pi, Equality, Abs, GreaterThan
from sympy.functions.elementary.exponential import exp
from sympy.functions.elementary.piecewise import Piecewise
from sympy.core.add import Add
from sympy.core.mul import Mul
from sympy.core.power import Pow
from sympy.functions.special.gamma_functions import gamma
from symfit.core.objectives import LogLikelihood



if __name__ == "__main__":
    output=pd.read_csv("/home/hag007/Desktop/emp_fdr/SOC/emp_diff_SOC_jactivemodules_greedy_emp_full_500.tsv", sep='\t')
    output=output.loc[output["filtered_pval"].sort_values(ascending=False).index, :]
    output = output.dropna()
    output = output.loc[output['dist_n_samples'].values != '0',:]
    counter=0
    oob_params_counter=0


    for index, cur in output.iterrows():
        if "cell cycle phase transition" != cur["GO name"] : # or "neuro" not in cur["GO name"]: #
            counter += 1
            continue


        data=np.array([float(x) for x in cur["dist_n_samples"][1:-1].split(", ")])
        data=data[data!=0]
        # pval = pval-np.min(pval)

        if np.size(data) == 0:
            data=np.array([0])

        chi2_params = scipy.stats.chi2.fit(data)

        ###
        # k = Parameter("k", values=1.0, max=100, min=0.01)
        sigma2 = Parameter("sigma2", value=1, min=0.5, max=5.0)
        mu = Parameter("mu", value=max(np.max(data) - 1, 5.0), min=5.0, max=max(np.max(data), 5.0))
        V = Parameter("V", value=0.5, min=0.0001, max=1)
        A = Parameter("A", value=1.0, min=0.0001, max=10)
        B = Parameter("B", value=0.0, min=0.0, max=max(min(data) + 1, 1))

        x = Variable()

        model = Add(Mul(V, Piecewise(
            (Mul(exp(-Mul(Add(x, -B), 1 / A)), 1 / A), GreaterThan(Add(x, -B), 0)), (1e-09, True))),
                    Mul(Add(1, -V), Piecewise(((1 / (sqrt(2 * pi * np.abs(sigma2)))) * exp(
                        -(x - mu) ** 2 / (2 * np.abs(sigma2))), GreaterThan(Add(x, -B), 0)), (1e-09, True))))

        ###



        param_counter = 0
        opt_objective_value = np.inf
        optional_mixture_params = [{"v": 0.5}, {"v": 0.75}, {"v": 1.0}, {"v": 0.25}, {"v": 0.0001}]
        mixture_result = None
        opt_mixture_params = None
        manager = multiprocessing.Manager()
        return_list = manager.list()
        prcs = []
        for i, cur_params in enumerate(optional_mixture_params):
            prcs.append(subprocess.Popen("{}/../python27/bin/python {} --data {} --v {}"
                                    .format(constants.dir_path, "mixture_dist_sa.py", ",".join(data.astype(np.str)), cur_params["v"]), shell=True,
                             stdout=subprocess.PIPE, cwd=os.path.join(constants.dir_path, "utils")))
        for i, prc in enumerate(prcs):

            output=prc.stdout.read()

            output=output.split("\n")[-1]
            if output== "None":
                continue
            objective_value=float(output.split("#")[1])
            mixture_params = {cur.split(":")[0] : float(cur.split(":")[1])  for cur in output.split("#")[0].split(",")}

            # mixture_params['k']< k.min or
            if mixture_params['A']< A.min or \
                mixture_params['B']< B.min or mixture_params['sigma2']< sigma2.min or \
                mixture_params['mu'] < mu.min:
                print "params out of bounds: {}".format(mixture_params)
                oob_params_counter+=1
                continue


            print("cur params: {}, cur objective: {}".format(mixture_params, objective_value))
            if opt_objective_value > np.abs(objective_value):
                opt_objective_value = np.abs(objective_value)
                opt_mixture_params = mixture_params

            fig, ax = plt.subplots(figsize=(17, 10))
            ls=np.linspace(0, 100, 1000)
            plt.plot(ls, scipy.stats.chi2(*chi2_params).pdf(ls))

            plt.plot(ls, model(ls,**mixture_params) )


            n_bins='sqrt'
            sns.distplot(data, norm_hist=True, kde=False, bins=n_bins)
            mixture_pval=(1 - mixture_params['V']) * scipy.stats.norm(
                mixture_params['mu'], mixture_params['sigma2']).sf(cur["filtered_pval"]) + \
                                                        mixture_params['V'] * scipy.stats.expon(loc=mixture_params['B'],
                                                                                                scale=mixture_params[
                                                                                                    'A']).sf(
                cur["filtered_pval"])
            # beta_pval=1 - scipy.stats.beta(*beta_params).cdf(cur["filtered_pval"])
            chi2_pval=scipy.stats.chi2(*chi2_params).sf(cur["filtered_pval"])
            plot_title = "terms: {}\n".format(cur["GO name"])
            plot_title+="mixture_params :{}\nobjective:{}\n".format(mixture_params, objective_value)
            plot_title+="chi2 pval: {}, hg pval: {}, mixture pval: {}".format(chi2_pval, cur["filtered_pval"], mixture_pval)
            plt.title(plot_title)
            plt.legend()
            plt.savefig(os.path.join("/media/hag007/Data/bnet/output", "go_pval_dist_SOC_{}_{}.png".format("|".join(cur["GO name"].split("/")), i)))
            plt.clf()
            print(cur["GO name"])
            print("loc={}, scale={}".format(*chi2_params))
            print("loc={}, scale={}, mu={}, sigma2={}, V={}".format(mixture_params['B'],
                                                                          mixture_params['A'], mixture_params['mu'],
                                                                          mixture_params['sigma2'],
                                                                          mixture_params['V']))

            counter+=1

        print "total oob so far: {}".format(oob_params_counter)
            # if counter >50:
            #     break

