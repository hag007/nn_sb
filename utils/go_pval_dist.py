from matplotlib import style

style.use("ggplot")
import seaborn as sns
sns.set(color_codes=True)
import logging
sh = logging.StreamHandler()
logger = logging.getLogger("log")
logger.addHandler(sh)
from infra import *
import pandas as pd
import scipy
import mixture_dist
from symfit import Parameter, Variable, Fit, sqrt, pi, Equality, Abs, GreaterThan
from sympy.functions.elementary.exponential import exp
from sympy.functions.elementary.miscellaneous import Max
from sympy.core.add import Add
from sympy.core.mul import Mul
from sympy.core.power import Pow
from sympy.functions.elementary.piecewise import Piecewise
# from sympy.functions.combinatorial.factorials import factorial
from sympy import *
from sympy.functions import *
from symfit.core.objectives import LogLikelihood, LeastSquares
import matplotlib.pyplot as plt


if __name__ == "__main__":
    algo="bionet"
    dataset="ERS_1"
    output=pd.read_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR,"emp_fdr","{}_MAX".format(dataset),"emp_diff_{}_{}.tsv".format(dataset, algo)), sep='\t')
    output=output.loc[output["filtered_pval"].sort_values(ascending=False).index, :]
    output = output.dropna()
    output = output.loc[output['dist_n_samples'].values != '0',:]
    counter=0
    for index, cur in output.iterrows():
        if  counter > 5 or counter < 5:
            counter += 1
            continue

        pval=np.array([float(x) for x in cur["dist_n_samples"][1:-1].split(", ")])
        pval=pval[pval!=0]
        # pval = pval-np.min(pval)+1E-0

        if np.size(pval) == 0:
            pval=np.array([0])
        n_bins='sqrt' #int(np.max(pval))

        # # beta_params = scipy.stats.beta.fit(pval)
        # chi2_params = scipy.stats.chi2.fit(pval)
        ext_params = scipy.stats.genextreme.fit(pval)

        # print beta_params

        # model, fit_result_params, objective_value = mixture_dist.fit(pval)
        k = Parameter("k", values=0.06828227772281918)
        # sigma2 = Parameter("sigma2", value=1)
        # mu = Parameter("mu", value=np.max(data) - 1)
        # V = Parameter("V", value=v)
        A = Parameter("A", value=6.6130080008380685, fixed=True)
        B = Parameter("B", value=11.067416534379872)
        x = Variable()
        # model = Piecewise((Mul(Mul(exp(-Pow((1 - Mul(k, Mul(Add(x, -B), 1 / A))), 1 / k)),
        #                        Pow((1 - Mul(k, Mul(Add(x, -B), 1 / A))), (1 / k - 1))),1/A), GreaterThan(1 / k, Mul(Add(x, -B), 1 / A))),
        #                   (1e-09, True))
        model, fit_result_params, objective_value=mixture_dist.fit(pval)

        print ext_params
        print fit_result_params

        # ls=np.linspace(0, 100, 1000)
        # plt.plot(ls, scipy.stats.chi2(*chi2_params).pdf(ls))
        # print chi2_params
        # # plt.hist(data, normed=True, bins=50)
        # plt.plot(ls, model(ls,**fit_result.params) )
        # print fit_result.params
        # plt.title("args: {}".format(fit_result.params))
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

        ls=np.linspace(0, 100, 1000)
        plt.plot(ls, scipy.stats.genextreme(*ext_params).pdf(ls), color='blue') # sse_chi2
        plt.plot(ls, model(ls,**fit_result_params) ) # **{"k":ext_params[0], "A":ext_params[2], "B":ext_params[1]})) #   #

        # # plot_title="terms: {}\npval: {} (real: {})\nk: {},mu: {}, sigma2: {}, A: {}, B: {}"\
        # #     .format(cur["GO name"], fit_result.params['A'] * scipy.stats.norm(fit_result.params['mu'], fit_result.params['sigma2']).sf(cur["filtered_pval"]) + \
        # fit_result.params['B'] * scipy.stats.expon(fit_result.params['k']).sf(cur["filtered_pval"]),
        #              cur["filtered_pval"], fit_result.params['k'], fit_result.params['mu'], fit_result.params['sigma2'], fit_result.params['A'], fit_result.params['B'])
        plot_title = "terms: {}".format(cur["GO name"])
        # # print plot_title
        sns.distplot(pval, norm_hist=True, kde=False, bins=n_bins)
        # plt.title("GO pval dist.")
        real_n_term=np.size(pval)
        plt.legend()
        print cur["GO name"]
        plt.title(plot_title)
        plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "go_pval_dist_{}_{}_{}.png".format(dataset,algo, counter)))
        plt.clf()
        counter+=1
        # if counter >50:
        #     break

