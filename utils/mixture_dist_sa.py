import numpy as np
import sys
from symfit import Parameter, Variable, Fit, sqrt, pi, Equality, Abs, GreaterThan
from sympy.functions.elementary.exponential import exp
from sympy.functions.elementary.piecewise import Piecewise
from sympy.core.add import Add
from sympy.core.mul import Mul
from sympy.core.power import Pow
from sympy.functions.special.gamma_functions import gamma
from symfit.core.objectives import LogLikelihood
import argparse



if __name__ == '__main__':


    parser = argparse.ArgumentParser(description='args')
    parser.add_argument('--datasets', dest='datasets', default="SOC")
    parser.add_argument('--prefix', dest='prefix', default="GE")
    parser.add_argument('--algos', dest='algos', default="jactivemodules_greedy")
    parser.add_argument('--network', dest='network', default="dip")
    parser.add_argument('--n_start', help="number of iterations (total n permutation is pf*(n_end-n_start))", dest='n_start', default=0)
    parser.add_argument('--n_end', help="number of iterations (total n permutation is pf*(n_end-n_start))", dest='n_end', default=1000)
    parser.add_argument('--data', help="data", dest='data', default="0,1")
    parser.add_argument('--v', help="v", dest='v', default="1.0")


    args = parser.parse_args()


    data=np.array(args.data.split(",")).astype(np.float)
    v=float(args.v)
    return_dict=None


    # k = Parameter("k", values=1.0, max=100, min=0.01)
    # sigma2 = Parameter("sigma2", value=1, min=0.5, max=5.0)
    # mu = Parameter("mu", value=np.max(data) - 1, min=5.0, max=np.max(data))
    # V = Parameter("V", value=v, min=0.0001, max=1)
    # A = Parameter("A", value=1.0, min=0.2, max=10)
    # B = Parameter("B", value=0.0, min=0.0, max=min(data)+1)

    # k = Parameter("k", values=1.0, max=100, min=0.01)
    sigma2 = Parameter("sigma2", value=1, min=0.5, max=5.0)
    mu = Parameter("mu", value=max(np.max(data) - 1, 5.0), min=5.0, max=max(np.max(data), 5.0))
    V = Parameter("V", value=v, min=0.0001, max=1)
    A = Parameter("A", value=1.0, min=0.2, max=10)
    B = Parameter("B", value=0.0, min=0.0, max=max(min(data) + 1, 1))

    x = Variable()


    # model = Add(Mul(V, Piecewise((Mul((1 / (2 ** (k / float(2)) * gamma(k / 2)) * Pow(Mul(Add(x, -B), 1 / A),
    #                                                                                   (k / 2 - 1)) * exp(
    #     -Mul(Add(x, -B), 1 / A) / 2)), 1 / A), GreaterThan(Add(x, -B), 0)), (1e-09, True))),
    #             Mul(Add(1, -V), Piecewise(((1 / (sqrt(2 * pi * np.abs(sigma2)))) * exp(
    #                 -(x - mu) ** 2 / (2 * np.abs(sigma2))), GreaterThan(Add(x, -B), 0)), (1e-09, True))))
    # model = Add(Mul(V, Piecewise((Mul(Mul(1 / k , exp(-Mul(Mul(Add(x,-B),1/A) ,1/ k))),1/A), GreaterThan(Add(x,-B),0)),(1e-09  , True))),
    #         Mul(Add(1,-V), Piecewise(((1 / (sqrt(2 * pi * np.abs(sigma2)))) * exp(-(x - mu) ** 2 / (2 * np.abs(sigma2))),GreaterThan(Add(x,-B),0)), (1e-09  , True))))
    model = Add(Mul(V, Piecewise(
        (Mul(exp(-Mul(Add(x, -B), 1 / A)), 1 / A), GreaterThan(Add(x, -B), 0)), (1e-09, True))),
                Mul(Add(1, -V), Piecewise(((1 / (sqrt(2 * pi * np.abs(sigma2)))) * exp(
                    -(x - mu) ** 2 / (2 * np.abs(sigma2))), GreaterThan(Add(x, -B), 0)), (1e-09, True))))

    # Do the fitting!   ##
    fit = Fit(model, data, objective=LogLikelihood, constraints=[
        GreaterThan(A, 0.2),
        GreaterThan(10, A),
        GreaterThan(B, 0),
        GreaterThan(max(min(data) + 1, 1), B),
        GreaterThan(mu, 5),
        GreaterThan(max(np.max(data), 5.0), mu),
        GreaterThan(sigma2, 0.5),
        GreaterThan(5, sigma2),
        # GreaterThan(k, 0.01),
        # GreaterThan(100, k), #2
        GreaterThan(1.0, V),
        GreaterThan(V, 0.0001), # 0.00001
    ])
    fit_result = None
    try:
        fit_result = fit.execute(tol=1e-09)  # tol=1e-09
    except:
        print "None"
        sys.stdout.write("None")
        exit(0)

    if fit_result is not None:
        # return_dict.append((model, fit_result.params, fit_result.gof_qualifiers["objective_value"]))
        sys.stdout.write(",".join(["{}:{}".format(k,v) for k, v in fit_result.params.iteritems()]) +"#"+str(fit_result.gof_qualifiers["objective_value"]))
