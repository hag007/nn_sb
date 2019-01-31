import numpy as np

from symfit import Parameter, Variable, Fit, sqrt, pi, Equality, Abs, GreaterThan
from sympy.functions.elementary.exponential import exp
from sympy.functions.elementary.piecewise import Piecewise
from sympy.core.add import Add
from sympy.core.mul import Mul
from sympy.core.power import Pow
from sympy.functions.special.gamma_functions import gamma
from symfit.core.objectives import LogLikelihood


def fit(data, v=1.0, return_dict=None):




    k = Parameter("k", values=0.06704382317644017)
    # l = Parameter("l", values=1.0)
    # sigma2 = Parameter("sigma2", value=1)
    # mu = Parameter("mu", value=np.max(data) - 1)
    # V = Parameter("V", value=0.99, fixed=True)
    A = Parameter("A", value=0.0)
    B = Parameter("B", value=1.0)
    # C = Parameter("C", value=1.0)
    # D = Parameter("D", value=0.0)

    x = Variable()


    # model = Add(Mul(V, Piecewise((Mul((1 / (2 ** (k / float(2)) * gamma(k / 2)) * Pow(Mul(Add(x, -B), 1 / A),
    #                                                                                   (k / 2 - 1)) * exp(
    #     -Mul(Add(x, -B), 1 / A) / 2)), 1 / A), GreaterThan(Add(x, -B), 0)), (1e-09, True))),
    #             Mul(Add(1, -V), Piecewise(((1 / (sqrt(2 * pi * np.abs(sigma2)))) * exp(
    #                 -(x - mu) ** 2 / (2 * np.abs(sigma2))), GreaterThan(Add(x, -B), 0)), (1e-09, True))))
    # model = Add(Mul(V, Piecewise((Mul(Mul(1 / k , exp(-Mul(Mul(Add(x,-B),1/A) ,1/ k))),1/A), GreaterThan(Add(x,-B),0)),(1e-09  , True))),
    #         Mul(Add(1,-V), Piecewise(((1 / (sqrt(2 * pi * np.abs(sigma2)))) * exp(-(x - mu) ** 2 / (2 * np.abs(sigma2))),GreaterThan(Add(x,-B),0)), (1e-09  , True))))
    #
    # model = Add(Mul(V, Piecewise(
    #                 (exp(-(1 - k * Mul(Add(x, -B), 1 / A)) ** (1 / k)) * (1 - k * Mul(Add(x, -B), 1 / A)) ** (
    #                             1 / k - Mul(Add(x, -B), 1 / A)), GreaterThan(k, 0)),
    #                 (Mul(exp(-exp(Mul(Add(x,-B),1/A))), exp(Mul(Add(x,-B),1/A))) , True))),
    #
    #             Mul(V, Piecewise(
    #                 (exp(-(1 - k * Mul(Add(x, -B), 1 / A)) ** (1 / k)) * (1 - k * Mul(Add(x, -B), 1 / A)) ** (
    #                         1 / k - Mul(Add(x, -B), 1 / A)), GreaterThan(k, 0)),
    #                 (Mul(exp(-exp(Mul(Add(x, -B), 1 / A))), exp(Mul(Add(x, -B), 1 / A))), True)))
    #             )

    # model = Piecewise(
    #     ((exp(-(1 - k * Mul(Add(x, -B), 1 / A)) ** (1 / k)) , (1 - k * Mul(Add(x, -B), 1 / A)) ** (
    #             1 / k - Mul(Add(x, -B), 1 / A))), GreaterThan(k, 0)),
    #     (Mul(exp(-exp(Mul(Add(x, -B), 1 / A))), exp(Mul(Add(x, -B), 1 / A))), True))

    # model = Mul(exp(-exp(Mul(Add(x, -B), 1 / A))), exp(Mul(Add(x, -B), 1 / A)))

    # model = Add(Mul(V,Piecewise((Mul(Mul(exp(-Pow((1 - Mul(k, Mul(Add(x, -B), 1 / A))), 1 / k)),
    #                            Pow((1 - Mul(k, Mul(Add(x, -B), 1 / A))), (1 / k - 1))), 1 / A),
    #                    GreaterThan(1 / k, Mul(Add(x, -B), 1 / A))),
    #                   (1e-09, True))),
    #             Mul(1-V, Piecewise((Mul(Mul(exp(-Pow((1 - Mul(k, Mul(Add(x, -D), 1 / C))), 1 / l)),
    #                                       Pow((1 - Mul(l, Mul(Add(x, -D), 1 / C))), (1 / k - 1))), 1 / C),
    #                               GreaterThan(1 / k, Mul(Add(x, -D), 1 / C))),
    #                              (1e-09, True)))
    #             )

    # model = Piecewise((Mul(Mul(exp(-Pow((1 - Mul(k, Mul(Add(x, -B), 1 / A))), 1 / k)),
    #                                       Pow((1 - Mul(k, Mul(Add(x, -B), 1 / A))), (1 / k - 1))), 1 / A),
    #                               GreaterThan(1 / k, Mul(Add(x, -B), 1 / A))),
    #                              (1e-09, True))

    model = Piecewise((Mul(Mul(exp(-Pow((1 - Mul(k, Mul(Add(x, -B), 1 / A))), 1 / k)),
                           Pow((1 - Mul(k, Mul(Add(x, -B), 1 / A))), (1 / k - 1))),1/A), GreaterThan(1 / k, Mul(Add(x, -B), 1 / A))),
                      (1e-09, True))

    # Do the fitting!   ##
    fit = Fit(model, data, objective=LogLikelihood, constraints=[
        GreaterThan(A, 0.1),
        GreaterThan(20, A),
        GreaterThan(B, 0),
        GreaterThan(20, B),

        # GreaterThan(C, 0.01),
        # GreaterThan(20, C),
        # GreaterThan(D, -1),
        # GreaterThan(13, D),

        # GreaterThan(mu, 5),
        # GreaterThan(np.max(data), mu),
        # GreaterThan(sigma2, 0.5),
        # GreaterThan(5, sigma2),
        GreaterThan(k, -10.0),
        GreaterThan(10.0, k), #2
        # GreaterThan(l, 0.01),
        # GreaterThan(10.0, l),  # 2
        # GreaterThan(1.0, V),
        # GreaterThan(V, 0.01), # 0.00001
    ])
    fit_result = None
    try:
        fit_result = fit.execute(tol=1e-09)  # tol=1e-09
    except ValueError:
        print(ValueError)
        return None

    if return_dict is not None:
        return_dict.append((model, fit_result.params, fit_result.gof_qualifiers["objective_value"]))
    return model, fit_result.params, fit_result.gof_qualifiers["objective_value"]

