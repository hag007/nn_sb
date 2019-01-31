
import rpy2.robjects as robjects
# import rpy2.robjects.numpy2ri  as numpy2ri
# numpy2ri.activate()

from rpy2.robjects import r, pandas2ri

pandas2ri.activate()


def run_rscript(script, output_vars = ["result"], **kwargs):
    """Call edgeR in R and organize the resulting differential expressed genes."""

    for k, v in kwargs.iteritems():
        robjects.globalenv[k.replace("_",".")] = v

    results = {}
    robjects.r(script)
    for cur_var in output_vars:
        results[cur_var] = robjects.r[cur_var.replace("_",".")]

    return results


