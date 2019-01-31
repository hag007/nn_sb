#!/usr/bin/env python
"""Calculate differentially expressed genes using EdgeR from bioconductor.
http://bioconductor.org/packages/2.5/bioc/html/edgeR.html
Usage:
    count_diffexp.py <count_file>
"""
import os
import pandas as pd

import numpy as np
# import rpy2.robjects.numpy2ri  as numpy2ri
# numpy2ri.activate()

from rpy2.robjects import pandas2ri
pandas2ri.activate()


import constants
import infra
import shutil
import time
import random
import utils.add_t_test_to_ge as t_test
from utils.r_runner import run_rscript


def run_FDR(pvals):

    script = file(os.path.join(constants.REPO_DIR,"r","scripts","robust_fdr.r")).read()
    return run_rscript(script=script, pvals=pvals)

