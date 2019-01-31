import numpy as np
import scipy
from scipy.stats import hypergeom
import logging
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('-n', '--N', type=int, default=8490)
parser.add_argument('-b', '--B', type=int, default=13)

args = parser.parse_args()
N = args.N
B = args.B
print "about to start HGT for N={} and B={}".format(N,B)

logging.basicConfig(filename='hgtlog.log',level=logging.INFO)
logger = logging.getLogger("log")
logger.info("about to start HGT for N={} and B={}\n".format(N,B))
HGTs = None


for n in range(N+1):
     print "n1: {}".format(n)
     #logger.info("n1: {}\n".format(n))
     b_tails = hypergeom.sf(np.arange(0, B + 1), N, B, n)
     # b_tails = np.add(hypergeom.sf(np.arange(0, B + 1), N, B, n), hypergeom.pmf(np.arange(0, B + 1), N, B, n))
     if type(HGTs)!=type(None):
         HGTs = np.c_[HGTs, b_tails]
     else:
         HGTs = b_tails
HGTs[0][0] = 1.0
np.save(os.path.join("HGTs_out_{}_{}".format(N,B)), HGTs)
