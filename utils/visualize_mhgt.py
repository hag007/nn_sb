import numpy as np
import scipy
from scipy.stats import hypergeom
import logging
import argparse
import os
import matplotlib.pyplot as plt
from matplotlib import colors


BASE_OUTPUT_DIR = "c:\\users\hagai\\desktop\\"

def visualize(hgt_preprocessing_file_name):
    HGTs =  np.load(os.path.join(BASE_OUTPUT_DIR, hgt_preprocessing_file_name))
    HGTs[HGTs == 0] = -1
    B,N = np.shape(HGTs)
    B -=1
    N -= 1
    mHGT = 0.0002
    left_tails = [(hypergeom.sf(i, N, B, i) + hypergeom.pmf(i, N, B, i)) for i in range(0, B + 1)] #
    top_tails = [(hypergeom.sf(B, N, B, i) + hypergeom.pmf(B, N, B, i)) for i in range(0, N + 1)] #

    for i, cur in enumerate(left_tails):
        if cur <= mHGT:
            left_edge = (i,i)
            break

    for i, cur in enumerate(top_tails):
        if cur >= mHGT:
            top_edge = (i-1,B)
            break
    slope = (left_edge[1] - top_edge[1]) / ((left_edge[0] - top_edge[0])*1.0)
    constant = top_edge[1] - slope*top_edge[0]

    cmap = colors.ListedColormap(["red", 'gray', 'skyblue', "red"])
    bounds = [-1, -0.1, mHGT, 1, 1]
    norm = colors.BoundaryNorm(bounds, cmap.N)

    fig, ax = plt.subplots()
    ax.imshow(HGTs, cmap=cmap, norm=norm)

    # draw gridlines
    ax.grid(which='minor', axis='both', linestyle='-', color='k', linewidth=1)
    ax.invert_yaxis()
    # ax.set_xticks(np.arange(-.5, 10, 1));
    ax.set_yticks(np.arange(0, 101, 8000));
    plt.plot([left_edge[0], top_edge[0]], [left_edge[1], top_edge[1]], "green")
    plt.show()
    x=1

visualize("HGTs_out_1000_100.npy")