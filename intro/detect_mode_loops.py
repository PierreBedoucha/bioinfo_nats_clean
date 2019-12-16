import seaborn as sns
from tempfile import mkstemp
from shutil import move
import os
import matplotlib.pyplot as plt
import math
import numpy as np
import pandas as pd

def is_outlier(points, thresh=3.5):
    """
    Returns a boolean array with True if points are outliers and False
    otherwise.

    Parameters:
    -----------
        points : An numobservations by numdimensions array of observations
        thresh : The modified z-score to use as a threshold. Observations with
            a modified z-score (based on the median absolute deviation) greater
            than this value will be classified as outliers.

    Returns:
    --------
        mask : A numobservations-length boolean array.

    References:
    ----------
        Boris Iglewicz and David Hoaglin (1993), "Volume 16: How to Detect and
        Handle Outliers", The ASQC Basic References in Quality Control:
        Statistical Techniques, Edward F. Mykytka, Ph.D., Editor.
    """
    if len(points.shape) == 1:
        points = points[:,None]
    median = np.median(points, axis=0)
    diff = np.sum((points - median)**2, axis=-1)
    diff = np.sqrt(diff)
    med_abs_deviation = np.median(diff)

    modified_z_score = 0.6745 * diff / med_abs_deviation

    return modified_z_score > thresh

def read_pdb_starts():
    file_path = os.path.join("../data/input/etc", "pdb_starts.txt")
    pdb_starts_dict = {}
    with open(file_path) as f1:
        for line in f1:
            if not line.startswith("#") and not line.startswith("\n"):
                line_array = line.split(',')
                pdb_starts_dict[line[0:4]] = int(line_array[1])
    return pdb_starts_dict

def read_pdb_loops():
    file_path = os.path.join("../data/input/etc", "pdb_loops.txt")
    pdb_loops_dict = {}
    with open(file_path) as f1:
        for line in f1:
            if not line.startswith("#") and not line.startswith("\n"):
                line_array = line.split(',')
                pdb_loops_dict[line[0:4]] = [int(line_array[1].split(":")[1].split("-")[0]),
                                              int(line_array[1].split(":")[1].split("-")[1]),
                                              int(line_array[2].split(":")[1].split("-")[0]),
                                              int(line_array[2].split(":")[1].split("-")[1])]
    return pdb_loops_dict


def facet_span(pdbid, y, **kwargs):
    plt.axvspan(read_pdb_loops()[pdbid.iloc[0]][0], read_pdb_loops()[pdbid.iloc[0]][1], facecolor="gray", alpha=0.2)
    plt.axvspan(read_pdb_loops()[pdbid.iloc[0]][2], read_pdb_loops()[pdbid.iloc[0]][3], facecolor="gray", alpha=0.2)
    plt.axhline(np.mean(y), linestyle="--", color="gray")


pdb_file_list = ["4kvm", "5k18", "4u9v", "3tfy", "5hh0", "5icv", "5wjd"]


if __name__ == '__main__':
    directory_list = []
    dir_name_list = []


    for root, dirs, files in os.walk("../data/input/etc", topdown=False):
        for name in dirs:
            if name in pdb_file_list:
                directory_list.append(os.path.join(root, name))
                dir_name_list.append(name)

    train = pd.DataFrame(columns=['pdbid', 'resid', 'mode', 'fluct_score', 'is_outliers', 'is_loops'])
    resid_list = []
    fluct_list = []
    mode_list = []
    pdbid_list = []
    for root, name in zip(directory_list, dir_name_list):
        # files_list.extend([x for x in os.listdir(dir) if os.path.isfile(x)])
        for x in os.listdir(root):
            if os.path.isfile(os.path.join(root, x)) and not x.startswith("."):
                with open(os.path.join(root, x), "r") as f1:
                    for line in f1.readlines():
                        resid_list.append(int(line.split(" ")[0].strip()) + (read_pdb_starts()[root.split("/")[-1]]-1))
                        fluct_list.append(float(line.split(" ")[2].strip()))
                        mode_list.append(int(''.join([s for s in x.split("_")[2] if s.isdigit()])))
                        pdbid_list.append(root.split("/")[-1])

    train['resid'] = resid_list
    train['fluct_score'] = fluct_list
    train['mode'] = mode_list
    mask = is_outlier(np.array(fluct_list))
    train['is_outliers'] = mask
    train['pdbid'] = pdbid_list

    # grouped = train.groupby(["pdbid"])
    # train = grouped.apply(lambda x: x.sort_values(["resid"], ascending = True)).reset_index(drop=True)

    modes=[7,8,9,10,11,12]
    # modes = [7]
    for m in modes:
        train_mode = train.loc[(train['mode'] == m)]

        h = sns.FacetGrid(train_mode, col="pdbid", palette='seismic', sharey=False, sharex=False, col_wrap=4,
                          height=2, aspect=1)

        # h.map(sns.scatterplot, x='resid', y='fluct_score', hue='is_outliers', data=train_mode)
        h.map(plt.plot, 'resid', 'fluct_score')
        h.map(facet_span, "pdbid", "fluct_score")
        plt.show()

