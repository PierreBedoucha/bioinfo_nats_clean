import seaborn as sns
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

'''
The script analyses the normalized squared fluctuation data for the first 6 normal modes (successively). It plots the
fluctuation graphs for the 7 protein structures (see data), adds the mean value for baseline description, and highlights
the 2 loop regions. The modes of interest are the ones with the highest fluctuation values for the 2 functional loops. 
'''


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
        points = points[:, None]
    median = np.median(points, axis=0)
    diff = np.sum((points - median) ** 2, axis=-1)
    diff = np.sqrt(diff)
    med_abs_deviation = np.median(diff)

    modified_z_score = 0.6745 * diff / med_abs_deviation

    return modified_z_score > thresh


def read_pdb_starts():
    """
    Reads at which index each pdb sequence is starting from the pdb_starts.txt file from ../data/input/etc
    :return: Dictionary. Keys: structure pdb id (str), Values: starting index (ind)
    :rtype: dict
    """
    file_path = os.path.join("../data/input/etc", "pdb_starts.txt")
    pdb_starts_dict = {}
    with open(file_path) as pdb_starts:
        for pdb_start_line in pdb_starts:
            if not pdb_start_line.startswith("#") and not pdb_start_line.startswith("\n"):
                line_array = pdb_start_line.split(',')
                pdb_starts_dict[pdb_start_line[0:4]] = int(line_array[1])
    return pdb_starts_dict


def read_pdb_loops():
    """
    Reading loop locations for the different 7 protein structures.
    Reads pdb_loops.txt file from ../data/input/etc containing simple data obtained from pdb structures.
    :return: Dictionary. Keys: structure pdb id (str), Values: loop 1 start (int), loop 1 end (int),
    loop 2 start (int), loop 2 end (int)
    :rtype: dict
    """
    file_path = os.path.join("../data/input/etc", "pdb_loops.txt")
    pdb_loops_dict = {}
    with open(file_path) as pdb_loops:
        for pdb_loop_line in pdb_loops:
            if not pdb_loop_line.startswith("#") and not pdb_loop_line.startswith("\n"):
                line_array = pdb_loop_line.split(',')
                pdb_loops_dict[pdb_loop_line[0:4]] = [int(line_array[1].split(":")[1].split("-")[0]),
                                                      int(line_array[1].split(":")[1].split("-")[1]),
                                                      int(line_array[2].split(":")[1].split("-")[0]),
                                                      int(line_array[2].split(":")[1].split("-")[1])]
    return pdb_loops_dict


def facet_span(pdbid, y):
    """
    Facet span function to apply to the matplotlib facetGrid
    Plots the 2 loop areas and the median dotted-line in gray
    :param pdbid: pdbids column used as index in the dataframe
    :type pdbid: str
    :param y: fluctuation values as column in the dataframe
    :type y: float
    """
    plt.axvspan(read_pdb_loops()[pdbid.iloc[0]][0], read_pdb_loops()[pdbid.iloc[0]][1], facecolor="gray", alpha=0.2)
    plt.axvspan(read_pdb_loops()[pdbid.iloc[0]][2], read_pdb_loops()[pdbid.iloc[0]][3], facecolor="gray", alpha=0.2)
    plt.axhline(np.mean(y), linestyle="--", color="gray")


# Global variables (Ugly)
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
        for x in os.listdir(root):
            if os.path.isfile(os.path.join(root, x)) and not x.startswith("."):
                with open(os.path.join(root, x), "r") as f1:
                    for line in f1.readlines():
                        resid_list.append(
                            int(line.split(" ")[0].strip()) + (read_pdb_starts()[root.split("/")[-1]] - 1))
                        fluct_list.append(float(line.split(" ")[2].strip()))
                        mode_list.append(int(''.join([s for s in x.split("_")[2] if s.isdigit()])))
                        pdbid_list.append(root.split("/")[-1])

    train['resid'] = resid_list
    train['fluct_score'] = fluct_list
    train['mode'] = mode_list
    mask = is_outlier(np.array(fluct_list))
    train['is_outliers'] = mask
    train['pdbid'] = pdbid_list

    # In case the grooups over pdbids are need to evaluate features across common columns
    # grouped = train.groupby(["pdbid"])
    # train = grouped.apply(lambda x: x.sort_values(["resid"], ascending = True)).reset_index(drop=True)

    modes = [7, 8, 9, 10, 11, 12]
    for m in modes:
        train_mode = train.loc[(train['mode'] == m)]

        h = sns.FacetGrid(train_mode, col="pdbid", palette='seismic', sharey=False, sharex=False, col_wrap=4,
                          height=2, aspect=1)

        # In case the outliers are needed to detect local maximum values, the hue can be used
        # h.map(sns.scatterplot, x='resid', y='fluct_score', hue='is_outliers', data=train_mode)
        h.map(plt.plot, 'resid', 'fluct_score')
        h.map(facet_span, "pdbid", "fluct_score")
        plt.show()
