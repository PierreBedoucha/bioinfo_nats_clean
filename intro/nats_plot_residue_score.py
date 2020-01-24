import os
import re
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns
import Bio.AlignIO as al

"""
This script plot the energy score (from pyROSETTA API) over all the structure files and per aligned residues.
The structures scores are listed in the relaxed pdb files and the script will isolate it for all the residues.
The final results is summed through all the current directory structure files.
"""


def read_msa_fasta():
    """
    Reads multiple structure alignment from MUSTANG. It determines the structurally aligned core of the proteins.
    Note: here, only the aligned regions are of interest, gaps are removed.
    :return: Dictionary. Keys: structure pdb id, Values: aligned indices
    :rtype: dict
    """
    pdb_align_dict = {'3tfy': [], '5isv': [], '4pv6': [], '2z0z': [], '1s7l': [], '2x7b': [], '3igr': [], '5k18': [],
                      '2cns': [],
                      '5hh0': [], '5wjd': [], '5icv': [], '4kvm': [], '4u9v': [], }
    file_path = os.path.join("../data/input/etc", "nats_alignment.afasta")
    records = al.read(open(file_path), "fasta")
    tlist = list(zip(*records))
    for i in range(0, records.get_alignment_length()):
        if '-' not in [y for y in tlist][i]:
            for rec in records:
                if not rec.id[0:4] == '4ua3':
                    ls = [i for i, e in enumerate(rec.seq) if e != '-']
                    res_cpt = ls.index(i)
                    pdb_align_dict[rec.id[0:4]].append(res_cpt + read_pdb_starts()[rec.id[0:4]])
    return pdb_align_dict


def read_pdb_starts():
    """
    Reads at which index each pdb sequence is starting from the pdb_starts.txt file from ../data/input/etc
    :return: Dictionary. Keys: structure pdb id, Values: starting index
    :rtype: dict
    """
    file_path = os.path.join("../data/input/etc", "pdb_starts.txt")
    pdb_starts_dict = {}
    with open(file_path) as f1_start:
        for line in f1_start:
            if not line.startswith("#") and not line.startswith("\n"):
                line_array = line.split(',')
                pdb_starts_dict[line[0:4]] = int(line_array[1])
    return pdb_starts_dict


ca_align_dict = read_msa_fasta()

if __name__ == '__main__':

    res_cpt_list = []
    res_scores_list = []
    df = pd.DataFrame(columns=['res_cpt', 'res_scores'])

    for file in os.listdir("."):
        if 'minimized' in file:
            res_to_be_aligned = ca_align_dict[file.split("_")[-8]]
            search = re.compile(r'[^A-Z]*_[0-9]+\s').search
            res_count = 1
            with open(file) as f1:
                for mini_line in f1:
                    if search(mini_line):
                        mini_line_array = mini_line.split(' ')
                        res_cpt_list.append(res_count)
                        res_scores_list.append(float(mini_line_array[-1].rstrip()))
                        df = df.append({'res_cpt': int(res_count), 'res_scores': float(mini_line_array[-1].rstrip())},
                                       ignore_index=True)
                        df['res_cpt'] = df['res_cpt'].astype(int)

                        res_count += 1

    ax = sns.boxplot(x="res_cpt", y="res_scores", data=df)
    ax.get_xaxis().set_major_formatter(
        ticker.FuncFormatter(lambda x, p: format(int(x))))
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45)
    ind_list = list(range(0, 200, 10))
    for ind, label in enumerate(ax.get_xticklabels()):
        if int(label._text) in ind_list:
            label.set_visible(True)
        else:
            label.set_visible(False)

    plt.show()
