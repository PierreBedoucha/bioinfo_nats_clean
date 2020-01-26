"""This script plots the lines benchmarking the energy score (pyROSETTA API) per structure

    in the current directory and after their relaxation and per amplitudes of deformation.

    The RMSD is overlaid as continuous color scale for each datapoint. The energy score is normalized per structure residue.
"""

import pandas as pd
import seaborn as sns
import os
import matplotlib.pyplot as plt
import math
from scipy.spatial import distance
import numpy as np
import Bio.PDB
import Bio.AlignIO as al
from sklearn import preprocessing


def read_pfam_align():
    """Reads multiple sequence alignment profile from pfam.
    It determines the envelope (start and end) of the GNAT fold

    Reads values from pfam_env.txt file in ../data/input/etc

    :return: Dictionary. Keys: structure pdb id, Values:
    :rtype: dict
    """
    file_path = os.path.join("../data/input/etc", "pfam_env.txt")
    pdb_align_dict = {}
    with open(file_path) as f1:
        for line in f1:
            if not line.startswith("#") and not line.startswith("\n"):
                pdb_align_dict[line[0:4]] = (int(line[15:17]), int(line[21:24]))
    return pdb_align_dict


def read_pdb_starts():
    """Reads at which index each pdb sequence is starting from the pdb_starts.txt file from ../data/input/etc

    :return: Dictionary. Keys: structure pdb id, Values: starting index
    :rtype: dict
    """
    file_path = os.path.join("../data/input/etc", "pdb_starts.txt")
    pdb_starts_dict = {}
    with open(file_path) as f1:
        for line in f1:
            if not line.startswith("#") and not line.startswith("\n"):
                line_array = line.split(',')
                pdb_starts_dict[line[0:4]] = int(line_array[1])
    return pdb_starts_dict


def read_msa_fasta():
    """Reads multiple structure alignment from MUSTANG.
    It determines the structurally aligned core of the proteins.

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


def compute_rmsd(pdb_path1, pdb_path2, **kwargs):
    """Computes RMS distance between two pdb structures and only from start to end indices in their sequence

    :param pdb_path1: First pdb file name
    :type pdb_path1: str

    :param pdb_path2: First pdb file name
    :type pdb_path2: str

    :param kwargs: Keyword arguments with optional start and end in the pdb sequence instead of obtaining the
        boundaries from reading the msa data.

    :return: Rmsd value between pdb structures
    :rtype: float
    """
    sum_dist_sq = 0
    atom_cpt = 1
    file1_ref_array = pdb_path1.split('_')
    file2_ref_array = pdb_path2.split('_')
    if not os.path.exists(pdb_path1):
        file1_ref_array[-3] = "1"
        pdb_path1 = "_".join(file1_ref_array)
    if not os.path.exists(pdb_path2):
        file2_ref_array[-3] = "1"
        pdb_path2 = "_".join(file2_ref_array)
    with open(pdb_path1) as f1, open(pdb_path2) as f2:
        for line1, line2 in zip(f1, f2):
            if 'ATOM' in line1[0:6] and ' CA ' in line1[12:16]:
                if 'ATOM' in line2[0:6] and ' CA ' in line2[12:16]:
                    if 'start' in kwargs and 'end' in kwargs:
                        if (kwargs.get("start") <= int(line1[23:26].strip()) <= kwargs.get("end")) \
                                and (kwargs.get("start") <= int(line2[23:26].strip()) <= kwargs.get("end")):
                            distances = calc_distance(line1, line2)
                            sum_dist_sq += math.pow(distances[0][0], 2)
                            atom_cpt += 1
                    else:
                        pdb1_res_list = ca_align_dict[file1_ref_array[-8]]
                        pdb2_res_list = ca_align_dict[file2_ref_array[-8]]
                        if (int(line1[23:26].strip()) in pdb1_res_list) and \
                                (int(line2[23:26].strip()) in pdb2_res_list):
                            distances = calc_distance(line1, line2)
                            sum_dist_sq += math.pow(distances[0][0], 2)
                            atom_cpt += 1
    rmsd = math.sqrt(sum_dist_sq / atom_cpt)
    return rmsd


def calc_distance(pdbfile_line_1, pdbfile_line_2):
    """Calculate in line distance (Angstroms) between two atoms.

    :param pdbfile_line_1: Str line for atom line in first pdb file

    :param pdbfile_line_2: Str line for atom line in second pdb file

    :return: Distance array from distance.cdist method
    :rtype: float[][]
    """
    try:
        dist = distance.cdist(
            np.array([np.float64(val) for val in pdbfile_line_1[31:54].split()]).reshape(1, -1),
            np.array([np.float64(val) for val in pdbfile_line_2[31:54].split()]).reshape(1, -1))
    except ValueError as ex:
        dist = distance.cdist(np.array([np.float64(pdbfile_line_1[30:38]),
                                        np.float64(pdbfile_line_1[38:46]),
                                        np.float64(pdbfile_line_1[46:54])]).reshape(1, -1),
                              np.array([np.float64(pdbfile_line_2[30:38]),
                                        np.float64(pdbfile_line_2[38:46]),
                                        np.float64(pdbfile_line_2[46:54])]).reshape(1, -1))
    return dist


def f(x, y, z, **kwargs):
    """Add annotation on seaborn plot

    :param x: x value for the annotation position on plot

    :param y: y value for the annotation poistion on plot

    :param z: Annotated value

    :param kwargs: keyword arguments
    """
    ax = sns.pointplot(x, y, **kwargs)
    ax.axhline(5, alpha=0.5, color='grey')
    for i in range(len(x)):
        ax.annotate('{:6.2f}'.format(z.values[i]), xy=(i, z.values[i]), fontsize=8,
                    color=kwargs.get("color", "k"),
                    bbox=dict(pad=.9, alpha=1, fc='w', color='none'),
                    va='center', ha='center', weight='bold')


def facet_scatter(x, y, c, **kwargs):
    """Draw scatterplot with point colors from a faceted DataFrame columns.

    :param x: x axis data

    :param y: y axis data

    :param kwargs: keyword arguments
    """
    kwargs.pop("color")
    plt.scatter(x, y, c=c, **kwargs)


# GLOBAL VARIABLES
amplitude_max = 30
ca_align_list = []
ca_align_dict = read_msa_fasta()

if __name__ == '__main__':
    train = pd.read_csv('Workbook24.csv')

    dict_ref_SC = {}
    dict_ref_relax = {}

    for file in os.listdir("."):
        if file.endswith(".pdb"):
            if 'minimized' not in file:
                # start, end = read_pfam_align()[file[0:4]]
                file_ref_temp = file.replace(file.split("_")[-2], "0.00")
                file_ref = file_ref_temp.replace(file_ref_temp.split("_")[3], "a{0}.00".format(str(amplitude_max)))
                if file.split("_")[-2] != "0.00":
                    train.loc[train.pdb_filename == file + "_1", 'rmsd_init'] = compute_rmsd(file, file_ref)
                else:
                    train.loc[train.pdb_filename == file + "_1", 'rmsd_init'] = 0.00
            elif 'minimized' in file:
                # start, end = read_pfam_align()[file.split("_")[-8]]
                file_ref_temp = file.replace(file.split("_")[-2], "0.00")
                file_ref = file_ref_temp.replace(file_ref_temp.split("_")[-5], "a{0}.00".format(str(amplitude_max)))
                if file.split("_")[-2] != "0.00":
                    # Align Relaxed structures and use rmsd align (CA only and pfam sequences)
                    # select_ca_align(file, start, end)
                    res_to_be_aligned = ca_align_dict[file.split("_")[-8]]
                    pdb_parser = Bio.PDB.PDBParser(QUIET=True)
                    # Get the structures
                    try:
                        ref_structure = pdb_parser.get_structure("reference", file_ref)
                    except FileNotFoundError as err:
                        print("You chose the 0.00 file with tag 1 instead of 2. Retrying...")
                        file_ref_array = file_ref.split('_')
                        file_ref_array[-3] = "1"
                        file_ref = "_".join(file_ref_array)
                        ref_structure = pdb_parser.get_structure("reference", file_ref)
                    sample_structure = pdb_parser.get_structure("sample", file)
                    ref_model = ref_structure[0]
                    sample_model = sample_structure[0]
                    ref_atoms = []
                    sample_atoms = []
                    # Iterate of all chains in the model in order to find all residues
                    for ref_chain in ref_model:
                        # Iterate of all residues in each model in order to find proper atoms
                        for ref_res in ref_chain:
                            # Check if residue number ( .get_id() ) is in the list
                            if ref_res.get_id()[1] in res_to_be_aligned:
                                # Append CA atom to list
                                ref_atoms.append(ref_res['CA'])
                    # Do the same for the sample structure
                    for sample_chain in sample_model:
                        for sample_res in sample_chain:
                            if sample_res.get_id()[1] in res_to_be_aligned:
                                sample_atoms.append(sample_res['CA'])
                    super_imposer = Bio.PDB.Superimposer()
                    super_imposer.set_atoms(ref_atoms, sample_atoms)
                    super_imposer.apply(sample_model.get_atoms())

                    file = '_'.join(file.split('_')[-8:])
                    train.loc[train.pdb_filename == file + "_1", 'rmsd_relax'] = super_imposer.rms
                else:
                    file_array = file.split('_')
                    file = "_".join(file_array[-8:])
                    train.loc[train.pdb_filename == file + "_1", 'rmsd_relax'] = 0.00

    train.to_csv("test_out.csv", sep=';', encoding='utf-8')

    train['pdbid'] = train.pdb_filename.str[:4]
    # new data frame with split value columns
    new = train["pdb_filename"].str.split("_", n=7, expand=True)
    train['amplitude'] = new[6]
    train['amplitude'] = train['amplitude'].astype(float)

    column_names_to_normalize = ['score_init']

    min_max_scaler = preprocessing.MinMaxScaler()
    x = train[column_names_to_normalize].values
    x_scaled = min_max_scaler.fit_transform(x)
    df_temp = pd.DataFrame(x_scaled, columns=column_names_to_normalize, index=train.index)
    train[column_names_to_normalize] = df_temp

    col_name = ['pdb_filename', 'rmsd_init', 'score_init', 'rmsd_relax', 'score_relax', 'pdbid', 'amplitude']

    train = train.loc[(train['pdbid'] != "5isv")]
    train = train.loc[(train['pdbid'] != "2cns")]

    grouped = train.groupby(["pdbid"])
    train = grouped.apply(lambda x: x.sort_values(["amplitude"], ascending=True)).reset_index(drop=True)

    h = sns.FacetGrid(train, col="pdbid", palette='seismic', sharey=False, sharex=True, col_wrap=6, height=2, aspect=1)

    vmin = train['rmsd_relax'].min()
    vmax = train['rmsd_relax'].max()

    cmap = sns.light_palette("seagreen", as_cmap=True)

    h.map(facet_scatter, "amplitude", "score_relax", "rmsd_relax", s=100, vmin=vmin, vmax=vmax, cmap=cmap)

    # Make space for the colorbar
    # h.fig.subplots_adjust(right=.92)
    # plt.tight_layout()

    cax = h.fig.add_axes([.40, .0339, .2, .023])

    # Get a mappable object with the same colormap as the data
    points = plt.scatter([], [], c=[], vmin=vmin, vmax=vmax, cmap=cmap)
    h.map(plt.plot, "amplitude", "score_relax")

    # Draw the colorbar
    cbar = h.fig.colorbar(points, cax=cax, orientation='horizontal')
    cbar.ax.set_title('RMSD relax', fontsize=10)

    plt.show()
