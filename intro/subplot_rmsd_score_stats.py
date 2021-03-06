"""This script handles PROCHECK analyses for the structure files in current directory.

    Be careful, reference pdb files (amplitude 0.00 for both native and relaxed structures) are needed for comparison.

    The PROCHECK files are created and placed in the current directory for later use.
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
import glob
import re


def read_pdb_resolution():
    """Reads the resolution from data for the pdb files.

    :return: Dictionary. Keys: structure pdb id, Values: resolution (float)
    :rtype: dict
    """
    file_path = os.path.join("../data/input/etc", "pdb_starts.txt")
    pdb_res_dict = {}
    with open(file_path) as f1:
        for line in f1:
            if not line.startswith("#") and not line.startswith("\n"):
                line_array = line.split(',')
                pdb_res_dict[line[0:4]] = int(line_array[1])
    return pdb_res_dict


def read_pfam_align():
    """Reads multiple sequence alignment profile from pfam.

    It determines the envelope (start and end) of the GNAT fold

    Reads values from pfam_env.txt file in ../data/input/etc

    :return: Dictionary. Keys: structure pdb id, Values: envelope indices (tuple)
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

    param pdbfile_line_1: Str line for atom line in first pdb file

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


def facet_stairs(x, y, **kwargs):
    """Draw scatterplot with point colors from a faceted DataFrame columns.

    :param x: x axis data

    :param y: y axis data

    :param kwargs: keyword arguments
    """
    plt.scatter(x, y, **kwargs)


def facet_line():
    """Draw scatterplot with point colors from a faceted DataFrame columns.
    """
    plt.axhline(-2.0, linestyle="--",
                color='gray')
    t = plt.text(3, -2.0, -2.0, horizontalalignment='right',
                 verticalalignment='center', color='gray')


def procheck(sc_pdbpath, pdb_resolution):
    """Main PROCHECK computation function. The full atom (sc) structure pb file is required as well as its corresponding
    resolution (Angstroms). All the resulting files are created in the current directory for later use.

    :param sc_pdbpath: Full atom pdb file path
    :type sc_pdbpath: str

    :param pdb_resolution: Resolution of the structure (Angstroms)
    :type pdb_resolution: float

    :return: Dictionary. Keys: PROCHECK measurement types, Values: measures
    :rtype: dict
    """
    import subprocess
    os.environ['prodir'] = '~/Software/procheck/procheck'
    if sc_pdbpath[:-4] + ".sum" not in os.listdir("."):
        p = subprocess.run(
            "{0} {1} {2}".format('~/Software/procheck/procheck/procheck.scr', sc_pdbpath, pdb_resolution),
            shell=True,
            stderr=subprocess.PIPE)
        procheck_error = p.stderr.decode('utf-8')

    # Get the value from the corresponding .sum file
    #  | Ramachandran plot: 91.5 % core 8.5 % allow 0.0 % gener 0.0 % disall |
    chars = []
    output_dict = {}
    with open("./" + sc_pdbpath[:-4] + ".sum") as f1:
        for line in f1:
            if not line.startswith("#") and not line.startswith("\n"):
                if "disall" in line:
                    line_array = line.split(' ')
                    chars.extend(line_array[-3])
                    output_dict['disall'] = float(''.join(chars[:-1]))
                if 'Bad contacts' in line:
                    line_array = line.split(' ')
                    output_dict['bad_contacts'] = float(line_array[-2])
                if 'Bond len/angle' in line:
                    line_array = re.findall(r"[-+]?\d*\.\d+|\d+", line)
                    output_dict['bond_lenangle'] = float(line_array[0])
                if 'G-factors' in line:
                    line_array = line.split(' ')
                    output_dict['g_factors'] = float(line_array[-2])
                if 'M/c bond lengths' in line:
                    line_array = re.findall(r"[-+]?\d*\.\d+|\d+", line)
                    output_dict['bond_lengths_highlighted'] = float(line_array[1])
                    if 'off graph' in line:
                        output_dict['bond_lengths_off'] = float(line_array[2])
                if 'M/c bond angles' in line:
                    line_array = re.findall(r"[-+]?\d*\.\d+|\d+", line)
                    output_dict['bond_angles_highlighted'] = float(line_array[1])
                    if 'off graph' in line:
                        output_dict['bond_angles_off'] = float(line_array[2])
    return output_dict


# GLOBAL VARIABLES
amplitude_max = 0
ca_align_list = []
ca_align_dict = read_msa_fasta()

if __name__ == '__main__':
    train = pd.read_csv('Workbook24.csv').dropna()

    dict_ref_SC = {}
    dict_ref_relax = {}

    for file in os.listdir("."):
        if file.endswith(".pdb"):
            if 'minimized' not in file:
                # start, end = read_pfam_align()[file[0:4]]
                if file.endswith("mini.pdb"):
                    file_ref_temp = file.replace(file.split("_")[-3], "0.00")
                    amplitude = file.split("_")[-3]
                else:
                    file_ref_temp = file.replace(file.split("_")[-2], "0.00")
                    amplitude = file.split("_")[-2]
                file_ref = file_ref_temp.replace(file_ref_temp.split("_")[3], "a{0}.00".format(str(amplitude_max)))
                temp_filename_ls = file_ref.split("_")[:2]
                temp_filename_ls.append("*")
                temp_filename_ls.extend(file_ref.split("_")[3:])
                temp_filename_ls[5] = "0"
                filelist = [file for file in glob.glob("_".join(temp_filename_ls))]
                file_ref = filelist[0]
                if amplitude != "0.00":
                    if file.endswith("mini.pdb"):
                        train.loc[train.pdb_filename == "1_" + file, 'rmsd_init'] = compute_rmsd(file, file_ref)
                    else:
                        train.loc[train.pdb_filename == "1_" + file, 'rmsd_init'] = compute_rmsd(file, file_ref)
                else:
                    if file.endswith("mini.pdb"):
                        train.loc[
                            train.pdb_filename == "1_" + file, 'rmsd_init'] = 0.00
                    else:
                        train.loc[train.pdb_filename == "1_" + file, 'rmsd_init'] = 0.00
            elif 'minimized' in file:
                if file.endswith("mini.pdb"):
                    file_ref_temp = file.replace(file.split("_")[-3], "0.00")
                    amplitude = file.split("_")[-3]
                    file_ref = file_ref_temp.replace(file_ref_temp.split("_")[-6], "a{0}.00".format(str(amplitude_max)))
                    pdbid_count = -9
                    snap_count = -4
                else:
                    file_ref_temp = file.replace(file.split("_")[-2], "0.00")
                    amplitude = file.split("_")[-2]
                    file_ref = file_ref_temp.replace(file_ref_temp.split("_")[-5], "a{0}.00".format(str(amplitude_max)))
                    pdbid_count = -8
                    snap_count = -3
                # start, end = read_pfam_align()[file.split("_")[-8]]
                temp_filename_ls = file_ref.split("_")[:5]
                temp_filename_ls.append("*")
                temp_filename_ls.extend(file_ref.split("_")[7:])
                temp_filename_ls[8] = "0"
                filelist = [file for file in glob.glob("_".join(temp_filename_ls))]
                file_ref = filelist[0]
                if amplitude != "0.00":
                    # Align Relaxed structures and use rmsd align (CA only and pfam sequences)
                    # select_ca_align(file, start, end)
                    res_to_be_aligned = ca_align_dict[file.split("_")[pdbid_count]]
                    pdb_parser = Bio.PDB.PDBParser(QUIET=True)
                    # Get the structures
                    try:
                        ref_structure = pdb_parser.get_structure("reference", file_ref)
                    except FileNotFoundError as err:
                        print("You chose the 0.00 file with tag 1 instead of 2. Retrying...")
                        file_ref_array = file_ref.split('_')
                        file_ref_array[snap_count] = "1"
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

                    file_base = '_'.join(file.split('_')[pdbid_count:])
                    if file_base.endswith("mini.pdb"):
                        train.loc[
                            train.pdb_filename == "1_" + file_base, 'rmsd_relax'] = super_imposer.rms
                    else:
                        train.loc[train.pdb_filename == "1_" + file_base, 'rmsd_relax'] = super_imposer.rms
                else:
                    file_array = file.split('_')
                    file_base = "_".join(file_array[pdbid_count:])
                    if file_base.endswith("mini.pdb"):
                        train.loc[
                            train.pdb_filename == "1_" + file_base, 'rmsd_relax'] = 0.00
                    else:
                        train.loc[train.pdb_filename == "1_" + file_base, 'rmsd_relax'] = 0.00

                procheck_dict = procheck(file, read_pdb_resolution()[file_base.split("_")[pdbid_count]])
                train.loc[train.pdb_filename == "1_" + file_base, 'disall'] = procheck_dict['disall']
                train.loc[train.pdb_filename == "1_" + file_base, 'bad_contacts'] = procheck_dict['bad_contacts']
                train.loc[train.pdb_filename == "1_" + file_base, 'bond_lenangle'] = procheck_dict['bond_lenangle']
                train.loc[train.pdb_filename == "1_" + file_base, 'g_factors'] = procheck_dict['g_factors']
                train.loc[train.pdb_filename == "1_" + file_base, 'bond_lengths_highlighted'] = procheck_dict[
                    'bond_lengths_highlighted']
                if 'bond_lengths_off' in procheck_dict:
                    train.loc[train.pdb_filename == "1_" + file_base, 'bond_lengths_off'] = procheck_dict[
                        'bond_lengths_off']
                else:
                    train.loc[train.pdb_filename == "1_" + file_base, 'bond_lengths_off'] = 0.00
                train.loc[train.pdb_filename == "1_" + file_base, 'bond_angles_highlighted'] = procheck_dict[
                    'bond_angles_highlighted']
                if 'bond_angles_off' in procheck_dict:
                    train.loc[train.pdb_filename == "1_" + file_base, 'bond_angles_off'] = procheck_dict[
                        'bond_angles_off']
                else:
                    train.loc[train.pdb_filename == "1_" + file_base, 'bond_angles_off'] = 0.00

    train.to_csv("test_out.csv", sep=';', encoding='utf-8')

    if train["pdb_filename"].str.endswith("mini.pdb").any():
        if train["pdb_filename"].str.startswith("minimized").any():
            new = train["pdb_filename"].str.split("_", n=13, expand=True)
            pdbid_start = 4
            amplitude_start = 10
            repeat_start = 3
            mode_start = 6
        else:
            new = train["pdb_filename"].str.split("_", n=10, expand=True)
            pdbid_start = 1
            amplitude_start = 7
            repeat_start = 0
            mode_start = 3
    else:
        new = train["pdb_filename"].str.split("_", n=9, expand=True)
        pdbid_start = 1
        amplitude_start = 7
        repeat_start = 0
        mode_start = 3
    train['pdbid'] = new[pdbid_start]
    # new data frame with split value columns
    train['amplitude'] = new[amplitude_start]
    train['amplitude'] = train['amplitude'].astype(float)

    train['repeat'] = new[repeat_start]
    train['repeat'] = train['repeat'].astype(int)

    train['mode'] = new[mode_start].str.slice(1)
    train['mode'] = train['mode'].astype(int)

    column_names_to_normalize = ['score_init']

    min_max_scaler = preprocessing.MinMaxScaler()
    x = train[column_names_to_normalize].values
    x_scaled = min_max_scaler.fit_transform(x)
    df_temp = pd.DataFrame(x_scaled, columns=column_names_to_normalize, index=train.index)
    train[column_names_to_normalize] = df_temp

    col_name = ['pdb_filename', 'rmsd_init', 'score_init', 'rmsd_relax', 'score_relax', 'pdbid', 'amplitude']

    train = train.loc[(train['pdbid'] != "5isv")]
    train = train.loc[(train['pdbid'] != "2cns")]

    # Remove the other repeats for now
    train = train.loc[(train['repeat'] == 1)]

    grouped = train.groupby(["pdbid"])
    train = grouped.apply(lambda x: x.sort_values(["amplitude"], ascending=True)).reset_index(drop=True)

    train['mean'] = grouped['score_relax'].transform('mean')
    train['std'] = grouped['score_relax'].transform('std')

    modes = [7, 8, 9, 10, 11, 12]
    for m in modes:
        train_mode = train.loc[(train['mode'] == m)]

        h = sns.FacetGrid(train_mode, col="pdbid", palette='seismic', sharey=False, sharex=True, col_wrap=4, height=2,
                          aspect=1)

        vmin = train_mode['rmsd_relax'].min()
        vmax = train_mode['rmsd_relax'].max()

        cmap = sns.light_palette("seagreen", as_cmap=True)

        h.map(facet_scatter, "amplitude", "score_relax", "rmsd_relax", s=100, vmin=vmin, vmax=vmax, cmap=cmap)

        h.map(facet_line, "score_relax")

        # Define a new Axes where the colorbar will go, left bottom, width, height
        cax = h.fig.add_axes([.40, .0339, .2, .023])

        # Get a mappable object with the same colormap as the data
        points = plt.scatter([], [], c=[], vmin=vmin, vmax=vmax, cmap=cmap)
        h.map(plt.plot, "amplitude", "score_relax")

        # Draw the colorbar
        cbar = h.fig.colorbar(points, cax=cax, orientation='horizontal')
        cbar.ax.set_title("RMSD relax" + " ($\AA$)", fontsize=10)

        plt.show()

        meltCov = pd.melt(train_mode, id_vars=['pdb_filename', 'rmsd_init', 'score_init', 'rmsd_relax', 'score_relax',
                                               'pdbid', 'amplitude', 'repeat', 'mode', 'mean', 'std'],
                          var_name='procheck')

        group_pdbid = meltCov.groupby(["pdbid"])
        for i in range(0, len(list(group_pdbid))):
            result = list(group_pdbid)[i][1]

            g = sns.FacetGrid(result, col='procheck', palette='seismic', sharey=False, sharex=True, col_wrap=4,
                              height=2)
            g.map(plt.scatter, 'amplitude', 'value')
            g.map(plt.plot, 'amplitude', 'value')

            g.fig.suptitle('{0} - mode{1}'.format(result['pdbid'].iloc[-1], m))

            g.fig.subplots_adjust(top=.85)
            axes = g.axes.flatten()
            axes[0].set_title("Disallowed")
            axes[1].set_title("Bad contacts")
            axes[2].set_title("Bond len/angle")
            axes[3].set_title("G-factors")
            axes[4].set_title("Bond len hlghtd")
            axes[5].set_title("Bond len off")
            axes[6].set_title("Bond angle hlghtd")
            axes[7].set_title("Bond angle off")

            plt.subplots_adjust(wspace=0.4)

            plt.savefig("Fig_mode{0}_procheck_{1}.png".format(m, result['pdbid'].iloc[-1]),
                        bbox_inches='tight', pad_inches=0.4, dpi=300)
            # plt.show()
