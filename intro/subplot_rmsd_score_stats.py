import pandas as pd
import seaborn as sns
from tempfile import mkstemp
from shutil import move
import os
import matplotlib.pyplot as plt
import math
from scipy.spatial import distance
import numpy as np
import Bio.PDB
import Bio.AlignIO as al
from sklearn import preprocessing
import glob
from pathlib import Path
import re
from os.path import realpath, basename

def read_pdb_resolution():
    file_path = os.path.join("../data/input/etc", "pdb_starts.txt")
    pdb_res_dict = {}
    with open(file_path) as f1:
        for line in f1:
            if not line.startswith("#") and not line.startswith("\n"):
                line_array = line.split(',')
                pdb_res_dict[line[0:4]] = int(line_array[1])
    return pdb_res_dict

def read_pfam_align():
    file_path = os.path.join("../data/input/etc", "pfam_env.txt")
    pdb_align_dict = {}
    with open(file_path) as f1:
        for line in f1:
            if not line.startswith("#") and not line.startswith("\n"):
                pdb_align_dict[line[0:4]] = (int(line[15:17]), int(line[21:24]))
    return pdb_align_dict

def read_pdb_starts():
    file_path = os.path.join("../data/input/etc", "pdb_starts.txt")
    pdb_starts_dict = {}
    with open(file_path) as f1:
        for line in f1:
            if not line.startswith("#") and not line.startswith("\n"):
                line_array = line.split(',')
                pdb_starts_dict[line[0:4]] = int(line_array[1])
    return pdb_starts_dict

def read_msa_fasta():
    pdb_align_dict = {'3tfy':[],'5isv':[],'4pv6':[],'2z0z':[],'1s7l':[],'2x7b':[],'3igr':[],'5k18':[],'2cns':[],
                      '5hh0':[],'5wjd':[],'5icv':[],'4kvm':[],'4u9v':[],}
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

def compute_rmsd_align(pdb_path1, pdb_path2):
    sum_dist_sq = 0
    atom_cpt = 1
    if pdb_path1.endswith("mini.pdb") and pdb_path2.endswith("mini.pdb"):
        snap_count = -4
        pdbid_count = -9
    else:
        snap_count = -3
        pdbid_count = -8
    file1_ref_array = pdb_path1.split('_')
    file2_ref_array = pdb_path2.split('_')
    if not os.path.exists(pdb_path1):
        file1_ref_array[snap_count] = "1"
        pdb_path1 = "_".join(file1_ref_array)
    if not os.path.exists(pdb_path2):
        file2_ref_array[snap_count] = "1"
        pdb_path2 = "_".join(file2_ref_array)
    # pdb1_res_list = read_msa_fasta()[file1_ref_array[-8]]
    # pdb2_res_list = read_msa_fasta()[file2_ref_array[-8]]
    pdb1_res_list = ca_align_dict[file1_ref_array[pdbid_count]]
    pdb2_res_list = ca_align_dict[file2_ref_array[pdbid_count]]
    with open(pdb_path1) as f1, open(pdb_path2) as f2:
        for line1, line2 in zip(f1, f2):
            # if line1[21:22] == pdb_current_chain and 'ATOM' in line1[0:6]:
            if 'ATOM' in line1[0:6] and ' CA ' in line1[12:16]:
                if 'ATOM' in line2[0:6] and ' CA ' in line2[12:16]:
                    if (int(line1[23:26].strip()) in pdb1_res_list) and \
                            (int(line2[23:26].strip()) in pdb2_res_list):
                        try:
                            dist = distance.cdist(
                                np.array([np.float64(val) for val in line1[31:54].split()]).reshape(1, -1),
                                np.array([np.float64(val) for val in line2[31:54].split()]).reshape(1, -1))
                        except ValueError as ex:
                            dist = distance.cdist(np.array([np.float64(line1[30:38]),
                                                            np.float64(line1[38:46]),
                                                            np.float64(line1[46:54])]).reshape(1, -1),
                                                  np.array([np.float64(line2[30:38]),
                                                            np.float64(line2[38:46]),
                                                            np.float64(line2[46:54])]).reshape(1, -1))
                        sum_dist_sq += math.pow(dist[0][0], 2)
                        atom_cpt += 1
    rmsd = math.sqrt(sum_dist_sq / atom_cpt)
    return rmsd

def compute_rmsd(pdb_path1, pdb_path2, start, end):
    sum_dist_sq = 0
    atom_cpt = 1
    if not os.path.exists(pdb_path1):
        file_ref_array = pdb_path1.split('_')
        file_ref_array[-3] = "1"
        pdb_path1 = "_".join(file_ref_array)
    if not os.path.exists(pdb_path2):
        file_ref_array = pdb_path2.split('_')
        file_ref_array[-3] = "1"
        pdb_path2 = "_".join(file_ref_array)
    with open(pdb_path1) as f1, open(pdb_path2) as f2:
        for line1, line2 in zip(f1, f2):
            # if line1[21:22] == pdb_current_chain and 'ATOM' in line1[0:6]:
            if 'ATOM' in line1[0:6] and ' CA ' in line1[12:16]:
                if 'ATOM' in line2[0:6] and ' CA ' in line2[12:16]:
                    if (start <= int(line1[23:26].strip()) <= end) and (start <= int(line2[23:26].strip()) <= end):
                        try:
                            dist = distance.cdist(
                                np.array([np.float64(val) for val in line1[31:54].split()]).reshape(1, -1),
                                np.array([np.float64(val) for val in line2[31:54].split()]).reshape(1, -1))
                        except ValueError as ex:
                            dist = distance.cdist(np.array([np.float64(line1[30:38]),
                                                            np.float64(line1[38:46]),
                                                            np.float64(line1[46:54])]).reshape(1, -1),
                                                  np.array([np.float64(line2[30:38]),
                                                            np.float64(line2[38:46]),
                                                            np.float64(line2[46:54])]).reshape(1, -1))
                        # print(dist[0][0])
                        sum_dist_sq += math.pow(dist[0][0], 2)
                        atom_cpt += 1
    rmsd = math.sqrt(sum_dist_sq / atom_cpt)
    return rmsd

def select_CA_align(pdb_path, start, end):
    with open(pdb_path) as f1:
        for line in f1:
            if 'ATOM' in line[0:6] and ' CA ' in line[12:16]:
                if start <= int(line[23:26].strip()) <= end:
                    # Append Atom id or Resid???
                    # ca_align_list.append(int(line[6:11].strip())) # Atom id
                    ca_align_list.append(int(line[23:26].strip()))  # Resid

# Global variables (Ugly)
ca_align_list = []
ca_align_dict = read_msa_fasta()


def f(x,y,z, **kwargs):
    ax = sns.pointplot(x,y,**kwargs)
    ax.axhline(5, alpha=0.5, color='grey')
    for i in range(len(x)):
        ax.annotate('{:6.2f}'.format(z.values[i]), xy=(i, z.values[i]),fontsize=8,
                    color=kwargs.get("color","k"),
                    bbox=dict(pad=.9,alpha=1, fc='w',color='none'),
                    va='center', ha='center',weight='bold')

def facet_scatter(x, y, c, **kwargs):
    """Draw scatterplot with point colors from a faceted DataFrame columns."""
    kwargs.pop("color")
    plt.scatter(x, y, c=c, **kwargs)

def facet_stairs(x, y, **kwargs):
    """Draw scatterplot with point colors from a faceted DataFrame columns."""
    plt.scatter(x, y, **kwargs)

def facet_line(y, **kwargs):
    # plt.axhline(y.mean(), linestyle="--",
    #             color='gray')
    # plt.axhspan(y.mean() + y.std(), y.min(), facecolor='gray', alpha=0.1)

    # p_25, p_75 = np.percentile(y, [25, 75])
    # iqr = p_75 - p_25
    # upper_bound = p_75 + 1.5 * iqr
    # lower_bound = p_25 - 1.5 * iqr
    # plt.axhline(y.median(), linestyle="--",
    #             color="gray")
    # plt.axhspan(p_75, y.min(), facecolor='gray', alpha=0.1)
    #
    # t = plt.text(3, p_75, round(p_75, 2), horizontalalignment='right',
    #              verticalalignment='center', color='gray')

    plt.axhline(-2.0, linestyle="--",
                color='gray')
    t = plt.text(3, -2.0, -2.0, horizontalalignment='right',
                 verticalalignment='center', color='gray')

def procheck(sc_pdbpath, pdb_resolution):
    import shutil
    import subprocess
    os.environ['prodir'] = '~/Software/procheck/procheck'
    if sc_pdbpath[:-4] + ".sum" not in os.listdir("."):
        p = subprocess.run("{0} {1} {2}".format('~/Software/procheck/procheck/procheck.scr', sc_pdbpath, pdb_resolution),
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
                    output_dict['disall'] =  float(''.join(chars[:-1]))
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

if __name__ == '__main__':
    train = pd.read_csv('Workbook31.csv').dropna()

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
                    # train.loc[train.pdb_filename == file, 'rmsd_init'] = compute_rmsd(file, file_ref, start, end)
                    if file.endswith("mini.pdb"):
                        train.loc[train.pdb_filename == "1_" + file, 'rmsd_init'] = compute_rmsd_align(file, file_ref)
                    else:
                        train.loc[train.pdb_filename == "1_" + file, 'rmsd_init'] = compute_rmsd_align(file, file_ref)
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
                    # res_to_be_aligned = ca_align_list
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
                    # train.loc[train.pdb_filename == file_base, 'rmsd_relax'] = compute_rmsd(file_base, file_ref, start, end)
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

    # train.set_index('pdb_filename', inplace=True)
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
    train = grouped.apply(lambda x: x.sort_values(["amplitude"], ascending = True)).reset_index(drop=True)

    # sns.set_style("whitegrid", {'axes.grid': False, 'axes.edgecolor': 'none'})

    train['mean'] = grouped['score_relax'].transform('mean')
    train['std'] = grouped['score_relax'].transform('std')

    # modes = [7]
    modes = [7, 8, 9, 10, 11, 12]
    # modes = [10, 11, 12]
    for m in modes:
        train_mode = train.loc[(train['mode'] == m)]

        # h = sns.FacetGrid(train, col="pdbid", hue='pdbid', col_wrap=7, sharey='row', sharex='col', margin_titles=True)
        # h = sns.FacetGrid(train, col="pdbid", palette = 'seismic', gridspec_kws={"hspace":0.4}, sharey=False, sharex=True)
        h = sns.FacetGrid(train_mode, col="pdbid", palette='seismic', sharey=False, sharex=True, col_wrap=4, height=2, aspect=1)
        # h.map(f, "amplitude", "score_init", "rmsd_init", scale=.7, markers="")

        vmin = train_mode['rmsd_relax'].min()
        vmax = train_mode['rmsd_relax'].max()
        # vmin = train['rmsd_init'].min()
        # vmax = train['rmsd_init'].max()

        # cmap = sns.diverging_palette(150, 275, s=80, l=55, center="light", as_cmap=True)
        # cmap = sns.light_palette((44,162,95), input="husl", as_cmap=True)
        cmap = sns.light_palette("seagreen",  as_cmap=True)

        # h.map(plt.plot, "amplitude", "score_relax", marker="o")
        # h.map(facet_scatter, "amplitude", "score_relax", "rmsd_relax", s=100, alpha=0.5, vmin=vmin, vmax=vmax, cmap=cmap)
        h.map(facet_scatter, "amplitude", "score_relax", "rmsd_relax", s=100, vmin=vmin, vmax=vmax, cmap=cmap)

        h.map(facet_line, "score_relax")

        # Make space for the colorbar
        # h.fig.subplots_adjust(right=.92)
        # plt.tight_layout()

        # Define a new Axes where the colorbar will go, left bottom, width, height
        # cax = h.fig.add_axes([.94, .25, .02, .6])
        cax = h.fig.add_axes([.40, .0339, .2, .023])

        # Get a mappable object with the same colormap as the data
        points = plt.scatter([], [], c=[], vmin=vmin, vmax=vmax, cmap=cmap)
        h.map(plt.plot, "amplitude", "score_relax")

        # Draw the colorbar
        cbar = h.fig.colorbar(points, cax=cax, orientation='horizontal')
        cbar.ax.set_title("RMSD relax" + " ($\AA$)", fontsize=10)

        plt.show()
        # plt.savefig("Fig_mode{0}_rmsd.png".format(m), bbox_inches='tight', pad_inches=0.4, dpi=150)

        # i = sns.FacetGrid(train_mode, col="pdbid", palette='seismic', sharey=False, sharex=False, col_wrap=4, height=2,
        #                   aspect=1)
        # i.map(facet_stairs, "rmsd_relax", "disall")
        # i.map(facet_stairs, "rmsd_relax", "bad_contacts")
        # i.map(facet_stairs, "rmsd_relax", "bond_lenangle")
        # i.map(facet_stairs, "rmsd_relax", "g_factors")
        # i.map(facet_stairs, "rmsd_relax", "bond_lengths_highlighted")
        # i.map(facet_stairs, "rmsd_relax", "bond_lengths_off")
        # i.map(facet_stairs, "rmsd_relax", "bond_angles_highlighted")
        # i.map(facet_stairs, "rmsd_relax", "bond_angles_off")
        # plt.show()

        meltCov = pd.melt(train_mode, id_vars=['pdb_filename', 'rmsd_init', 'score_init', 'rmsd_relax', 'score_relax'
                                               ,'pdbid', 'amplitude', 'repeat', 'mode', 'mean', 'std'], var_name='procheck')
        # g = sns.FacetGrid(meltCov, col='pdbid', hue='procheck')
        # g = sns.FacetGrid(meltCov, col='pdbid', hue='procheck', sharey=False, sharex=True)

        group_pdbid = meltCov.groupby(["pdbid"])
        for i in range(0, len(list(group_pdbid))):
            result = list(group_pdbid)[i][1]

            g = sns.FacetGrid(result, col='procheck', palette='seismic', sharey=False, sharex=True, col_wrap=4, height=2)
            g.map(plt.scatter, 'amplitude', 'value')
            g.map(plt.plot, 'amplitude', 'value')
            # g.fig.tight_layout()
            g.fig.suptitle('{0} - mode{1}'.format(result['pdbid'].iloc[-1], m))
            # plt.title('{0} - mode{1}'.format(result['pdbid'].iloc[-1], m), loc='left')
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
            # plt.subplots_adjust(hspace=0.4, wspace=0.4)
            plt.subplots_adjust(wspace=0.4)
            # g.map(sns.lineplot, 'amplitude', 'value', markers=True)
            # g.set_xticklabels(rotation=45)
            # g.add_legend()
            # plt.legend(loc='lower left')
            # g.fig.get_axes()[0].legend(loc='lower left')
            # plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
            #           ncol=4, mode="expand", borderaxespad=0.)
            # plt.legend(bbox_to_anchor=(0., 1.02, 0.5, .102), loc='lower left',
            #            ncol=2, mode="expand", borderaxespad=0.)
            plt.savefig("Fig_mode{0}_procheck_{1}.png".format(m, result['pdbid'].iloc[-1]),
                        bbox_inches='tight', pad_inches=0.4, dpi=300)
            # plt.show()

