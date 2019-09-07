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
    file1_ref_array = pdb_path1.split('_')
    file2_ref_array = pdb_path2.split('_')
    if not os.path.exists(pdb_path1):
        file1_ref_array[-3] = "1"
        pdb_path1 = "_".join(file1_ref_array)
    if not os.path.exists(pdb_path2):
        file2_ref_array[-3] = "1"
        pdb_path2 = "_".join(file2_ref_array)
    # pdb1_res_list = read_msa_fasta()[file1_ref_array[-8]]
    # pdb2_res_list = read_msa_fasta()[file2_ref_array[-8]]
    pdb1_res_list = ca_align_dict[file1_ref_array[-8]]
    pdb2_res_list = ca_align_dict[file2_ref_array[-8]]
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


# GLOBAL VARIABLES
amplitude_max = 30

if __name__ == '__main__':
    train = pd.read_csv('Workbook20.csv')

    dict_ref_SC = {}
    dict_ref_relax = {}

    for file in os.listdir("."):
        if file.endswith(".pdb"):
            if 'minimized' not in file:
                # start, end = read_pfam_align()[file[0:4]]
                file_ref_temp = file.replace(file.split("_")[-2], "0.00")
                file_ref = file_ref_temp.replace(file_ref_temp.split("_")[3], "a{0}.00".format(str(amplitude_max)))
                if file.split("_")[-2] != "0.00":
                    # train.loc[train.pdb_filename == file, 'rmsd_init'] = compute_rmsd(file, file_ref, start, end)
                    train.loc[train.pdb_filename == "1_" + file, 'rmsd_init'] = compute_rmsd_align(file, file_ref)
                else:
                    train.loc[train.pdb_filename == "1_" + file, 'rmsd_init'] = 0.00
            elif 'minimized' in file:
                # start, end = read_pfam_align()[file.split("_")[-8]]
                file_ref_temp = file.replace(file.split("_")[-2], "0.00")
                file_ref = file_ref_temp.replace(file_ref_temp.split("_")[-5], "a{0}.00".format(str(amplitude_max)))
                if file.split("_")[-2] != "0.00":
                    # Align Relaxed structures and use rmsd align (CA only and pfam sequences)
                    # select_CA_align(file, start, end)
                    # res_to_be_aligned = ca_align_list
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
                    # train.loc[train.pdb_filename == file, 'rmsd_relax'] = compute_rmsd(file, file_ref, start, end)
                    train.loc[train.pdb_filename == "1_" + file, 'rmsd_relax'] = super_imposer.rms
                else:
                    file_array = file.split('_')
                    file = "_".join(file_array[-8:])
                    train.loc[train.pdb_filename == "1_" + file, 'rmsd_relax'] = 0.00

    # train.set_index('pdb_filename', inplace=True)
    train.to_csv("test_out.csv", sep=';', encoding='utf-8')

    new = train["pdb_filename"].str.split("_", n=9, expand=True)
    train['pdbid'] = new[1]
    # new data frame with split value columns
    train['amplitude'] = new[7]
    train['amplitude'] = train['amplitude'].astype(float)

    train['repeat'] = new[0]
    train['repeat'] = train['repeat'].astype(int)

    train['mode'] = new[3].str.slice(1)
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

    modes=[9,10,12]
    for m in modes:
        train_mode = train.loc[(train['mode'] == m)]

        # h = sns.FacetGrid(train, col="pdbid", hue='pdbid', col_wrap=7, sharey='row', sharex='col', margin_titles=True)
        # h = sns.FacetGrid(train, col="pdbid", palette = 'seismic', gridspec_kws={"hspace":0.4}, sharey=False, sharex=True)
        h = sns.FacetGrid(train_mode, col="pdbid", palette='seismic', sharey=False, sharex=True, col_wrap=6, height=2, aspect=1)
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
        cbar.ax.set_title('RMSD relax', fontsize=10)

        plt.show()
        # plt.savefig("Fig_mode{0}_rmsd.png".format(m), bbox_inches='tight', pad_inches=0.4, dpi=150)
