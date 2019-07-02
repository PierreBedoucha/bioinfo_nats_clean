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
from sklearn import preprocessing


def read_pfam_align():
    import os
    file_path = os.path.join("../data/input/etc", "pfam_env.txt")
    pdb_align_dict = {}
    with open(file_path) as f1:
        for line in f1:
            if not line.startswith("#") and not line.startswith("\n"):
                pdb_align_dict[line[0:4]] = (int(line[15:17]), int(line[21:24]))
    return pdb_align_dict


def compute_rmsd(pdb_path1, pdb_path2, start, end):
    sum_dist_sq = 0
    atom_cpt = 1
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


if __name__ == '__main__':
    train = pd.read_csv('Workbook1.csv')

    dict_ref_SC = {}
    dict_ref_relax = {}

    for file in os.listdir("."):
        if file.endswith(".pdb"):
            if 'minimized' not in file:
                start, end = read_pfam_align()[file[0:4]]
                file_ref_temp = file.replace(file.split("_")[-2], "0.00")
                file_ref = file_ref_temp.replace(file_ref_temp.split("_")[3], "a20.00")
                if file.split("_")[-2] != "0.00":
                    train.loc[train.pdb_filename == file, 'rmsd_init'] = compute_rmsd(file, file_ref, start, end)
                else:
                    train.loc[train.pdb_filename == file, 'rmsd_init'] = 0.00
            elif 'minimized' in file:
                start, end = read_pfam_align()[file.split("_")[1]]
                file_ref_temp = file.replace(file.split("_")[-2], "0.00")
                file_ref = file_ref_temp.replace(file_ref_temp.split("_")[4], "a20.00")
                if file.split("_")[-2] != "0.00":
                    # Align Relaxed structures and use rmsd align (CA only and pfam sequences)
                    select_CA_align(file, start, end)
                    res_to_be_aligned = ca_align_list
                    pdb_parser = Bio.PDB.PDBParser(QUIET=True)
                    # Get the structures
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

                    # train.loc[train.pdb_filename == file, 'rmsd_init'] = compute_rmsd(file, file_ref, start, end)
                    file = file.replace(file.split('_')[0] + '_', '')
                    train.loc[train.pdb_filename == file, 'rmsd_relax'] = super_imposer.rms
                else:
                    file = file.replace(file.split('_')[0] + '_', '')
                    train.loc[train.pdb_filename == file, 'rmsd_relax'] = 0.00

    # train.set_index('pdb_filename', inplace=True)
    train.to_csv("test_out.csv", sep=';', encoding='utf-8')

    # train['pdbid'] = train["pdb_filename"].str.split('_')[0]
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
    train = grouped.apply(lambda x: x.sort_values(["amplitude"], ascending = True)).reset_index(drop=True)

    # sns.set_style("whitegrid", {'axes.grid': False, 'axes.edgecolor': 'none'})

    # h = sns.FacetGrid(train, col="pdbid", hue='pdbid', col_wrap=7, sharey='row', sharex='col', margin_titles=True)
    # h = sns.FacetGrid(train, col="pdbid", palette = 'seismic', gridspec_kws={"hspace":0.4}, sharey=False, sharex=True)
    h = sns.FacetGrid(train, col="pdbid", palette='seismic', sharey=False, sharex=True, col_wrap=6, height=2, aspect=1)
    # h.map(f, "amplitude", "score_init", "rmsd_init", scale=.7, markers="")

    vmin = train['rmsd_relax'].min()
    vmax = train['rmsd_relax'].max()
    # vmin = train['rmsd_init'].min()
    # vmax = train['rmsd_init'].max()
    cmap = sns.diverging_palette(240, 10, l=65, center="dark", as_cmap=True)

    # h.map(plt.plot, "amplitude", "score_relax", marker="o")
    h.map(facet_scatter, "amplitude", "score_relax", "rmsd_relax", s=100, alpha=0.5, vmin=vmin, vmax=vmax, cmap=cmap)

    # Make space for the colorbar
    # h.fig.subplots_adjust(right=.92)
    # plt.tight_layout()

    # Define a new Axes where the colorbar will go
    # cax = h.fig.add_axes([.94, .25, .02, .6])
    cax = h.fig.add_axes([.40, .0339, .2, .023])

    # Get a mappable object with the same colormap as the data
    points = plt.scatter([], [], c=[], vmin=vmin, vmax=vmax, cmap=cmap)
    h.map(plt.plot, "amplitude", "score_relax")

    # Draw the colorbar
    cbar = h.fig.colorbar(points, cax=cax, orientation='horizontal')
    cbar.ax.set_title('RMSD relax', fontsize=10)

    plt.show()
