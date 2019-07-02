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

    # for file in os.listdir("."):
    #     if file.endswith(".pdb") and 'minimized' not in file:
    #         # X = train.loc[(train['pdb_filename'].split('_')[0] == file.split('_')[0])].copy()
    #
    #         X = train[train['pdb_filename'].str.contains(file.split('_')[0])]
    #         X.sort_values(by=['rmsd_init'])
    #
    #         # g = sns.FacetGrid(X, row='rmsd_init', col='score_init')
    #         # g.map(plt.scatter, "x", "y")
    #         ac = sns.lineplot(x="rmsd_init", y="score_init", markers="o", data=X)
    #         ac.set(xlabel='RMSD (A)', ylabel='Score')
    #         plt.title('Initial Rosetta score')
    #         plt.show()
    #         print("fin test")

    # file_list = [f for f in os.listdir(".") if f.endswith(".pdb") and 'minimized' not in f]

    # melted = train.melt(id_vars=['pdb_filename', 'rmsd_init'], value_vars=['score_init'])
    # g = sns.FacetGrid(melted, col='pdb_filename', hue='pdb_filename', row='variable', sharey='row', margin_titles=True)
    # g.map(plt.plot, 'rmsd_init', 'value')
    # plt.show()

    # train['pdbid'] = train["pdb_filename"].str.split('_')[0]
    train['pdbid'] = train.pdb_filename.str[:4]
    # new data frame with split value columns
    new = train["pdb_filename"].str.split("_", n=7, expand=True)
    train['amplitude'] = new[6]
    train = train.groupby(['pdbid']).sor

    g = sns.FacetGrid(train, col="pdbid", row='variable', height=1.5, hue='pdbid', sharey='row', margin_titles=True)
    g = g.map(plt.plot, "rmsd_init", "score_init", marker=".")
    plt.show()

    # trainB = train.groupby(train['pdbid'], as_index=False, sort=False)
    # groups = dict(list(trainB))
    print(train)