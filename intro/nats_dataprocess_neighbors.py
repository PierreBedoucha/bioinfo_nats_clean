import Bio.PDB
import Bio.AlignIO as al
import pandas as pd
import os

"""
This script summarizes the energy score of specific residues in the available structures (relaxed pdb file) 
in the current directory.
The structures scores are listed in the relaxed pdb files and the script will isolate only the the resiudes neighbouring
the catalytic sites. The resulting score is summed over this selection.
The resulting data is logged in a csv file 'pyrosetta_out.csv' for further analyses.
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
                    res_count = ls.index(i)
                    pdb_align_dict[rec.id[0:4]].append(res_count + read_pdb_starts()[rec.id[0:4]])
    return pdb_align_dict


def read_pdb_starts():
    """
    Reads at which index each pdb sequence is starting from the pdb_starts.txt file from ../data/input/etc
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


def read_pdb_catalytic():
    """
    Reads the list of catalytic residues listed in pdb_catalytic.txt file in ../data/input/etc
    :return: Dictionary. Keys: structure pdb id, Values: catalytic residue index
    :rtype: dict
    """
    file_path = os.path.join("../data/input/etc", "pdb_catalytic.txt")
    pdb_catalytic_dict = {}
    with open(file_path) as file_catalytic:
        for line_pdb in file_catalytic:
            if not line_pdb.startswith("#"):
                if not line_pdb.startswith('\n') or not line_pdb == '\n':
                    line_array_cat = line_pdb.split(':')
                    if len(line_array_cat) > 1 and (not line_array_cat[1].startswith('\n')
                                                    and not line_array_cat[1] == ''):
                        line_catalytic_array = line_array_cat[1].split(',')
                        pdb_catalytic_dict[line_pdb[0:4]] = [int(x) for x in line_catalytic_array]
    return pdb_catalytic_dict


def neighbor_res_select(pdbfilename, sel_chain, cutoff):
    """
    Selection of the catalytic residue neighbours from the structure pdb files. It uses NeighborSearch method from
    BioPDB. All the detected residues indices are returned in a list.
    :param pdbfilename: Pdb file name
    :type pdbfilename: str
    :param sel_chain: Identification letter of the structure chain
    :type sel_chain: str
    :param cutoff: Cutoff value in Angstroms for neighbouring search
    :type cutoff: float
    :return: List of detected residue indices
    :rtype: list
    """
    parser = Bio.PDB.PDBParser(QUIET=True)  # QUIET=True avoids comments on errors in the pdb.
    target_atoms = []
    res_list = []
    structures = parser.get_structure('1prot', pdbfilename)
    structure = structures[0]  # 'structures' may contain several proteins in this case only one.
    pdbid = pdbfilename.split('_')[-8]
    if pdbid in read_pdb_catalytic().keys():
        for res in read_pdb_catalytic()[pdbid]:
            try:
                target_atoms.append(structure[sel_chain.upper()][res]['CA'])
            except KeyError as err:
                print(err.args)
                target_atoms.append(structure['A'][res]['CA'])

    atoms = Bio.PDB.Selection.unfold_entities(structure, 'A')
    ns = Bio.PDB.NeighborSearch(atoms)

    for at in target_atoms:
        res_list.extend(ns.search(at.coord, cutoff))

    res_list = [x.parent.id[1] for x in res_list]
    res_list = list(set(res_list))
    return res_list


# Global variables
pdbfile_list = []
score_init_list = []
score_relax_list = []
score_relax_dict = {}
nb_of_repeats = 5

if __name__ == '__main__':
    import csv
    import numpy as np
    import re

    res_cpt_list = []
    res_scores_list = []
    res_score_dict = {}
    df = pd.DataFrame(columns=['res_cpt', 'res_scores'])

    for file in os.listdir("."):
        if 'minimized' in file:
            search = re.compile(r'[^A-Z]*_[0-9]+\s').search
            res_cpt = read_pdb_starts()[file.split('_')[-8]]
            with open(file) as pdb_f1:
                for pdb_line in pdb_f1:
                    if search(pdb_line):
                        pdb_line_array = pdb_line.split(' ')
                        res_cpt_list.append(res_cpt)
                        res_scores_list.append(float(pdb_line_array[-1].rstrip()))
                        res_score_dict[res_cpt] = float(pdb_line_array[-1].rstrip())

                        res_cpt += 1
            # Filter the residue at cutoff distance from active site
            # Read active site for atom selection
            chain = file.split('_')[-7]
            res_list_cutoff = neighbor_res_select(file, chain, 5)
            list_add = []
            for res_cut in res_list_cutoff:
                list_add.append(res_score_dict[res_cut])
            if list_add:
                score_final = np.sum(list_add) / len(list_add)
                score_relax_dict[file] = score_final

    wtr = csv.writer(open('pyrosetta_out.csv', 'w+'), delimiter=',', lineterminator='\n', quoting=csv.QUOTE_NONE,
                     escapechar='\\')  #
    wtr.writerow(['pdb_filename', 'score_relax'])

    score_init_list.append('0.0')

    padded_list = [0] * (nb_of_repeats * len(score_init_list))  # [0,0,0,0,0,0,0,0,0,0,0]
    padded_list[::nb_of_repeats] = score_init_list
    padded_list_rmsd = [0] * (nb_of_repeats * len(score_init_list))

    list_write = [(x for x in score_relax_dict.keys()), (x for x in score_relax_dict.values())]
    wtr.writerows([i for i in zip(*list_write)])
