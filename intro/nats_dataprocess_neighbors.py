
from pyrosetta import *
from rosetta.core.scoring import *
from rosetta.protocols.relax import *
from rosetta.core.pose import *
from rosetta.protocols.constraint_movers import *
import Bio.PDB
import Bio.AlignIO as al

import pandas as pd

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

def read_pdb_starts():
    file_path = os.path.join("../data/input/etc", "pdb_starts.txt")
    pdb_starts_dict = {}
    with open(file_path) as f1:
        for line in f1:
            if not line.startswith("#") and not line.startswith("\n"):
                line_array = line.split(',')
                pdb_starts_dict[line[0:4]] = int(line_array[1])
    return pdb_starts_dict

def read_pdb_catalytic():
    file_path = os.path.join("../data/input/etc", "pdb_catalytic.txt")
    pdb_catalytic_dict = {}
    with open(file_path) as f1:
        for line in f1:
            if not line.startswith("#"):
                if not line.startswith('\n') or not line == '\n':
                    line_array = line.split(':')
                    if len(line_array) > 1 and (not line_array[1].startswith('\n') and not line_array[1] == ''):
                        line_catalytic_array = line_array[1].split(',')
                        pdb_catalytic_dict[line[0:4]] = [int(x) for x in line_catalytic_array]
    return pdb_catalytic_dict

def neighbor_res_select(pdbfilename, chain, cutoff):
    parser = Bio.PDB.PDBParser(QUIET=True) # QUIET=True avoids comments on errors in the pdb.
    target_atoms = []
    res_list = []
    structures = parser.get_structure('1prot', pdbfilename)
    structure = structures[0] # 'structures' may contain several proteins in this case only one.
    pdbid = pdbfilename.split('_')[-8]
    if pdbid in read_pdb_catalytic().keys():
        for res in read_pdb_catalytic()[pdbid]:
            try:
                target_atoms.append(structure[chain][res]['CA'])
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
            with open(file) as f1:
                for line in f1:
                    if search(line):
                        line_array = line.split(' ')
                        res_cpt_list.append(res_cpt)
                        res_scores_list.append(float(line_array[-1].rstrip()))
                        res_score_dict[res_cpt] = float(line_array[-1].rstrip())

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
    # wtr.writerow(['pdb_filename', 'score_init', 'score_relax'])
    wtr.writerow(['pdb_filename', 'score_relax'])

    score_init_list.append('0.0')

    padded_list = [0] * (nb_of_repeats * len(score_init_list))  # [0,0,0,0,0,0,0,0,0,0,0]
    padded_list[::nb_of_repeats] = score_init_list
    padded_list_rmsd = [0] * (nb_of_repeats * len(score_init_list))


    # l = [(x for x in score_relax_dict.keys()), (x for x in score_init_list), (x for x in score_relax_dict.values())]
    # l = [(x for x in score_relax_dict.keys()), (x for x in padded_list), (x for x in score_relax_dict.values())]
    l = [(x for x in score_relax_dict.keys()), (x for x in score_relax_dict.values())]
    # l = [(x for x in score_relax_dict.keys()), (x for x in padded_list_rmsd), (x for x in padded_list),
    #      (x for x in padded_list_rmsd), (x for x in score_relax_dict.values())]
    wtr.writerows([i for i in zip(*l)])