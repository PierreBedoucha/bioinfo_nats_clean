"""Energy scoring of the NATSs protein models

    This script computes the energy score of available structures (Pose object) in the current directory using pyROSETTA API

    The scoring function parameters are detailed in the score_proteins function. The structures are scored before and
    after a relaxation step done with the FastRelax algorithm.

    The resulting data is logged in a csv file 'pyrosetta_out.csv' for further analyses.
"""

from tempfile import mkstemp
from shutil import move
from os import fdopen, remove
from pyrosetta import *
from rosetta.core.scoring import *
from rosetta.protocols.relax import *
from rosetta.core.pose import *
from rosetta.protocols.constraint_movers import *


def read_pdb_chains():
    """Read the selected chains for the protein dataset.

    The data is parsed from pdb_chains.txt file in ../data/input/etc.

    :return: Dictionary. Keys: structure pdb id, Values: selected chain letter
    :rtype: dict
    """
    file_path = os.path.join("../data/input/etc", "pdb_chains.txt")
    pdb_chains_dict = {}
    with open(file_path) as f1:
        for line in f1:
            if not line.startswith("#") and not line.startswith("\n"):
                line_array = line.split(':')
                if len(line_array[1].split(',')) > 1:
                    chains_array = [x.rstrip().upper() for x in line_array[1].split(',')]
                    pdb_chains_dict[line[0:4].lower()] = chains_array
                else:
                    pdb_chains_dict[line[0:4].lower()] = [line_array[1].rstrip().upper()]
    return pdb_chains_dict


def score_proteins(pdb_filename):
    """Main structure scoring function.

    Describes the scoring parameters and set the values of two global lists containing
    the initial score before relaxation on one hand, and the final score after relaxation on the other.

    :param pdb_filename: pdb file name
    """
    init(extra_options="-constant_seed")
    # scorefxn = ScoreFunction()
    # scorefxn.set_weight(fa_atr, 0.800)
    # scorefxn.set_weight(fa_rep, 0.440)  # full-atom repulsive score
    # scorefxn.set_weight(fa_sol, 0.750)  # full-atom solvation score
    # scorefxn.set_weight(fa_intra_rep, 0.004)  # f.a. intraresidue rep. score
    # scorefxn.set_weight(fa_elec, 0.700)  # full-atom electronic score
    # scorefxn.set_weight(pro_close, 1.000)  # proline closure
    # scorefxn.set_weight(hbond_sr_bb, 1.170)  # short-range hbonding
    # scorefxn.set_weight(hbond_lr_bb, 1.170)  # long-range hbonding
    # scorefxn.set_weight(hbond_bb_sc, 1.170)  # backbone-sidechain hbonding
    # scorefxn.set_weight(hbond_sc, 1.100)  # sidechain-sidechain hbonding
    # scorefxn.set_weight(dslf_fa13, 1.000)  # disulfide full-atom score
    # scorefxn.set_weight(rama, 0.200)  # ramachandran score
    # scorefxn.set_weight(omega, 0.500)  # omega torsion score
    # scorefxn.set_weight(fa_dun, 0.560)  # fullatom Dunbrack rotamer score
    # scorefxn.set_weight(p_aa_pp, 0.320)
    # scorefxn.set_weight(ref, 1.000)  # reference identity score
    scorefxn = pyrosetta.rosetta.core.scoring.ScoreFunctionFactory.create_score_function('ref2015_cst')

    # Scorefunction constraints setup -  Already done with talaris2014_cst
    # score_manager = pyrosetta.rosetta.core.scoring.ScoreTypeManager()
    # constraint = score_manager.score_type_from_name('atom_pair_constraint')
    # scorefxn.set_weight(constraint, 5)

    # create a pose from the desired PDB file
    # -create an empty Pose object
    pose = Pose()
    # -load the data from pdb_file into the pose
    pose_from_file(pose, pdb_filename)
    # default to the median residue number
    pose_score = scorefxn(pose)
    score_init_list.append(pose_score / pyrosetta.rosetta.core.pose.Pose.total_residue(pose))

    # === Bolean for Relax here ===
    # relax = ClassicRelax()
    # relax.set_scorefxn(scorefxn)
    # relax.apply(pose)
    # =============================

    # Pose constraints setup
    # constraints = pyrosetta.rosetta.protocols.constraint_movers.ConstraintSetMover()
    # constraints.constraint_file('constraints.cst')
    # constraints.add_constraints(True)
    # constraints.apply(pose)

    relax = FastRelax(standard_repeats=10)
    relax.set_scorefxn(scorefxn)

    # ------
    # relax.repeats(10)
    # relax.nstruct(10)

    relax.constrain_relax_to_start_coords(True)
    relax.ramp_down_constraints(False)

    relax.apply(pose)

    pose.dump_pdb("minimized_fast_cst_" + pdb_filename)
    pose_score_2 = scorefxn(pose)

    score_relax_list.append(pose_score_2 / pyrosetta.rosetta.core.pose.Pose.total_residue(pose))


def replace(file_path, pattern, subst):
    """Small helper function to replace a str pattern by another one in a given file.

    :param file_path: File path of the file to consider

    :param pattern: str pattern to replace

    :param subst: str pattern for replacing
    """
    # Create temp file
    fh, abs_path = mkstemp()
    with fdopen(fh, 'w') as new_file:
        with open(file_path) as old_file:
            for line in old_file:
                new_file.write(line.replace(pattern, subst))
    # Remove original file
    remove(file_path)
    # Move new file
    move(abs_path, file_path)


def pdb_occupancy():
    """Cleans the pdb files in the current directory by quickly replacing its fixed version and launches the scoring

    Each cleaned pdb filename is appended to a list to later log the data.
    """
    import os
    # Read pdbid // chain mapping
    pdb_chains_dict = read_pdb_chains()
    chains = None
    for file in os.listdir("."):
        if file.endswith(".pdb") and 'minimized' not in file:
            fh, abs_path = mkstemp()
            if file[0:4] in pdb_chains_dict.keys():
                chains = pdb_chains_dict[file[0:4]]
            with fdopen(fh, 'w') as new_file:
                with open(file) as f1:
                    for line in f1:
                        if line[0:6] == 'ATOM  ':
                            temp = line[50:60].replace(line[54:60], '{:6.2f}'.format(1.0))
                            line = line.replace(line[50:60], temp)
                            # if line[12:16] == ' O  ' or line[12:16] == 'OCT1':
                            #     line = line.replace(' O  ', ' OT1')
                            # elif line[12:16] == ' OXT' or line[12:16] == 'OCT2':
                            #     old = line[12:16]
                            #     line = line.replace(old, ' OT2')
                            # if line[17:20] == 'ILE' and line[12:16] == ' CD1':
                            #     line = line.replace('CD1', 'CD ')
                            if line[17:20] == 'HSD' or line[17:20] == 'HSE':
                                line = line.replace(line[17:20], 'HIS')
                            if line[17:20] == 'SER' and line[12:16] == ' O  ':
                                line = line.replace(' O  ', ' OXT')
                            if line[17:20] == 'CYX' or line[17:20] == ' CYM':
                                line = line.replace(line[17:20], 'CYS ')
                            if chains is None:
                                new_file.write("%s" % line)
                            elif line[21:22] in chains:
                                new_file.write("%s" % line)
            # Remove original file
            remove(file)
            # Move new file
            move(abs_path, file)

            pdbfile_list.append(file)
            score_proteins(file)


# Global variables
pdbfile_list = []
score_init_list = []
score_relax_list = []

if __name__ == '__main__':
    pdb_occupancy()
    # score_proteins()
    import csv

    wtr = csv.writer(open('pyrosetta_out.csv', 'w+'), delimiter=',', lineterminator='\n', quoting=csv.QUOTE_NONE,
                     escapechar='\\')  #
    wtr.writerow(['pdb_filename', 'score_init', 'score_relax'])

    list_write = [(x for x in pdbfile_list), (x for x in score_init_list), (x for x in score_relax_list)]
    wtr.writerows([i for i in zip(*list_write)])
