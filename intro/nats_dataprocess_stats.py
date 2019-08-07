from tempfile import mkstemp
from shutil import move
from os import fdopen, remove

from pyrosetta import *
from rosetta.core.scoring import *
from rosetta.protocols.relax import *
from rosetta.core.pose import *
from rosetta.protocols.constraint_movers import *


def read_pdb_chains():
    file_path = os.path.join("../data/input/etc", "pdb_chains.txt")
    pdb_chains_dict = {}
    with open(file_path) as f1:
        for line in f1:
            if not line.startswith("#") and not line.startswith("\n"):
                chains_array = []
                line_array = line.split(':')
                if len(line_array[1].split(',')) > 1:
                    chains_array = [x.rstrip().upper() for x in line_array[1].split(',')]
                    pdb_chains_dict[line[0:4].lower()] = chains_array
                else:
                    pdb_chains_dict[line[0:4].lower()] = [line_array[1].rstrip().upper()]
    return pdb_chains_dict


def score_proteins(pdb_filename):
    init(extra_options = "-constant_seed")

    scorefxn = pyrosetta.rosetta.core.scoring.ScoreFunctionFactory.create_score_function('ref2015_cst')


    # create a pose from the desired PDB file
    # -create an empty Pose object
    pose = Pose()
    # -load the data from pdb_file into the pose
    pose_from_file(pose, pdb_filename)
    pose_score = scorefxn(pose)
    score_init_list.append(pose_score / pyrosetta.rosetta.core.pose.Pose.total_residue(pose))

    relax = FastRelax(standard_repeats=5)
    relax.set_scorefxn(scorefxn)

    relax.constrain_relax_to_start_coords(True)
    relax.ramp_down_constraints(False)

    for i in list(range(1, 501, 1)):
        relax.apply(pose)
        pose.dump_pdb("minimized_fast_cst_" + str(i) + "_" + pdb_filename)
        pose_score_2 = scorefxn(pose)

        # score_relax_list.append(pose_score_2 / pyrosetta.rosetta.core.pose.Pose.total_residue(pose))
        score_relax_dict[pdb_filename + "_" + str(i)] = pose_score_2 / pyrosetta.rosetta.core.pose.Pose.total_residue(pose)


def pdb_occupancy():
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
            break


# Global variables
pdbfile_list = []
score_init_list = []
score_relax_list = []
score_relax_dict = {}

if __name__ == '__main__':
    pdb_occupancy()
    import csv
    import numpy as np
    wtr = csv.writer(open('pyrosetta_out.csv', 'w+'), delimiter=',', lineterminator='\n', quoting=csv.QUOTE_NONE,
                     escapechar='\\')  #
    # wtr.writerow(['pdb_filename', 'score_init', 'score_relax'])
    wtr.writerow(['pdb_filename', 'score_relax'])

    padded_list = [0] * (10 * len(score_init_list))  # [0,0,0,0,0,0,0,0,0,0,0]
    padded_list[::10] = score_init_list

    # l = [(x for x in score_relax_dict.keys()), (x for x in score_init_list), (x for x in score_relax_dict.values())]
    # l = [(x for x in score_relax_dict.keys()), (x for x in padded_list), (x for x in score_relax_dict.values())]
    l = [(x for x in score_relax_dict.keys()), (x for x in score_relax_dict.values())]
    wtr.writerows([i for i in zip(*l)])