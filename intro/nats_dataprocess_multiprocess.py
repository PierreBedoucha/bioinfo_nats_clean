"""Energy scoring of the NATSs protein models. MULTIPROCESS

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
import multiprocessing as mp
from functools import partial
import psutil


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


def score_proteins(matches, pdb_filename):
    """Base structure scoring function.

    Here only the relaxation with FastRelax algorithm is performed and the relaxed pose is scored and passed through
    the mp channel to handle replicates.

    :param matches: shared list (channel)
    :param pdb_filename: pdb file name
    """
    scorefxn = pyrosetta.rosetta.core.scoring.ScoreFunctionFactory.create_score_function('ref2015_cst')
    pose = Pose()
    pose_from_file(pose, '_'.join(pdb_filename.split("_")[1:]))

    relax = FastRelax(standard_repeats=5)
    relax.set_scorefxn(scorefxn)
    relax.apply(pose)
    pose.dump_pdb("minimized_fast_CST_" + pdb_filename)
    pose_score_2 = scorefxn(pose)

    matches.append((pdb_filename, pose_score_2 / pyrosetta.rosetta.core.pose.Pose.total_residue(pose)))


def pdb_occupancy():
    """Cleans the pdb files in the current directory by quickly replacing with its fixed version.

    Each cleaned pdb filename is appended to a list to later log the data.
    """
    import os
    # Read pdbid // chain mapping
    pdb_chains_dict = read_pdb_chains()
    chains = None
    for pdb_file in os.listdir("."):
        if pdb_file.endswith(".pdb") and 'minimized' not in pdb_file:
            fh, abs_path = mkstemp()
            if pdb_file[0:4] in pdb_chains_dict.keys():
                chains = pdb_chains_dict[pdb_file[0:4]]
            with fdopen(fh, 'w') as new_file:
                with open(pdb_file) as f1:
                    for line in f1:
                        if line[0:6] == 'ATOM  ':
                            temp = line[50:60].replace(line[54:60], '{:6.2f}'.format(1.0))
                            line = line.replace(line[50:60], temp)
                            if not line.endswith('\n'):
                                line = line + '\n'
                            if line[17:20] == 'HSD' or line[17:20] == 'HSE':
                                line = line.replace(line[17:20], 'HIS')
                            if line[17:20] == 'SER' and line[12:16] == ' O  ':
                                line = line.replace(' O  ', ' OXT')
                            if line[17:20] == 'CYX' or line[17:20] == ' CYM':
                                line = line.replace(line[17:20], 'CYS ')
                            if chains is None:
                                new_file.write("%s" % line)
                            elif line[21:22] in chains and len(chains) == 1:
                                # line = line.replace(line[21:22], "{0}".format(chains[0]))
                                new_file.write("%s" % line)
                            elif len(chains) == 1:
                                line = line.replace(line[21:30], chains[0] + line[22:30])
                                new_file.write("%s" % line)
            # Remove original file
            remove(pdb_file)
            # Move new file
            move(abs_path, pdb_file)

            pdbfile_list.append(pdb_file)


# Global variables
rosetta_options = ["-ignore_unrecognized_res false",
                   "-ex1",
                   "-ex2",
                   "-ex1aro",
                   "-ex2aro",
                   "-use_input_sc",
                   "-flip_HNQ",
                   "-no_optH false",
                   "-relax:constrain_relax_to_start_coords",
                   # "-relax:coord_constrain_sidechains",
                   "-relax:ramp_constraints false",
                   "-constant_seed",
                   "-no_his_his_pairE",
                   # "-linmem_ig 10",
                   "-nblist_autoupdate true",
                   "-relax:coord_cst_stdev 0.5"]
pdbfile_list = []
score_init_list = []
score_relax_list = []
score_relax_dict = {}
nb_of_repeats = 1

if __name__ == '__main__':
    pdb_occupancy()

    list_files_relaxed = [x for x in os.listdir(".") if x.endswith(".pdb") and 'minimized' in x]
    list_files = [x for x in os.listdir(".") if x.endswith(".pdb") and 'minimized' not in x
                  and ("minimized_fast_CST_{}_".format(str(nb_of_repeats)) + x) not in list_files_relaxed]
    if list_files:
        init(extra_options=" ".join(rosetta_options))
        scorefxn = pyrosetta.rosetta.core.scoring.ScoreFunctionFactory.create_score_function('ref2015_cst')
        for file in list_files:
            pose = Pose()
            pose_from_file(pose, file)
            pose_score = scorefxn(pose)
            score_init_list.append(pose_score / pyrosetta.rosetta.core.pose.Pose.total_residue(pose))

        manager = mp.Manager()  # create SyncManager
        shared_matches = manager.list()  # create a shared list here
        link_matches = partial(score_proteins, shared_matches)  # create one arg callable to

        # Repeating the relax/scoring
        list_files_rep = []
        for file in list_files:
            for i in list(range(1, nb_of_repeats + 1, 1)):
                list_files_rep.append(str(i) + "_" + file)

        pool = mp.Pool(psutil.cpu_count(logical=False))
        pool.map(link_matches, list_files_rep)  # apply partial to files list
        pool.close()
        pool.join()

        import csv

        wtr = csv.writer(open('pyrosetta_out.csv', 'w+'), delimiter=',', lineterminator='\n', quoting=csv.QUOTE_NONE,
                         escapechar='\\')  #
        wtr.writerow(
            ['pdb_filename', 'rmsd_init', 'score_init', 'rmsd_relax', 'disall', 'bad_contacts', 'bond_lenangle',
             'g_factors', 'bond_lengths_highlighted', 'bond_lengths_off', 'bond_angles_highlighted',
             'bond_angles_off',
             'score_relax'])

        padded_list = [0] * (nb_of_repeats * len(score_init_list))  # [0,0,0,0,0,0,0,0,0,0,0]
        padded_list[::nb_of_repeats] = score_init_list
        padded_list_rmsd = [0] * (nb_of_repeats * len(score_init_list))

        list_write = [(x[0] for x in shared_matches), (x for x in padded_list_rmsd), (x for x in padded_list),
                      (x for x in padded_list_rmsd), (x for x in padded_list_rmsd),
                      (x for x in padded_list_rmsd), (x for x in padded_list_rmsd),
                      (x for x in padded_list_rmsd), (x for x in padded_list_rmsd),
                      (x for x in padded_list_rmsd), (x for x in padded_list_rmsd),
                      (x for x in padded_list_rmsd), (x[1] for x in shared_matches)]
        wtr.writerows([i for i in zip(*list_write)])
