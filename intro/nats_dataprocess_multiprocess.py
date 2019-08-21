from tempfile import mkstemp
from shutil import move
from os import fdopen, remove

from pyrosetta import *
from rosetta.core.scoring import *
from rosetta.protocols.relax import *
from rosetta.core.pose import *
from rosetta.protocols.constraint_movers import *
import multiprocessing as mp
from itertools import repeat
from functools import partial
import psutil


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

rosetta_options = ["-ignore_unrecognized_res false",
                   # "-ex1",
                   # "-ex2",
                   "-use_input_sc",
                   "-flip_HNQ",
                   "-no_optH false",
                   "-relax:constrain_relax_to_start_coords",
                   "-relax:coord_constrain_sidechains",
                   "-relax:ramp_constraints false",
                   "-constant_seed",
                   "-no_his_his_pairE",
                   "-linmem_ig 10",
                   "-nblist_autoupdate true",
                   "-relax:coord_cst_stdev 0.5"]

def score_proteins(matches, pdb_filename):
    # init(extra_options = "-constant_seed -ignore_unrecognized_res -ex2 -use_input_sc -no_his_his_pairE -no_optH false -flip_HNQ -relax:sc_cst_maxdist 3 -relax:coord_cst_stdev 0.1 -relax:coord_cst_width 1.0")
    scorefxn = pyrosetta.rosetta.core.scoring.ScoreFunctionFactory.create_score_function('ref2015_cst')
    pose = Pose()
    pose_from_file(pose, pdb_filename)

    relax = FastRelax(standard_repeats=5)
    relax.set_scorefxn(scorefxn)
    relax.apply(pose)
    pose.dump_pdb("minimized_fast_CST_" + pdb_filename)
    pose_score_2 = scorefxn(pose)

    matches.append(pose_score_2 / pyrosetta.rosetta.core.pose.Pose.total_residue(pose))


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


# Global variables
pdbfile_list = []
score_init_list = []
score_relax_list = []
score_relax_dict = {}
nb_of_repeats = 1

if __name__ == '__main__':
    pdb_occupancy()

    list_files = [x for x in os.listdir(".") if x.endswith(".pdb") and 'minimized' not in x]
    init(extra_options=" ".join(rosetta_options))
    scorefxn = pyrosetta.rosetta.core.scoring.ScoreFunctionFactory.create_score_function('ref2015_cst')
    for file in list_files:
        pose = Pose()
        pose_from_file(pose, file)
        pose_score = scorefxn(pose)
        score_init_list.append(pose_score / pyrosetta.rosetta.core.pose.Pose.total_residue(pose))

    # with mp.Pool(mp.cpu_count()) as pool:
    #     # results = pool.starmap(score_proteins, zip(list_file, repeat(scorefxn)))
    #     results = pool.starmap(score_proteins, list_file)
    # result_cpt = 1
    # for res, pdb in zip(results, list_file):
    #     score_relax_dict[pdb + "_" + str(result_cpt)] = res
    #     result_cpt += 1

    manager = mp.Manager()  # create SyncManager
    matches = manager.list()  # create a shared list here
    link_matches = partial(score_proteins, matches)  # create one arg callable to
    # pass to pool.map()
    # pool = mp.Pool(mp.cpu_count())
    pool = mp.Pool(psutil.cpu_count(logical = False))
    pool.map(link_matches, list_files)  # apply partial to files list
    pool.close()
    pool.join()
    print(matches)

    import csv
    wtr = csv.writer(open('pyrosetta_out.csv', 'w+'), delimiter=',', lineterminator='\n', quoting=csv.QUOTE_NONE,
                     escapechar='\\')  #
    wtr.writerow(['pdb_filename', 'rmsd_init', 'score_init', 'rmsd_relax', 'score_relax'])

    padded_list = [0] * (nb_of_repeats * len(score_init_list))  # [0,0,0,0,0,0,0,0,0,0,0]
    padded_list[::nb_of_repeats] = score_init_list
    padded_list_rmsd = [0] * (nb_of_repeats * len(score_init_list))

    l = [(x for x in list_files), (x for x in padded_list_rmsd), (x for x in padded_list),
         (x for x in padded_list_rmsd), (x for x in matches)]
    wtr.writerows([i for i in zip(*l)])