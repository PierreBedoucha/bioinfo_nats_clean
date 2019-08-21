from tempfile import mkstemp
from shutil import move
from os import fdopen, remove

from pyrosetta import *
from rosetta.core.scoring import *
from rosetta.protocols.relax import *
from rosetta.core.pose import *
from rosetta.protocols.constraint_movers import *
from pyrosetta.mpi import mpi_init
from PyRosetta_TACC_MPI import *
from mpi4py import MPI


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
    # init(extra_options = "-constant_seed -ignore_unrecognized_res -ex2 -use_input_sc -no_his_his_pairE -no_optH false -flip_HNQ -relax:sc_cst_maxdist 3 -relax:coord_cst_stdev 0.5")
    #
    # scorefxn = pyrosetta.rosetta.core.scoring.ScoreFunctionFactory.create_score_function('ref2015_cst')
    #
    # # create a pose from the desired PDB file
    # # -create an empty Pose object
    # pose = Pose()
    # # -load the data from pdb_file into the pose
    # pose_from_file(pose, pdb_filename)
    # pose_score = scorefxn(pose)
    # score_init_list.append(pose_score / pyrosetta.rosetta.core.pose.Pose.total_residue(pose))
    #
    # relax = FastRelax(standard_repeats=5)
    # relax.set_scorefxn(scorefxn)
    #
    # relax.constrain_relax_to_start_coords(True)
    # relax.coord_constrain_sidechains(True)
    # relax.ramp_down_constraints(False)

    for i in list(range(1, nb_of_repeats+1, 1)):
        # relax.apply(pose)
        # pose.dump_pdb("minimized_fast_CST_" + str(i) + "_" + pdb_filename)
        # pose_score_2 = scorefxn(pose)
        #
        # score_relax_dict[pdb_filename + "_" + str(i)] = pose_score_2 / pyrosetta.rosetta.core.pose.Pose.total_residue(pose)
        _main(pdb_filename)


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


# Global variables
pdbfile_list = []
score_init_list = []
score_relax_list = []
score_relax_dict = {}
nb_of_repeats = 1


def _main(start_pdb):
    out_pdb = start_pdb

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    n_start_poses = None

    if n_start_poses:
        n_start_poses = int(n_start_poses)
    else:
        n_start_poses = size

    rosetta_options = ["-ignore_unrecognized_res false",
                       "-ex1",
                       "-ex2",
                       "-use_input_sc",
                       "-flip_HNQ",
                       "-no_optH false",
                       "-relax:constrain_relax_to_start_coords",
                       "-relax:coord_constrain_sidechains",
                       "-relax:ramp_constraints false",
                       "-constant_seed",
                       "-no_his_his_pairE",
                       # "-linmem_ig 10",
                       "-nblist_autoupdate true",
                       # "-relax:sc_cst_maxdist 3",
                       # "-relax:coord_cst_width 1.0",
                       "-relax:coord_cst_stdev 0.5"]

    mpi_init(extra_options=" ".join(rosetta_options))

    n_local_jobs = None

    if n_start_poses > size:
        job_div = int(np.floor(n_start_poses / size))
        job_rem = n_start_poses - (job_div * size)
        job_div_procs = size - job_rem
        if rank <= job_div_procs:
            n_local_jobs = job_div
        else:
            n_local_jobs = job_div + 1
    else:
        job_div = 1
        if rank <= n_start_poses:
            n_local_jobs = 1
        else:
            n_local_jobs = 0

    for i in range(0, n_local_jobs):
        pack_id = "%d-%d" % (rank, n_local_jobs)
        outfile_name = out_pdb + ".pack" + pack_id + ".pdb"
        packer_job = FastRelaxPoseJob(start_pdb, rank, scorefn="ref2015_cst")
        packer_job.pack_pose()
        print("POSE %d in PROCESS %d COMPLETE, WRITING TO %s" % (i, rank, outfile_name))
        packer_job.dump_pose(outfile_name)

        score_relax_dict[start_pdb + "_" + str(i)] = packer_job.final_score / pyrosetta.rosetta.core.pose.Pose.total_residue(
            packer_job.pose)

if __name__ == '__main__':

    pdb_occupancy()
    import csv
    import numpy as np

    wtr = csv.writer(open('pyrosetta_out.csv', 'w+'), delimiter=',', lineterminator='\n', quoting=csv.QUOTE_NONE,
                     escapechar='\\')
    wtr.writerow(['pdb_filename', 'rmsd_init', 'score_init', 'rmsd_relax', 'score_relax'])

    padded_list = [0] * (nb_of_repeats * len(score_init_list))
    padded_list[::nb_of_repeats] = score_init_list
    padded_list_rmsd = [0] * (nb_of_repeats * len(score_init_list))

    l = [(x for x in score_relax_dict.keys()), (x for x in padded_list_rmsd), (x for x in padded_list),
         (x for x in padded_list_rmsd), (x for x in score_relax_dict.values())]
    wtr.writerows([i for i in zip(*l)])