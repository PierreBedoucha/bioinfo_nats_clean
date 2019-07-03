import pandas as pd
import seaborn as sns
from tempfile import mkstemp
from shutil import move
from os import fdopen, remove
# import matplotlib.pyplot as plt

# train = pd.read_csv('out.csv')
#
# sns.barplot(x="pdbid", y="rmsd", data=train)
# plt.title('Ca RMSD in A')
# plt.show()

from pyrosetta import *
from rosetta.core.scoring import *
from rosetta.protocols.relax import *

def score_proteins(pdb_filename):
    init()
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
    # print(scorefxn)
    scorefxn = pyrosetta.rosetta.core.scoring.ScoreFunctionFactory.create_score_function('talaris2014_cst')
    scorefxn.show()


    # create a pose from the desired PDB file
    # -create an empty Pose object
    pose = Pose()
    # -load the data from pdb_file into the pose
    pose_from_file(pose, pdb_filename)
    # default to the median residue number
    residues = range(1, pose.total_residue() + 1)
    pose_score = scorefxn(pose)
    # print(pose_score)
    score_init_list.append(pose_score)

    # === Bolean for Relax here ===
    # relax = ClassicRelax()
    # relax.set_scorefxn(scorefxn)
    # relax.apply(pose)
    # =============================

    relax = FastRelax()
    relax.set_scorefxn(scorefxn)
    relax.apply(pose)

    pose.dump_pdb("minimized_fast" + pdb_filename)
    pose_score_2 = scorefxn(pose)
    # print(pose_score_2)

    score_relax_list.append(pose_score_2)

def replace(file_path, pattern, subst):
    #Create temp file
    fh, abs_path = mkstemp()
    with fdopen(fh,'w') as new_file:
        with open(file_path) as old_file:
            for line in old_file:
                new_file.write(line.replace(pattern, subst))
    #Remove original file
    remove(file_path)
    #Move new file
    move(abs_path, file_path)

def pdb_occupancy():
    import os
    for file in os.listdir("."):
        if file.endswith(".pdb") and 'minimized' not in file:
            fh, abs_path = mkstemp()
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

    l = [(x for x in pdbfile_list), (x for x in score_init_list), (x for x in score_relax_list)]
    wtr.writerows([i for i in zip(*l)])