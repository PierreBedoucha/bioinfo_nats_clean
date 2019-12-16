import os
import re
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns
import Bio.AlignIO as al


if __name__ == '__main__':

    prot_cpt_list = []
    prot_scores_list = []
    df = pd.DataFrame(columns=['pdb_filename', 'rmsd_init', 'score_init', 'rmsd_relax', 'score_relax'])

    for file in os.listdir("."):
        if 'minimized' in file:
            search = re.compile(r'^VRT_[0-9]+\s').search
            with open(file) as f1:
                for line in f1:
                    if search(line):
                        line_array = line.split(' ')
                        res_count = int(line_array[0].split('_')[1])
                        prot_scores_list.append(float(line_array[-1].rstrip()))
                        df = df.append({'pdb_filename':file, 'rmsd_init':0.0, 'score_init':0.0, 'rmsd_relax':0.0,
                                        'score_relax':float(line_array[-1].rstrip()) / res_count}, ignore_index=True)
                        # df['res_cpt'] = df['res_cpt'].astype(int)

    df.to_csv("pyrosetta_out_fromfiles.csv", sep=';', encoding='utf-8')
    # wtr.writerow(['pdb_filename', 'rmsd_init', 'score_init', 'rmsd_relax', 'score_relax'])