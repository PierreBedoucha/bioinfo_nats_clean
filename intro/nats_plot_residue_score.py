import os
import re
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns
import Bio.AlignIO as al

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


ca_align_dict = read_msa_fasta()


if __name__ == '__main__':

    res_cpt_list = []
    res_scores_list = []
    df = pd.DataFrame(columns=['res_cpt', 'res_scores'])

    for file in os.listdir("."):
        if 'minimized' in file:
            res_to_be_aligned = ca_align_dict[file.split("_")[-8]]
            search = re.compile(r'[^A-Z]*_[0-9]+\s').search
            res_cpt = 1
            with open(file) as f1:
                for line in f1:
                    if search(line):
                        line_array = line.split(' ')
                        res_cpt_list.append(res_cpt)
                        res_scores_list.append(float(line_array[-1].rstrip()))
                        df = df.append({'res_cpt': int(res_cpt), 'res_scores': float(line_array[-1].rstrip())}, ignore_index=True)
                        df['res_cpt'] = df['res_cpt'].astype(int)

                        res_cpt += 1
            # df = df.assign(res_cpt=pd.Series(res_cpt_list))
            # df = df.assign(res_scores=pd.Series(res_scores_list))

    ax = sns.boxplot(x="res_cpt", y="res_scores", data=df)
    ax.get_xaxis().set_major_formatter(
        ticker.FuncFormatter(lambda x, p: format(int(x))))
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45)
    ind_list = list(range(0, 200, 10))
    for ind, label in enumerate(ax.get_xticklabels()):
        if int(label._text) in ind_list:
            label.set_visible(True)
        else:
            label.set_visible(False)

    plt.show()


