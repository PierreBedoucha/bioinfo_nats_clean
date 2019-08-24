

import pandas as pd
import seaborn as sns
from tempfile import mkstemp
from shutil import move
import os
import matplotlib.pyplot as plt
import math
from scipy import stats
from scipy.stats import *
import numpy as np



if __name__ == '__main__':
    train = pd.read_csv('Workbook17.csv')

    new = train["pdb_filename"].str.split("_", n=8, expand=True)
    train['amplitude'] = new[6]
    train['amplitude'] = train['amplitude'].astype(float)

    train['repeat'] = new[8]

    train['pdbid'] = train.pdb_filename.str[:4]
    grouped = train.groupby(["pdbid"])
    train = grouped.apply(lambda x: x.sort_values(["amplitude"], ascending=True)).reset_index(drop=True)

    # amplitude_list = [0.00, 3.11, 7.33, 11.56, 15.78, 20.00]
    # for j in amplitude_list:
    #     iter_list = []
    #     score_list = train.loc[(train.amplitude == j) & (train['repeat'].notnull()), 'score_relax']
    #     for i in score_list:
    #         iter_list.append(i)
    #         if len(iter_list) > 1:
    #             db = iter_list[:-1]
    #         else:
    #             db = [iter_list[0]]
    #         # ttest_sp, p = ttest_ind(db, iter_list)
    #         mean_test = np.mean(db)
    #         ttest, p = ttest_1samp(iter_list, mean_test)
    #         # n = len(db) + len(da)
    #         # df = n - 2
    #         print(1 - p)


    # j = 15.78
    # iter_list = []
    # score_list = train.loc[(train.amplitude == j) & (train['repeat'].notnull()), 'score_relax']
    # for i in score_list:
    #     iter_list.append(i)
    #     if len(iter_list) > 1:
    #         db = iter_list[:-1]
    #         # -----------
    #         iter_list.pop(0)
    #     else:
    #         db = [iter_list[0]]
    #     # ttest_sp, p = ttest_ind(db, iter_list)
    #     mean_test = np.mean(db)
    #     # ttest, p = ttest_1samp(iter_list, mean_test)
    #     ttest, p = ttest_rel(iter_list, db)
    #     # n = len(db) + len(da)
    #     # df = n - 2
    #     print(p)


    # # score_list = np.random.randn(500) + 1
    # mu = 1
    # sigma = 9
    # score_list = np.random.normal(mu, sigma, 500)
    # iter_list = []
    # db = []
    # for i in score_list:
    #     iter_list.append(i)
    #     if len(iter_list) > 1:
    #         db = iter_list[:-1]
    #     else:
    #         db = [iter_list[0]]
    #     mean_test = np.mean(db)
    #     ttest, p = ttest_1samp(iter_list, 0)
    #     # ttest_sp, p = ttest_rel(db, iter_list)
    #     # ttest_sp, p = wilcoxon(db, iter_list)
    #     # n = len(db) + len(da)
    #     # df = n - 2
    #     print(p)

    # j = 15.78
    # iter_list = []
    # score_list = train.loc[(train.amplitude == j) & (train['repeat'].notnull()), 'score_relax']
    # for i in score_list:
    #     iter_list.append(i)

    # train.to_csv("stats_out.csv", sep=';', encoding='utf-8')

    # sns.distplot(iter_list)

    # sns.catplot(x="amplitude", y="score_relax", kind="swarm", data=train)
    # sns.catplot(x="amplitude", y="score_relax", jitter=False, data=train)

    ax = sns.boxplot(x="amplitude", y="score_relax", data=train)

    plt.show()