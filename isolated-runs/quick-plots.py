#!/usr/bin/env python3

import pandas as pd
import numpy as np
import scipy as sp
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import ticker
import sys
import os
import seaborn as sns
from scipy import stats
import sqlite3
from mpl_toolkits.axes_grid1 import make_axes_locatable

dir = 'sylvain-martin-collab/trial-3/gathered-analyses'
DBWALLS_PATH = os.path.join(
    '/Volumes/Yadgah/{0}/walls-shannon/mean-walls-shannon.sqlite'.format(dir))
DBEXT_PATH = os.path.join(
    '/Volumes/Yadgah/{0}/extinctions/mean-extinctions.sqlite'.format(dir))


DBSIM_PATH = os.path.join(
    '/Volumes/Yadgah/sylvain-martin-collab/8_MOI3/sweep_db.sqlite')
DBEXT_PATH = os.path.join(
    '/Volumes/Yadgah/sylvain-martin-collab/8_MOI3/gathered-analyses/extinctions/extinctions.sqlite')
DBESC_PATH = os.path.join(
    '/Volumes/Yadgah/sylvain-martin-collab/8_MOI3/gathered-analyses/escape-count/escape-count.sqlite')
conSim = sqlite3.connect(DBSIM_PATH)
curSim = conSim.cursor()
conExt = sqlite3.connect(DBEXT_PATH)
curExt = conExt.cursor()
conEsc = sqlite3.connect(DBESC_PATH)
curEsc = conEsc.cursor()

micMutProb = 1
viralMutProb = 1.04e-6
viralDecay = .05
failureProb = 0
carryCap = 316228.0
micDeath = .025
comboSpace = pd.read_sql_query(
    "SELECT combo_id, max_allele \
FROM param_combos WHERE microbe_mutation_prob = {0} AND viral_mutation_rate = {1} \
        AND viral_decay_rate = {2} AND crispr_failure_prob = {3} AND microbe_carrying_capacity = {4} \
        AND microbe_death_rate = {5} \
ORDER BY combo_id".format(micMutProb, viralMutProb, viralDecay, failureProb, carryCap, micDeath), conSim)
runIDs = pd.read_sql_query(
    "SELECT combo_id, run_id FROM runs \
        WHERE combo_id in ({})"
    .format(', '.join(map(str, comboSpace['combo_id']))), conSim)
maxAlleles = [i for (i,) in conSim.execute("SELECT DISTINCT max_allele \
                                           FROM param_combos\
                                           ORDER BY max_allele").fetchall()]


print('SQLite Query: mean microbial wall occurences per parameter combination')
numEsc = pd.read_sql_query("SELECT num_escapes, run_id FROM escape_count WHERE run_id in ({})"
                           .format(', '.join(map(str, runIDs['run_id']))), conEsc)
numEsc = numEsc.merge(runIDs, on=['run_id']).merge(comboSpace, on=['combo_id'])\
    .sort_values(by=['max_allele', 'run_id'])
meanNumEsc = numEsc.groupby(['combo_id', 'max_allele'])\
    .agg(mean=('num_escapes', 'sum')).reset_index()
numRuns = numEsc.groupby(['combo_id', 'max_allele'])\
    .agg(n=('num_escapes', 'size')).reset_index()
meanNumEsc = meanNumEsc.merge(numRuns, on=['combo_id', 'max_allele'])
meanNumEsc['mean'] = np.array(meanNumEsc['mean'])/np.array(meanNumEsc['n'])


print('SQLite Query: mean microbial wall occurences per parameter combination')
numExt = pd.read_sql_query("SELECT viruses, run_id FROM extinction_occurrence WHERE run_id in ({})"
                           .format(', '.join(map(str, runIDs['run_id']))), conExt)
numExt = numExt.merge(runIDs, on=['run_id']).merge(comboSpace, on=['combo_id'])\
    .sort_values(by=['max_allele', 'run_id'])
meanNumExt = numExt.groupby(['combo_id', 'max_allele'])\
    .agg(mean=('viruses', 'mean')).reset_index()
# numRuns = numExt.groupby(['combo_id', 'max_allele'])\
#     .agg(n=('viruses', 'size')).reset_index()
# meanNumExt = meanNumExt.merge(numRuns, on=['combo_id', 'max_allele'])
# meanNumExt['mean'] = np.array(meanNumExt['mean'])/np.array(meanNumExt['n'])

timeExt = pd.read_sql_query("SELECT virus_end_time, run_id FROM simulation_end_time WHERE run_id in ({})"
                            .format(', '.join(map(str, runIDs['run_id']))), conExt)
timeExt = timeExt.merge(numExt[numExt.viruses == 1], on=['run_id'])
meanTimeExt = timeExt.groupby(['combo_id', 'max_allele'])\
    .agg(mean=('virus_end_time', 'mean')).reset_index()
numRuns = timeExt.groupby(['combo_id', 'max_allele'])\
    .agg(n=('virus_end_time', 'size')).reset_index()
# meanTimeExt = meanTimeExt.merge(numRuns, on=['combo_id', 'max_allele'])
# meanTimeExt['mean'] = np.array(meanTimeExt['mean'])/np.array(meanTimeExt['n'])

# meanTimeExt.loc[len(meanTimeExt.index)] = [37, 135, 0, 50]
# check = pd.read_sql_query("SELECT microbes, run_id FROM extinction_occurrence WHERE run_id in ({})"
#                            .format(', '.join(map(str, runIDs['run_id']))), conExt)

fig, ax = plt.subplots(3, sharex=True, figsize=(10, 7.5))
ax[0].bar(meanNumEsc['max_allele'], meanNumEsc['mean'], color='darkred',
          width=10)
ax[0].set_xticks(maxAlleles)
ax[0].set_xticklabels(maxAlleles)
ax[1].bar(meanNumExt['max_allele'], meanNumExt['mean'], color='darkgreen',
          width=10)
ax[1].set_xticks(maxAlleles)
ax[1].set_xticklabels(maxAlleles)
ax[2].bar(meanTimeExt['max_allele'], meanTimeExt['mean'], color='darkblue',
          width=10)
ax[2].set_xticks(maxAlleles)
ax[2].set_xticklabels(maxAlleles)
ax[0].set_ylabel(ylabel='Expected Number of\nEscapes',
                 labelpad=15, fontsize=10)
ax[1].set_ylabel(ylabel='Probability of\nExtinction', labelpad=15, fontsize=10)
ax[2].set_ylabel(ylabel='Expected Time to\nExtinction',
                 labelpad=15, fontsize=10)
ax[2].set_xlabel(xlabel='Max Allele', fontsize=10)
fig.suptitle('Expectations @ host fit. mut = {0}, viral mut = {1},\nviral decay = {2}, host mort = {3},'.format(micMutProb, viralMutProb,
                                                                                                                viralDecay, micDeath), fontsize=20)


wallOccurrences = pd.read_sql_query(
    "SELECT wall_occurrence, combo_id FROM microbial_wall_statistics ORDER BY combo_id", conWalls)
wallOccurrences = wallOccurrences.merge(comboSpace, on=['combo_id']).drop(columns=['combo_id'])\
    .sort_values(by=['max_allele'])
numWalls = pd.read_sql_query(
    "SELECT expected_num_walls, combo_id FROM microbial_wall_statistics ORDER BY combo_id", conWalls)
numWalls = numWalls.merge(comboSpace, on=['combo_id']).drop(columns=['combo_id'])\
    .sort_values(by=['max_allele'])
conExt = sqlite3.connect(DBEXT_PATH)
curExt = conExt.cursor()
extOccurrences = pd.read_sql_query(
    "SELECT microbes,viruses,combo_id FROM mean_extinction_occurrences ORDER BY combo_id", conExt)
extOccurrences = extOccurrences.merge(comboSpace, on=['combo_id']).drop(columns=['combo_id'])\
    .sort_values(by=['max_allele'])
simEndTimes = pd.read_sql_query(
    "SELECT microbe_end_time,virus_end_time,combo_id FROM mean_simulation_end_times ORDER BY combo_id", conExt)
simEndTimes = simEndTimes .merge(comboSpace, on=['combo_id']).drop(columns=['combo_id'])\
    .sort_values(by=['max_allele'])

fig, ax = plt.subplots(4, sharex=True)
ax[0].bar(wallOccurrences['max_allele'], wallOccurrences['wall_occurrence'], color='darkred',
          width=2)
ax[0].set_xticks([30, 35, 45, 60, 80, 105, 135])
ax[0].set_xticklabels([30, 35, 45, 60, 80, 105, 135])
ax[1].bar(numWalls['max_allele'], numWalls['expected_num_walls'], color='darkorange',
          width=2)
ax[1].set_xticks([30, 35, 45, 60, 80, 105, 135])
ax[1].set_xticklabels([30, 35, 45, 60, 80, 105, 135])
ax[2].bar(extOccurrences['max_allele'], extOccurrences['viruses'], color='darkgreen',
          width=2)
ax[2].set_xticks([30, 35, 45, 60, 80, 105, 135])
ax[2].set_xticklabels([30, 35, 45, 60, 80, 105, 135])
ax[3].bar(simEndTimes['max_allele'], simEndTimes['virus_end_time'], color='darkblue',
          width=2)
ax[3].set_xticks([30, 35, 45, 60, 80, 105, 135])
ax[3].set_xticklabels([30, 35, 45, 60, 80, 105, 135])
ax[0].set_ylabel(ylabel='Expected Wall Occurrence', labelpad=15, fontsize=7)
ax[1].set_ylabel(ylabel='Expected Number of Walls', labelpad=15, fontsize=7)
ax[2].set_ylabel(ylabel='Expected Extinction Occurrence',
                 labelpad=15, fontsize=7)
ax[3].set_ylabel(ylabel='Expected Time to Viral Extinction',
                 labelpad=15, fontsize=7)
ax[3].set_xlabel(xlabel='Max Allele', fontsize=7)


def fgm(a, centerA, maxA, maxFit):
    minA = 2*centerA-maxA
    scale = -1*maxFit/((centerA-minA)*(centerA - maxA))
    return -1*(a-minA)*(a-maxA)*scale


a = list(np.linspace(-30, 30, 100))
centerA = 0
maxFit = 1
for maxA in maxAlleles:
    fig, ax = plt.subplots(1, sharex=True)
    ax.plot(a, list(map(lambda x: fgm(x, centerA, maxA, maxFit), a)), label='max allele: {0}'.format(maxA),
            color='darkblue', linewidth=1.5)
    ax.set_ylim(0, 1)
    ax.legend()
    fig.savefig(os.path.join(
        '/Users/armun/Desktop/fgm{0}.png'.format(maxA)), dpi=500)

# fig, ax = plt.subplots(1,sharex=True)
# ax.plot([0,1e-5,1e-3,1e-2],[7/30,7/30,8/30,11/30],color='darkred',linewidth=2)
# ax.set_xticks([0,1e-5,1e-3,1e-2])
# ax.set_xticklabels([0,1e-5,1e-3,1e-2])
# ax.set_ylabel(ylabel ='Proportion of Realizations',labelpad=15,fontsize=7)
#
# fig, ax = plt.subplots(1,sharex=True)
# ax.plot([0,1e-5,1e-3,1e-2],[5/30,6/30,4/30,1/30],color='darkblue',linewidth=2)
# ax.set_xticks([0,1e-5,1e-3,1e-2])
# ax.set_xticklabels([0,1e-5,1e-3,1e-2])
# ax.set_ylabel(ylabel ='Proportion of Realizations',labelpad=15,fontsize=7)
