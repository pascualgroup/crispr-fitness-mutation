import pandas as pd
import numpy as np
import scipy as sp
import math
import matplotlib.pyplot as plt
from matplotlib import ticker
import sys
import os
import seaborn as sns
from scipy import stats
import sqlite3
from random import sample
from matplotlib.lines import Line2D

DBSIM_PATH = os.path.join('/Volumes/Yadgah/sylvain-martin-collab/12_MOI3/sweep_db_gathered.sqlite')
DBEXT_PATH = os.path.join('/Volumes/Yadgah/sylvain-martin-collab/12_MOI3/gathered-analyses/extinctions/extinctions.sqlite')
DBESC_PATH = os.path.join('/Volumes/Yadgah/sylvain-martin-collab/12_MOI3/gathered-analyses/escape-count/escape-count.sqlite')
####
conSim = sqlite3.connect(DBSIM_PATH)
curSim = conSim.cursor()
conExt = sqlite3.connect(DBEXT_PATH)
curExt = conExt.cursor()
conEsc = sqlite3.connect(DBESC_PATH)
curEsc = conEsc.cursor()
micMutProb = 0
viralMutProb = 1.04e-6
viralDecay = .1
failureProb = 0
carryCap = 316228.0
micDeath = 0.025
comboSpace = pd.read_sql_query(
"SELECT combo_id, evofunctionScale \
FROM param_combos WHERE microbe_mutation_prob = {0} AND viral_mutation_rate = {1} \
        AND viral_decay_rate = {2} AND crispr_failure_prob = {3} AND microbe_carrying_capacity = {4} \
        AND microbe_death_rate = {5} \
        ORDER BY combo_id" \
            .format(micMutProb,viralMutProb, viralDecay, failureProb, carryCap, micDeath), \
                conSim)
runIDs = pd.read_sql_query(
    "SELECT combo_id, run_id FROM runs \
        WHERE combo_id in ({})"\
        .format(', '.join(map(str, comboSpace['combo_id']))),conSim)


escapes = pd.DataFrame({'run_id':[], 'avg_escapes_per_vstrain':[], 'num_escaped_vstrains':[]})
for runID in runIDs['run_id']:
    subEscapes = pd.read_sql_query("SELECT DISTINCT vstrain_id, parent_vstrain_id, run_id FROM vstrain_escapes \
                            WHERE run_id in ({})"\
                            .format(runID), conEsc)
    if subEscapes.empty == False:
        subEscapes = subEscapes.groupby(['parent_vstrain_id','run_id'])\
                                .agg(num_escapes=('vstrain_id','size'))\
                                .reset_index()
        subEscapes = subEscapes.rename(columns={'parent_vstrain_id':'vstrain_id'})
        subEscapes = subEscapes.groupby(['run_id'])\
                            .agg(avg_escapes_per_vstrain=('num_escapes','mean'),\
                                    num_escaped_vstrains=('vstrain_id','size'))\
                            .reset_index()
        escapes = pd.concat([escapes,subEscapes], ignore_index=True)


escapes = escapes.merge(runIDs, on=['run_id']).merge(comboSpace, on=['combo_id'])\
    .sort_values(by=['evofunctionScale', 'run_id'])
meanEscapes = escapes.groupby(['combo_id', 'evofunctionScale'])\
    .agg(mean_escapes_per_vstrain=('avg_escapes_per_vstrain', 'mean'),
         std_escapes_per_vstrain=('avg_escapes_per_vstrain', 'std')).reset_index()

#
numExt = pd.read_sql_query("SELECT viruses, run_id FROM extinction_occurrence WHERE run_id in ({})"
                           .format(', '.join(map(str, runIDs['run_id']))), conExt)
numExt = numExt.merge(runIDs, on=['run_id']).merge(comboSpace, on=['combo_id'])\
    .sort_values(by=['evofunctionScale', 'run_id'])
meanNumExt = numExt.groupby(['combo_id', 'evofunctionScale'])\
    .agg(mean=('viruses', 'mean')).reset_index()
timeExt = pd.read_sql_query("SELECT virus_end_time, run_id FROM simulation_end_time WHERE run_id in ({})"
                           .format(', '.join(map(str, runIDs['run_id']))), conExt)
timeExt = timeExt.merge(numExt[numExt.viruses==1],on=['run_id'])
timeExt = timeExt.drop(columns=['viruses'])
timeExt = timeExt.rename(columns={'virus_end_time':'t_extinction'})
count = timeExt.groupby(['combo_id','evofunctionScale'])\
            .agg(n=('run_id', 'size')).reset_index()
maxRuns = 300
maxT = 2000
for comboID in np.unique(runIDs['combo_id']):
    if len(count[count.combo_id==comboID][count.n<maxRuns]) > 0:
        timeExt = pd.concat([timeExt, pd.DataFrame(
                    {'t_extinction': len(count[count.n<maxRuns])*[maxT],
                    'run_id': len(count[count.n<maxRuns])*[maxRuns+1],
                    'combo_id': count[count.n<maxRuns]['combo_id'],
                    'evofunctionScale': count[count.n<maxRuns]['evofunctionScale']
                    })],
                    ignore_index=True)
meanTimeExt = timeExt.groupby(['combo_id', 'evofunctionScale'])\
    .agg(mean=('t_extinction', 'mean'),std=('t_extinction', 'std')).reset_index()
numRuns = timeExt.groupby(['combo_id', 'evofunctionScale'])\
    .agg(n=('t_extinction', 'size')).reset_index()
meanTimeExt  = meanTimeExt.merge(numRuns,on=['combo_id','evofunctionScale'])


# fig, ax1 = plt.subplots(1,2 figsize=(10,10))
# ax = [ax1[0], ax[]]
fig = plt.figure(figsize=(20,8))
gs = fig.add_gridspec(1, 2, hspace=0, wspace=.35)
(ax1, ax2) = gs.subplots()
ax = [ax1, ax2]
# ax[1].fill_between(modMeanTimeExt['evofunctionScale'],
#                 modMeanTimeExt['mean'] -
#                 modMeanTimeExt['std'],
#                 modMeanTimeExt['mean'] +
#                 modMeanTimeExt['std'], color='darkseagreen', alpha=0.25)
# ax[0].fill_between(modMeanNumEsc['evofunctionScale'],
#                 modMeanNumEsc['meanProp'] -
#                 modMeanNumEsc['stdProp'],
#                 modMeanNumEsc['meanProp'] +
#                 modMeanNumEsc['stdProp'], color='darkslateblue', alpha=0.25)
# ax[0].scatter(modMeanNumEsc['evofunctionScale'], modMeanNumEsc['meanProp'], color='darkslateblue',
#           s=30)
# ax[0].plot(modMeanNumEsc['evofunctionScale'], modMeanNumEsc['meanProp'], color='darkslateblue',
#           linewidth=2,alpha=0.1)
ax[0].errorbar(meanEscapes['evofunctionScale'], meanEscapes['mean_escapes_per_vstrain'], 
               yerr=meanEscapes['std_escapes_per_vstrain'], fmt='o',
               color = 'darkslateblue',capsize=2)
# ax[0].scatter(modMeanNumEsc['evofunctionScale'], modMeanNumEsc['meanProp'], color='darkslateblue',
#           s=.2)
# ax[0].set_xticks(modMeanNumEsc['evofunctionScale'])
# ax[0].set_xticklabels(modMeanNumEsc['evofunctionScale'])
ax[0].set_ylim(0, 100)
# ax[0].set_yticks([np.round(i,4) for i in np.arange(.4,1.01,.1)])
# ax[1].scatter(modMeanTimeExt['evofunctionScale'], modMeanTimeExt['mean'], color='darkseagreen',
#           s=30)
# ax[1].plot(modMeanTimeExt['evofunctionScale'], modMeanTimeExt['mean'], color='darkseagreen',
#           linewidth=2,alpha=0.1)
ax[1].errorbar(meanTimeExt['evofunctionScale'], meanTimeExt['mean'], 
               yerr=meanTimeExt['std'], fmt='o',
               color = 'darkseagreen',capsize=2)
# ax[1].set_xticks(modMeanTimeExt['evofunctionScale'])
# ax[1].set_xticklabels(modMeanTimeExt['evofunctionScale'])
lim = ax[1].get_ylim()
ax[1].set_ylim(100, 901)
ax[1].set_yticks([np.round(i,4) for i in np.arange(100,901,200)])
# ax[1].set_yticklabels(modMeanTimeExt['evofunctionScale'])
s = ['Expected Proportion of\n',r'Host Strains Escaped $p_{esc}$']
ax[0].set_ylabel(ylabel=''.join(s), labelpad=15, fontsize=15)
s = ['Expected Time to\n',r'Viral Extinction $\tau_{ext}$']
# ax[1].set_ylabel(ylabel=r''.join(s),
#                  labelpad=45, fontsize=15,rotation=270)
ax[1].set_ylabel(ylabel=r''.join(s),
                 labelpad=15, fontsize=15)
ax[0].set_xlabel(xlabel=r'Intensity of Selection $\sigma$', labelpad=15, fontsize=15)
ax[1].set_xlabel(xlabel=r'Intensity of Selection $\sigma$', labelpad=15, fontsize=15)
# ax[0].set_axisbelow(True)
# ax[0].yaxis.grid(color='gray',linestyle='dashed')
# ax[1].set_axisbelow(True)
# ax[1].yaxis.grid(color='gray',linestyle='dashed')
# fig.tight_layout()
# ax[0].legend(loc="lower right", fontsize=15,
#     handles=[
#         Line2D([], [], c='darkslateblue', label=r"$p_{esc}$"),
#         Line2D([], [], c='darkseagreen', label=r"$\tau_{ext}$")])
ax[0].tick_params(axis='x', labelsize=15)
ax[0].tick_params(axis='y', labelsize=15)
ax[1].tick_params(axis='x', labelsize=15)
ax[1].tick_params(axis='y', labelsize=15)
# fig.tight_layout()
plt.show()






conSim = sqlite3.connect(DBSIM_PATH)
curSim = conSim.cursor()
conExt = sqlite3.connect(DBEXT_PATH)
curExt = conExt.cursor()
conEsc = sqlite3.connect(DBESC_PATH)
curEsc = conEsc.cursor()
micMutProb = 0
viralMutProb = 1.04e-6
viralDecay = .1
failureProb = 0
carryCap = 316228.0
micDeath = 0.025
comboSpace = pd.read_sql_query(
"SELECT combo_id, evofunctionScale \
FROM param_combos WHERE microbe_mutation_prob = {0} AND viral_mutation_rate = {1} \
        AND viral_decay_rate = {2} AND crispr_failure_prob = {3} AND microbe_carrying_capacity = {4} \
        AND microbe_death_rate = {5} \
        ORDER BY combo_id" \
            .format(micMutProb,viralMutProb, viralDecay, failureProb, carryCap, micDeath), \
                conSim)
runIDs = pd.read_sql_query(
    "SELECT combo_id, run_id FROM runs \
        WHERE combo_id in ({})"\
        .format(', '.join(map(str, comboSpace['combo_id']))),conSim)


escapes = pd.DataFrame({'run_id':[], 'num_unique_escaped_vstrains':[]})
for runID in runIDs['run_id']:
    subEscapes = pd.read_sql_query("SELECT DISTINCT parent_vstrain_id, run_id FROM vstrain_escapes \
                            WHERE run_id in ({})"\
                            .format(runID), conEsc)
    if subEscapes.empty == False:
        subEscapes = subEscapes.groupby(['run_id'])\
                                .agg(num_unique_escaped_vstrains=('parent_vstrain_id','size'))\
                                .reset_index()
        escapes = pd.concat([escapes,subEscapes], ignore_index=True)


escapes = escapes.merge(runIDs, on=['run_id']).merge(comboSpace, on=['combo_id'])\
    .sort_values(by=['evofunctionScale', 'run_id'])
meanEscapes = escapes.groupby(['combo_id', 'evofunctionScale'])\
    .agg(mean_num_unique_escaped_vstrains=('num_unique_escaped_vstrains', 'mean'),
         std_num_unique_escaped_vstrains=('num_unique_escaped_vstrains', 'std')).reset_index()

#
numExt = pd.read_sql_query("SELECT viruses, run_id FROM extinction_occurrence WHERE run_id in ({})"
                           .format(', '.join(map(str, runIDs['run_id']))), conExt)
numExt = numExt.merge(runIDs, on=['run_id']).merge(comboSpace, on=['combo_id'])\
    .sort_values(by=['evofunctionScale', 'run_id'])
meanNumExt = numExt.groupby(['combo_id', 'evofunctionScale'])\
    .agg(mean=('viruses', 'mean')).reset_index()
timeExt = pd.read_sql_query("SELECT virus_end_time, run_id FROM simulation_end_time WHERE run_id in ({})"
                           .format(', '.join(map(str, runIDs['run_id']))), conExt)
timeExt = timeExt.merge(numExt[numExt.viruses==1],on=['run_id'])
timeExt = timeExt.drop(columns=['viruses'])
timeExt = timeExt.rename(columns={'virus_end_time':'t_extinction'})
count = timeExt.groupby(['combo_id','evofunctionScale'])\
            .agg(n=('run_id', 'size')).reset_index()
maxRuns = 300
maxT = 2000
for comboID in np.unique(runIDs['combo_id']):
    if len(count[count.combo_id==comboID][count.n<maxRuns]) > 0:
        timeExt = pd.concat([timeExt, pd.DataFrame(
                    {'t_extinction': len(count[count.n<maxRuns])*[maxT],
                    'run_id': len(count[count.n<maxRuns])*[maxRuns+1],
                    'combo_id': count[count.n<maxRuns]['combo_id'],
                    'evofunctionScale': count[count.n<maxRuns]['evofunctionScale']
                    })],
                    ignore_index=True)
meanTimeExt = timeExt.groupby(['combo_id', 'evofunctionScale'])\
    .agg(mean=('t_extinction', 'mean'),std=('t_extinction', 'std')).reset_index()
numRuns = timeExt.groupby(['combo_id', 'evofunctionScale'])\
    .agg(n=('t_extinction', 'size')).reset_index()
meanTimeExt  = meanTimeExt.merge(numRuns,on=['combo_id','evofunctionScale'])


# fig, ax1 = plt.subplots(1,2 figsize=(10,10))
# ax = [ax1[0], ax[]]
fig = plt.figure(figsize=(20,8))
gs = fig.add_gridspec(1, 2, hspace=0, wspace=.35)
(ax1, ax2) = gs.subplots()
ax = [ax1, ax2]
# ax[1].fill_between(modMeanTimeExt['evofunctionScale'],
#                 modMeanTimeExt['mean'] -
#                 modMeanTimeExt['std'],
#                 modMeanTimeExt['mean'] +
#                 modMeanTimeExt['std'], color='darkseagreen', alpha=0.25)
# ax[0].fill_between(modMeanNumEsc['evofunctionScale'],
#                 modMeanNumEsc['meanProp'] -
#                 modMeanNumEsc['stdProp'],
#                 modMeanNumEsc['meanProp'] +
#                 modMeanNumEsc['stdProp'], color='darkslateblue', alpha=0.25)
# ax[0].scatter(modMeanNumEsc['evofunctionScale'], modMeanNumEsc['meanProp'], color='darkslateblue',
#           s=30)
# ax[0].plot(modMeanNumEsc['evofunctionScale'], modMeanNumEsc['meanProp'], color='darkslateblue',
#           linewidth=2,alpha=0.1)
ax[0].errorbar(meanEscapes['evofunctionScale'], meanEscapes['mean_num_unique_escaped_vstrains'], 
               yerr=meanEscapes['std_num_unique_escaped_vstrains'], fmt='o',
               color = 'darkslateblue',capsize=2)
# ax[0].scatter(modMeanNumEsc['evofunctionScale'], modMeanNumEsc['meanProp'], color='darkslateblue',
#           s=.2)
# ax[0].set_xticks(modMeanNumEsc['evofunctionScale'])
# ax[0].set_xticklabels(modMeanNumEsc['evofunctionScale'])
ax[0].set_ylim(0, 100)
# ax[0].set_yticks([np.round(i,4) for i in np.arange(.4,1.01,.1)])
# ax[1].scatter(modMeanTimeExt['evofunctionScale'], modMeanTimeExt['mean'], color='darkseagreen',
#           s=30)
# ax[1].plot(modMeanTimeExt['evofunctionScale'], modMeanTimeExt['mean'], color='darkseagreen',
#           linewidth=2,alpha=0.1)
ax[1].errorbar(meanTimeExt['evofunctionScale'], meanTimeExt['mean'], 
               yerr=meanTimeExt['std'], fmt='o',
               color = 'darkseagreen',capsize=2)
# ax[1].set_xticks(modMeanTimeExt['evofunctionScale'])
# ax[1].set_xticklabels(modMeanTimeExt['evofunctionScale'])
lim = ax[1].get_ylim()
ax[1].set_ylim(100, 901)
ax[1].set_yticks([np.round(i,4) for i in np.arange(100,901,200)])
# ax[1].set_yticklabels(modMeanTimeExt['evofunctionScale'])
s = ['Expected Proportion of\n',r'Host Strains Escaped $p_{esc}$']
ax[0].set_ylabel(ylabel=''.join(s), labelpad=15, fontsize=15)
s = ['Expected Time to\n',r'Viral Extinction $\tau_{ext}$']
# ax[1].set_ylabel(ylabel=r''.join(s),
#                  labelpad=45, fontsize=15,rotation=270)
ax[1].set_ylabel(ylabel=r''.join(s),
                 labelpad=15, fontsize=15)
ax[0].set_xlabel(xlabel=r'Intensity of Selection $\sigma$', labelpad=15, fontsize=15)
ax[1].set_xlabel(xlabel=r'Intensity of Selection $\sigma$', labelpad=15, fontsize=15)
# ax[0].set_axisbelow(True)
# ax[0].yaxis.grid(color='gray',linestyle='dashed')
# ax[1].set_axisbelow(True)
# ax[1].yaxis.grid(color='gray',linestyle='dashed')
# fig.tight_layout()
# ax[0].legend(loc="lower right", fontsize=15,
#     handles=[
#         Line2D([], [], c='darkslateblue', label=r"$p_{esc}$"),
#         Line2D([], [], c='darkseagreen', label=r"$\tau_{ext}$")])
ax[0].tick_params(axis='x', labelsize=15)
ax[0].tick_params(axis='y', labelsize=15)
ax[1].tick_params(axis='x', labelsize=15)
ax[1].tick_params(axis='y', labelsize=15)
# fig.tight_layout()
plt.show()

