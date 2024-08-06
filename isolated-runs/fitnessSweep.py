import pandas as pd
import numpy as np
import scipy as sp
import math
import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
from matplotlib import ticker
import sys
import os
import seaborn as sns
from scipy import stats
import sqlite3
from random import sample
from matplotlib.lines import Line2D

DBSIM_PATH = os.path.join('/Volumes/Yadgah/sylvain-martin-collab/12_MOI3/sweep_db_gathered.sqlite')
# DBSIM_PATH = os.path.join('/Volumes/Yadgah/sylvain-martin-collab/12_MOI3/sweep_db.sqlite')
DBEXT_PATH = os.path.join('/Volumes/Yadgah/sylvain-martin-collab/12_MOI3/gathered-analyses/extinctions/extinctions.sqlite')
DBESC_PATH = os.path.join('/Volumes/Yadgah/sylvain-martin-collab/12_MOI3/gathered-analyses/escape-count/escape-count.sqlite')
DBSIM_PATHnull = os.path.join('/Volumes/Yadgah/sylvain-martin-collab/14_MOI3/sweep_db_gathered.sqlite')
DBEXT_PATHnull = os.path.join('/Volumes/Yadgah/sylvain-martin-collab/14_MOI3/gathered-analyses/extinctions/extinctions.sqlite')
DBESC_PATHnull = os.path.join('/Volumes/Yadgah/sylvain-martin-collab/14_MOI3/gathered-analyses/escape-count/escape-count.sqlite')
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




numEsc = []
for runID in runIDs['run_id']:
    # print(runID)
    numEsc += [len(pd.read_sql_query("SELECT DISTINCT bstrain_id_escaped, run_id FROM vstrain_escapes \
                           WHERE run_id in ({}) AND t_creation <= 2000"\
                           .format(runID), conEsc)['bstrain_id_escaped'].values)]
    
numEsc = pd.DataFrame({'num_escapes':numEsc,'run_id':runIDs['run_id']})
div = []
for runID in runIDs['run_id']:
    # print(runID)
    div += [len(pd.read_sql_query("SELECT bstrain_id FROM bstrains WHERE run_id in ({}) \
                                AND t_creation <= 2000 ORDER BY bstrain_id"\
                            .format(runID), conSim)['bstrain_id'].values)]
div = pd.DataFrame({'div':div,'run_id':runIDs['run_id']})
numEsc = numEsc.merge(div,on=['run_id'])
numEsc['prop'] = numEsc['num_escapes']/numEsc['div']
numEsc = numEsc.merge(runIDs,on=['run_id']).merge(comboSpace,on=['combo_id'])\
                .sort_values(by=['evofunctionScale','run_id'])
meanNumEsc = numEsc.groupby(['combo_id','evofunctionScale'])\
                .agg(meanEsc=('num_escapes','mean'),meanDiv=('div','mean'), 
                     meanProp=('prop','mean'),stdProp=('prop','std')).reset_index()
numRuns = numEsc.groupby(['combo_id', 'evofunctionScale'])\
    .agg(n=('num_escapes', 'size')).reset_index()
meanNumEsc = meanNumEsc.merge(numRuns,on=['combo_id','evofunctionScale'])
# meanNumEsc['mean'] = np.array(meanNumEsc['mean'])/np.array(meanNumEsc['n'])
#
numExt = pd.read_sql_query("SELECT viruses, run_id FROM extinction_occurrence WHERE run_id in ({})"
                           .format(', '.join(map(str, runIDs['run_id']))), conExt)
numExt = numExt.merge(runIDs, on=['run_id']).merge(comboSpace, on=['combo_id'])\
    .sort_values(by=['evofunctionScale', 'run_id'])
meanNumExt = numExt.groupby(['combo_id', 'evofunctionScale'])\
    .agg(mean=('viruses', 'mean')).reset_index()
# numRuns = numExt.groupby(['combo_id', 'evofunctionScale'])\
#     .agg(n=('viruses', 'size')).reset_index()
# meanNumExt = meanNumExt.merge(numRuns, on=['combo_id', 'evofunctionScale'])
# meanNumExt['mean'] = np.array(meanNumExt['mean'])/np.array(meanNumExt['n'])

timeExt = pd.read_sql_query("SELECT virus_end_time, run_id FROM simulation_end_time WHERE run_id in ({})"
                           .format(', '.join(map(str, runIDs['run_id']))), conExt)
timeExt = timeExt.merge(numExt[numExt.viruses==1],on=['run_id'])
meanTimeExt = timeExt.groupby(['combo_id', 'evofunctionScale'])\
    .agg(mean=('virus_end_time', 'mean'),std=('virus_end_time', 'std')).reset_index()
numRuns = timeExt.groupby(['combo_id', 'evofunctionScale'])\
    .agg(n=('virus_end_time', 'size')).reset_index()
meanTimeExt  = meanTimeExt.merge(numRuns,on=['combo_id','evofunctionScale'])




conSim = sqlite3.connect(DBSIM_PATHnull)
curSim = conSim.cursor()
conExt = sqlite3.connect(DBEXT_PATHnull)
curExt = conExt.cursor()
conEsc = sqlite3.connect(DBESC_PATHnull)
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

numEsc2 = []
for runID in runIDs['run_id']:
    # print(runID)
    numEsc2 += [len(pd.read_sql_query("SELECT DISTINCT bstrain_id_escaped, run_id FROM vstrain_escapes \
                           WHERE run_id in ({}) AND t_creation <= 2000"\
                           .format(runID), conEsc)['bstrain_id_escaped'].values)]
    
numEsc2 = pd.DataFrame({'num_escapes':numEsc2,'run_id':runIDs['run_id']})
div = []
for runID in runIDs['run_id']:
    # print(runID)
    div += [len(pd.read_sql_query("SELECT bstrain_id FROM bstrains WHERE run_id in ({}) \
                                AND t_creation <= 2000 ORDER BY bstrain_id"\
                            .format(runID), conSim)['bstrain_id'].values)]
div = pd.DataFrame({'div':div,'run_id':runIDs['run_id']})
numEsc2 = numEsc2.merge(div,on=['run_id'])
numEsc2['prop'] = numEsc2['num_escapes']/numEsc2['div']
numEsc2 = numEsc2.merge(runIDs,on=['run_id']).merge(comboSpace,on=['combo_id'])\
                .sort_values(by=['evofunctionScale','run_id'])
rmRuns=[]
for comboID in np.unique(runIDs['combo_id']):
    rmRuns += sample(list(runIDs[runIDs.combo_id==comboID]['run_id']),100)
#
numEsc2 = numEsc2[~numEsc2['run_id'].isin(rmRuns)]
meanNumEsc2 = numEsc2.groupby(['combo_id','evofunctionScale'])\
                .agg(meanEsc=('num_escapes','mean'),meanDiv=('div','mean'), 
                     meanProp=('prop','mean'),stdProp=('prop','std')).reset_index()
numRuns = numEsc2.groupby(['combo_id', 'evofunctionScale'])\
    .agg(n=('num_escapes', 'size')).reset_index()
meanNumEsc2 = meanNumEsc2.merge(numRuns,on=['combo_id','evofunctionScale'])
# meanNumEsc2['mean'] = np.array(meanNumEsc2['mean'])/np.array(meanNumEsc2['n'])
#
numExt2 = pd.read_sql_query("SELECT viruses, run_id FROM extinction_occurrence WHERE run_id in ({})"
                           .format(', '.join(map(str, runIDs['run_id']))), conExt)
numExt2 = numExt2.merge(runIDs, on=['run_id']).merge(comboSpace, on=['combo_id'])\
    .sort_values(by=['evofunctionScale', 'run_id'])
meanNumExt = numExt2.groupby(['combo_id', 'evofunctionScale'])\
    .agg(mean=('viruses', 'mean')).reset_index()
###
###
timeExt2all = pd.read_sql_query("SELECT virus_end_time, run_id FROM simulation_end_time WHERE run_id in ({})"
                           .format(', '.join(map(str, runIDs['run_id']))), conExt)
timeExt2all = timeExt2all.merge(numExt2[numExt2.viruses==1],on=['run_id'])
timeExt2 = timeExt2all[~timeExt2all.run_id.isin(rmRuns)]
#
meanTimeExt2 = timeExt2.groupby(['combo_id', 'evofunctionScale'])\
    .agg(mean=('virus_end_time', 'mean'),std=('virus_end_time', 'std')).reset_index()
numRuns = timeExt2.groupby(['combo_id', 'evofunctionScale'])\
    .agg(n=('virus_end_time', 'size')).reset_index()
meanTimeExt2  = meanTimeExt2.merge(numRuns,on=['combo_id','evofunctionScale'])


modMeanNumEsc = meanNumEsc.copy()
modMeanTimeExt = meanTimeExt.copy()
modMeanNumEsc = modMeanNumEsc[modMeanNumEsc.evofunctionScale!=.01]
modMeanTimeExt = modMeanTimeExt[modMeanTimeExt.evofunctionScale!=.01]
modMeanNumEsc = pd.concat([modMeanNumEsc,meanNumEsc2[meanNumEsc2.evofunctionScale==0]])\
                    .sort_values(by=['evofunctionScale'])
modMeanTimeExt = pd.concat([modMeanTimeExt,meanTimeExt2[meanNumEsc2.evofunctionScale==0]])\
                    .sort_values(by=['evofunctionScale'])
modMeanNumEsc = modMeanNumEsc[modMeanNumEsc.evofunctionScale<=3.1]
modMeanTimeExt = modMeanTimeExt[modMeanTimeExt.evofunctionScale<=3.1]


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
ax[0].errorbar(modMeanNumEsc['evofunctionScale'], modMeanNumEsc['meanProp'], 
               yerr=modMeanNumEsc['stdProp'], fmt='o',
               color = 'darkslateblue',capsize=2)
# ax[0].scatter(modMeanNumEsc['evofunctionScale'], modMeanNumEsc['meanProp'], color='darkslateblue',
#           s=.2)
# ax[0].set_xticks(modMeanNumEsc['evofunctionScale'])
# ax[0].set_xticklabels(modMeanNumEsc['evofunctionScale'])
ax[0].set_ylim(.4, 1)
ax[0].set_yticks([np.round(i,4) for i in np.arange(.4,1.01,.1)])
# ax[1].scatter(modMeanTimeExt['evofunctionScale'], modMeanTimeExt['mean'], color='darkseagreen',
#           s=30)
# ax[1].plot(modMeanTimeExt['evofunctionScale'], modMeanTimeExt['mean'], color='darkseagreen',
#           linewidth=2,alpha=0.1)
ax[1].errorbar(modMeanTimeExt['evofunctionScale'], modMeanTimeExt['mean'], 
               yerr=modMeanTimeExt['std'], fmt='o',
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

