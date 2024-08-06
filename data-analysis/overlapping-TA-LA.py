#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import sqlite3
import seaborn as sns
import matplotlib.cm as cm
from matplotlib.colors import Normalize
import matplotlib.colors as mc
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.ticker as ticker

# DBSIM_PATH = os.path.join('/Volumes/Yadgah/sylvain-martin-collab/',
#                             '15_MOI3/db-truncate-combinations.sqlite')
# DBTEMP_PATH = os.path.join('/Volumes/Yadgah/sylvain-martin-collab/15_MOI3',
#                             'temporal-adaptation2.sqlite')  # cluster
# DBLOC_PATH = os.path.join('/Volumes/Yadgah/sylvain-martin-collab/15_MOI3',
#                             'local-adaptation.sqlite')  # cluster

DBSIM_PATH = os.path.join('/Volumes/Yadgah/sylvain-martin-collab/',
                            '26_MOI3/sweep_db_gathered.sqlite')
DBTEMP_PATH = os.path.join('/Volumes/Yadgah/sylvain-martin-collab/26_MOI3',
                            'temporal-adaptation2.sqlite')  # cluster
DBHTEMP_PATH = os.path.join('/Volumes/Yadgah/sylvain-martin-collab/26_MOI3',
                            'host-temporal-adaptation.sqlite')  # cluster
DBLOC_PATH = os.path.join('/Volumes/Yadgah/sylvain-martin-collab/26_MOI3',
                            'local-adaptation.sqlite')  # cluster
simDir = '26_MOI3'

conSim = sqlite3.connect(DBSIM_PATH)
curSim = conSim.cursor()
conTemp = sqlite3.connect(DBTEMP_PATH)
curTemp = conTemp.cursor()
conHTemp = sqlite3.connect(DBHTEMP_PATH)
curHTemp = conHTemp.cursor()
conLoc = sqlite3.connect(DBLOC_PATH)
curLoc = conLoc.cursor()

scale = 10
hosts_per_strain = 100
viruses_per_strain = 100
micMutRep = 0
combos = [(1,0)]
# combos = [(2,0),(2,1)]
##

# tempIDs = pd.read_sql_query(
#     "SELECT run_id, combo_id \
#         FROM runs",conSim)\
#             .merge(pd.read_sql_query("SELECT DISTINCT run_id \
#                     FROM vtemporal_adaptation", conTemp),on=['run_id']).drop(columns=['run_id']).drop_duplicates()\
#             .merge(pd.read_sql_query(
#             "SELECT combo_id, init_bcomm_function, evofunctionScale, microbe_mutation_prob_per_spacer \
#                 FROM param_combos",conSim),on=['combo_id'])

tempAdaptMaster = pd.DataFrame({'combo_id': [], 'run_id': [], 'delay': [], 'mean': []})
localAdaptMaster = pd.DataFrame({'combo_id': [], 'run_id': [], 't': [], 'LA': []})
for i in range(0,len(combos)):
    (bcomm,micMutSpacer) = combos[i]
    if bcomm == 1:
        divtype = 'Monomorphic'
    if bcomm == 2:
        divtype = 'Polymorphic'
    comboSpace = pd.read_sql_query(
        "SELECT combo_id \
            FROM param_combos WHERE \
            microbe_mutation_prob_per_spacer = {0} \
            AND init_bcomm_function = {1} \
            AND evofunctionScale = {2} \
            AND n_particles_per_vstrain > 0 \
            ORDER BY combo_id"
        .format(micMutSpacer, bcomm,  scale),
        conSim)
    comboSpace0 = pd.read_sql_query(
        "SELECT combo_id \
            FROM param_combos WHERE \
            microbe_mutation_prob_per_spacer = 0 \
            AND init_bcomm_function = {0} \
            AND evofunctionScale = 0 \
            AND n_particles_per_vstrain > 0 \
            ORDER BY combo_id"
        .format(bcomm),
        conSim)
    cID = comboSpace['combo_id'].values[0]
    cID0 = comboSpace0['combo_id'].values[0] 
    print(cID)
    runIDs = pd.read_sql_query(
    "SELECT combo_id, run_id FROM runs \
        WHERE combo_id in ({})"
    .format(cID), conSim)
    runIDs0 = pd.read_sql_query(
    "SELECT combo_id, run_id FROM runs \
        WHERE combo_id in ({})"
    .format(cID0), conSim)
    ####
    if cID0 not in tempAdaptMaster['combo_id'].values:
        tempAdapt = pd.read_sql_query("SELECT t, delay, vfrequency, bfrequency, run_id \
                                    FROM vtemporal_adaptation \
                                    WHERE run_id in ({})".format(', '.join(map(str, runIDs0['run_id'].values))), conTemp)\
                                    .merge(runIDs0,on=['run_id'])
        tempAdapt = tempAdapt[tempAdapt['t'] > 10]
        tempAdapt['TA'] = tempAdapt['vfrequency'] * tempAdapt['bfrequency']
        tempAdapt = tempAdapt.groupby(['t','delay','run_id', 'combo_id'])\
            .agg(TA=('TA', 'sum')).reset_index()
        tempAdapt = tempAdapt.groupby(['run_id','combo_id','delay'])\
            .agg(mean=('TA', 'mean')).reset_index()\
            .dropna()
        tempAdaptMaster = pd.concat([tempAdaptMaster, tempAdapt], ignore_index=True)
    ##
    ##
    if cID0 not in localAdaptMaster['combo_id'].values:
        localAdapt = pd.DataFrame({'t': [], 'run_id': [], 'intraFitness': [], 'combo_id': []})
        interFitness = pd.DataFrame({'t': [], 'run_id': [], 'interFitness': [], 'combo_id': []})
        for runID in runIDs0['run_id'].values:
            print(runID)
            sumFitness = pd.read_sql_query("SELECT t, vfrequency, bfrequency, run_id \
                        FROM vlocal_adaptation WHERE run_id = {0} AND brun_id = {0}".format(runID), conLoc)\
                        .merge(runIDs0,on=['run_id'])
            sumFitness['intraFitness'] = sumFitness['vfrequency'] * sumFitness['bfrequency']
            sumFitness = sumFitness.groupby(['t','run_id', 'combo_id']).agg(intraFitness=('intraFitness','sum')).reset_index()
            localAdapt = pd.concat([localAdapt, sumFitness], ignore_index=True)
            sumFitness = pd.read_sql_query("SELECT t, vfrequency, bfrequency, brun_id, run_id \
                        FROM vlocal_adaptation WHERE run_id = {0} AND brun_id != {0}".format(runID), conLoc)\
                        .merge(runIDs0,on=['run_id'])
            sumFitness['interFitness'] = sumFitness['vfrequency'] * sumFitness['bfrequency']
            sumFitness = sumFitness.groupby(['t', 'run_id', 'brun_id', 'combo_id'])\
                            .agg(interFitness=('interFitness', 'sum'))\
                                .reset_index()
            sumFitness['interFitness'] = 1/(len(runIDs)-1)*sumFitness['interFitness']
            sumFitness = sumFitness.groupby(['t', 'run_id', 'combo_id'])\
                        .agg(interFitness=('interFitness', 'sum'))\
                            .reset_index()
            interFitness = pd.concat([interFitness, sumFitness], ignore_index=True)
        ##
        ##
        localAdapt = interFitness.merge(localAdapt,on=['combo_id','run_id','t'])
        localAdapt['LA'] = localAdapt['intraFitness'] - localAdapt['interFitness']
        localAdaptMaster = pd.concat([localAdaptMaster, localAdapt], ignore_index=True)
    ##
    ##
    tempAdapt = pd.read_sql_query("SELECT t, delay, vfrequency, bfrequency, run_id \
                                FROM vtemporal_adaptation \
                                WHERE run_id in ({})".format(', '.join(map(str, runIDs['run_id'].values))), conTemp)\
                                .merge(runIDs,on=['run_id'])
    tempAdapt = tempAdapt[tempAdapt['t'] > 10]
    tempAdapt['TA'] = tempAdapt['vfrequency'] * tempAdapt['bfrequency']
    tempAdapt = tempAdapt.groupby(['t','delay','run_id', 'combo_id'])\
        .agg(TA=('TA', 'sum')).reset_index()
    tempAdapt = tempAdapt.groupby(['run_id','combo_id','delay'])\
        .agg(mean=('TA', 'mean')).reset_index()\
        .dropna()
    tempAdaptMaster = pd.concat([tempAdaptMaster, tempAdapt], ignore_index=True)
    ##
    ##
    localAdapt = pd.DataFrame({'t': [], 'run_id': [], 'intraFitness': [], 'combo_id': []})
    interFitness = pd.DataFrame({'t': [], 'run_id': [], 'interFitness': [], 'combo_id': []})
    for runID in runIDs['run_id'].values:
        print(runID)
        sumFitness = pd.read_sql_query("SELECT t, vfrequency, bfrequency, run_id \
                    FROM vlocal_adaptation WHERE run_id = {0} AND brun_id = {0}".format(runID), conLoc)\
                    .merge(runIDs,on=['run_id'])
        sumFitness['intraFitness'] = sumFitness['vfrequency'] * sumFitness['bfrequency']
        sumFitness = sumFitness.groupby(['t','run_id', 'combo_id']).agg(intraFitness=('intraFitness','sum')).reset_index()
        localAdapt = pd.concat([localAdapt, sumFitness], ignore_index=True)
        sumFitness = pd.read_sql_query("SELECT t, vfrequency, bfrequency, brun_id, run_id \
                    FROM vlocal_adaptation WHERE run_id = {0} AND brun_id != {0}".format(runID), conLoc)\
                    .merge(runIDs,on=['run_id'])
        sumFitness['interFitness'] = sumFitness['vfrequency'] * sumFitness['bfrequency']
        sumFitness = sumFitness.groupby(['t', 'run_id', 'brun_id', 'combo_id'])\
                        .agg(interFitness=('interFitness', 'sum'))\
                            .reset_index()
        sumFitness['interFitness'] = 1/(len(runIDs)-1)*sumFitness['interFitness']
        sumFitness = sumFitness.groupby(['t', 'run_id', 'combo_id'])\
                    .agg(interFitness=('interFitness', 'sum'))\
                        .reset_index()
        interFitness = pd.concat([interFitness, sumFitness], ignore_index=True)
    ##
    ##
    localAdapt = interFitness.merge(localAdapt,on=['combo_id','run_id','t'])
    localAdapt['LA'] = localAdapt['intraFitness'] - localAdapt['interFitness']
    localAdaptMaster = pd.concat([localAdaptMaster, localAdapt], ignore_index=True)
    ##
    ##

scale = 10
combo = [1,0]
figsize = (14,7)
resolve = 500
c1 = 'darkorange'
a = 0.3
nThresh = 10 
localxmin = 0
localxmax = 0
tempxmin = 0
tempxmax = 0

fig, axes = plt.subplots(1, 2, figsize=figsize, sharex='col')
plt.subplots_adjust(wspace=0.3)
[bcomm,micMutSpacer] = combo
if bcomm == 1:
    if micMutSpacer == 0:
        c2 = 'darkblue'
    if micMutSpacer == 1:
        c2 = 'deepskyblue'
    c1 = 'darkorange'
    divtype = 'Monomorphic'


##
##
if bcomm == 2:
    if micMutSpacer == 0:
        c2 = 'darkred'
    if micMutSpacer == 1:
        c2 = 'm'
    c1 = 'darkorange'
    divtype = 'Polymorphic'


##
##
comboSpace = pd.read_sql_query(
    "SELECT combo_id \
        FROM param_combos WHERE \
        microbe_mutation_prob_per_spacer = {0} \
        AND init_bcomm_function = {1} \
        AND n_particles_per_vstrain > 0 \
        AND evofunctionScale = {2} \
        ORDER BY combo_id"
    .format(micMutSpacer, bcomm, scale),
    conSim)
comboSpace0 = pd.read_sql_query(
    "SELECT combo_id \
        FROM param_combos WHERE \
        microbe_mutation_prob_per_spacer = 0 \
        AND init_bcomm_function = {0} \
        AND n_particles_per_vstrain > 0 \
        AND evofunctionScale = 0 \
        ORDER BY combo_id"
    .format(bcomm),
    conSim)
cID = comboSpace['combo_id'].values[0]
cID0 = comboSpace0['combo_id'].values[0] 
localAdaptTrunc = localAdaptMaster[(localAdaptMaster.combo_id == cID)].drop_duplicates()
maxTime = 1000
localAdaptTrunc = localAdaptTrunc[localAdaptTrunc['t'] <= maxTime]
tempAdaptTrunc = tempAdaptMaster[(tempAdaptMaster.combo_id == cID)].drop_duplicates()
maxDelay = abs(min(tempAdaptTrunc['delay']))
maxDelay = 1000
tempAdaptTrunc = tempAdaptTrunc[tempAdaptTrunc['delay'] <= maxDelay]
tempAdaptTrunc = tempAdaptTrunc[tempAdaptTrunc['delay'] >= -1*maxDelay]
localAdaptTrunc0 = localAdaptMaster[(localAdaptMaster.combo_id == cID0)].drop_duplicates()
localAdaptTrunc0 = localAdaptTrunc0[localAdaptTrunc0['t'] <= maxTime]
tempAdaptTrunc0 = tempAdaptMaster[(tempAdaptMaster.combo_id == cID0)].drop_duplicates()
maxDelay = abs(min(tempAdaptTrunc0['delay']))
maxDelay = 1000
tempAdaptTrunc0 = tempAdaptTrunc0[tempAdaptTrunc0['delay'] <= maxDelay]
tempAdaptTrunc0 = tempAdaptTrunc0[tempAdaptTrunc0['delay'] >= -1*maxDelay]
for runID in np.unique(tempAdaptTrunc0['run_id']):
    if runID == np.unique(tempAdaptTrunc0['run_id'])[0]:
        axes[0].plot(tempAdaptTrunc0[tempAdaptTrunc0['run_id'] == runID]['delay'], 
                     tempAdaptTrunc0[tempAdaptTrunc0['run_id'] == runID]['mean'], 
                        linewidth=2, color=c1, alpha = a, label = r'$\sigma = 0$')
        axes[1].plot(localAdaptTrunc0['t'], localAdaptTrunc0['LA'],
                    linewidth=2, color=c1, alpha = a, label= r'$\sigma = 0$')
    else:
        axes[0].plot(tempAdaptTrunc0[tempAdaptTrunc0['run_id'] == runID]['delay'], 
                     tempAdaptTrunc0[tempAdaptTrunc0['run_id'] == runID]['mean'], 
                linewidth=2, color=c1, alpha = a)
        axes[1].plot(localAdaptTrunc0['t'], localAdaptTrunc0['LA'],
                    linewidth=2, color=c1, alpha = a)
        


for runID in np.unique(tempAdaptTrunc['run_id']):
    if runID == np.unique(tempAdaptTrunc['run_id'])[0]:
        axes[0].plot(tempAdaptTrunc[tempAdaptTrunc['run_id'] == runID]['delay'],
                     tempAdaptTrunc[tempAdaptTrunc['run_id'] == runID]['mean'],
                        linewidth=2, color=c2, alpha = a, label=''.join([r'$\sigma =$','{}'.format(scale)]))
        axes[1].plot(localAdaptTrunc['t'], localAdaptTrunc['LA'],
                        linewidth=2, color=c2, alpha = a, label = ''.join([r'$\sigma =$','{}'.format(scale)]))
    else:
        axes[0].plot(tempAdaptTrunc[tempAdaptTrunc['run_id'] == runID]['delay'],
                     tempAdaptTrunc[tempAdaptTrunc['run_id'] == runID]['mean'],
                linewidth=2, color=c2, alpha = a)
        axes[1].plot(localAdaptTrunc['t'], localAdaptTrunc['LA'],
                        linewidth=2, color=c2, alpha = a)
        


axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(100))
axes[1].xaxis.set_minor_locator(ticker.MultipleLocator(50))
axes[0].tick_params(axis='x', labelsize=15)
axes[0].tick_params(axis='y', labelsize=15)
axes[1].tick_params(axis='x', labelsize=15)
axes[1].tick_params(axis='y', labelsize=15)
#
axes[0].set_ylabel(ylabel ='Viral Temporal Adaptation (TA)',labelpad=15,fontsize=15)
axes[1].set_ylabel(ylabel ='Viral Local Adaptation (LA)',labelpad=15,fontsize=15)
#
axes[0].legend(loc='upper right', fontsize=15)
axes[0].legend(loc='upper right', fontsize=15)
#
axes[0].set_xlabel(xlabel = r'Delay $\tau$',fontsize=15,labelpad=15)
axes[1].set_xlabel(xlabel = 'Time t',fontsize=15,labelpad=15)
lim = axes[0].get_xlim()
axes[0].set_xlim(-1000, 1000)
lim = axes[1].get_xlim()
axes[1].set_xlim(0, 1000)
axes[0].margins(x=0)
axes[1].margins(x=0)
fig.savefig(os.path.join('/Volumes/Yadgah/sylvain-martin-collab',simDir,'adaptation-{0}-{1}-trunc-overlap.pdf'.format(divtype,int(micMutSpacer))),dpi=resolve)
plt.close('all')
plt.show()





plt.fill_between(X, Y-np.percentile(Y, 95),Y+np.percentile(Y, 95), color="k", alpha = alpha)
plt.fill_between(X, Y-np.percentile(Y, 80),Y+np.percentile(Y, 80), color="r", alpha = alpha)
plt.fill_between(X, Y-np.percentile(Y, 60),Y, color="b", alpha = alpha)


tempAdaptMaster = pd.DataFrame({'combo_id': [], 'run_id': [], 'delay': [], 'mean': []})
localAdaptMaster = pd.DataFrame({'combo_id': [], 'run_id': [], 't': [], 'LA': []}) 