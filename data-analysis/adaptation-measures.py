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

conSim.close()
conSim = None
conTemp.close()
conTemp = None
conHTemp.close()
conHTemp = None
conLoc.close()
conLoc = None

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

tempAdaptMaster = pd.DataFrame({'combo_id': [], 'delay': [], 'mean': [], 'std': [], 'n' : []})
localAdaptMaster = pd.DataFrame({'combo_id': [], 't': [], 'mean': [], 'std': [], 'n' : []})
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
        tempAdapt = tempAdapt.groupby(['delay','run_id','combo_id'])\
            .agg(mean=('TA', 'mean'),std=('TA','std')).reset_index()\
            .dropna()
        tempAdapt = tempAdapt.groupby(['combo_id', 'delay'])\
            .agg(mean=('mean', 'mean'),std=('mean','std'),n=('run_id','size')).reset_index()\
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
        localAdapt = interFitness.merge(localAdapt,on=['t','run_id', 'combo_id'])
        localAdapt['LA'] = localAdapt['intraFitness'] - localAdapt['interFitness']
        localAdapt = localAdapt.groupby(['combo_id', 't'])\
                        .agg(mean=('LA', 'mean'),std=('LA','std'),n=('LA','size')).reset_index()
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
    tempAdapt = tempAdapt.groupby(['delay','run_id', 'combo_id'])\
        .agg(mean=('TA', 'mean'),std=('TA','std')).reset_index()\
        .dropna()
    tempAdapt = tempAdapt.groupby(['combo_id', 'delay'])\
        .agg(mean=('mean', 'mean'),std=('mean','std'),n=('run_id','size')).reset_index()\
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
    localAdapt = interFitness.merge(localAdapt,on=['t','run_id', 'combo_id'])
    localAdapt['LA'] = localAdapt['intraFitness'] - localAdapt['interFitness']
    localAdapt = localAdapt.groupby(['combo_id', 't'])\
                    .agg(mean=('LA', 'mean'),std=('LA','std'),n=('LA','size')).reset_index()
    localAdaptMaster = pd.concat([localAdaptMaster, localAdapt], ignore_index=True)
    ##
    ##


tempAdaptMaster = tempAdaptMaster.drop_duplicates()
localAdaptMaster = localAdaptMaster.drop_duplicates()

#######
## HOST TEMPORAL ADAPTATION
#######

htempAdaptMaster = pd.DataFrame({'combo_id': [], 'delay': [], 'mean': [], 'std': [], 'n' : []})
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
    if cID0 not in htempAdaptMaster['combo_id'].values:
        tempAdapt = pd.read_sql_query("SELECT t, delay, vfrequency, bfrequency, run_id \
                                    FROM btemporal_adaptation \
                                    WHERE run_id in ({})".format(', '.join(map(str, runIDs0['run_id'].values))), conHTemp)\
                                    .merge(runIDs0,on=['run_id']).drop_duplicates()
        tempAdapt = tempAdapt[tempAdapt['t'] > 10]
        tempAdapt['TA'] = tempAdapt['vfrequency'] * tempAdapt['bfrequency']
        tempAdapt = tempAdapt.groupby(['t','delay','run_id', 'combo_id'])\
            .agg(TA=('TA', 'sum')).reset_index()
        tempAdapt = tempAdapt.groupby(['delay','run_id','combo_id'])\
            .agg(mean=('TA', 'mean'),std=('TA','std')).reset_index()\
            .dropna()
        tempAdapt = tempAdapt.groupby(['combo_id', 'delay'])\
            .agg(mean=('mean', 'mean'),std=('mean','std'),n=('run_id','size')).reset_index()\
            .dropna()
        htempAdaptMaster = pd.concat([htempAdaptMaster, tempAdapt], ignore_index=True)
    ##
    ##
    tempAdapt = pd.read_sql_query("SELECT t, delay, vfrequency, bfrequency, run_id \
                                FROM btemporal_adaptation \
                                WHERE run_id in ({})".format(', '.join(map(str, runIDs['run_id'].values))), conHTemp)\
                                .merge(runIDs,on=['run_id']).drop_duplicates()
    tempAdapt = tempAdapt[tempAdapt['t'] > 10]
    tempAdapt['TA'] = tempAdapt['vfrequency'] * tempAdapt['bfrequency']
    tempAdapt = tempAdapt.groupby(['t','delay','run_id', 'combo_id'])\
        .agg(TA=('TA', 'sum')).reset_index()
    tempAdapt = tempAdapt.groupby(['delay','run_id', 'combo_id'])\
        .agg(mean=('TA', 'mean'),std=('TA','std')).reset_index()\
        .dropna()
    tempAdapt = tempAdapt.groupby(['combo_id', 'delay'])\
        .agg(mean=('mean', 'mean'),std=('mean','std'),n=('run_id','size')).reset_index()\
        .dropna()
    htempAdaptMaster = pd.concat([htempAdaptMaster, tempAdapt], ignore_index=True)
    ##




stempAdaptMaster = pd.DataFrame({'combo_id': [], 'delay': [], 'mean': [], 'std': [], 'n' : []})
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
    if cID0 not in stempAdaptMaster['combo_id'].values:
        tempAdapt = pd.read_sql_query("SELECT t, delay, vfrequency, bfrequency, run_id \
                                    FROM btemporal_adaptation \
                                    WHERE run_id in ({})".format(', '.join(map(str, runIDs0['run_id'].values))), conHTemp)\
                                    .merge(runIDs0,on=['run_id']).drop_duplicates()
        tempAdapt = tempAdapt[tempAdapt['t'] > 10]
        tempAdapt['TA'] = (1 - np.array(tempAdapt['vfrequency'])) * np.array(tempAdapt['bfrequency'])
        tempAdapt = tempAdapt.groupby(['t','delay','run_id', 'combo_id'])\
            .agg(TA=('TA', 'sum')).reset_index()
        tempAdapt = tempAdapt.groupby(['delay','run_id','combo_id'])\
            .agg(mean=('TA', 'mean'),std=('TA','std')).reset_index()\
            .dropna()
        tempAdapt = tempAdapt.groupby(['combo_id', 'delay'])\
            .agg(mean=('mean', 'mean'),std=('mean','std'),n=('run_id','size')).reset_index()\
            .dropna()
        stempAdaptMaster = pd.concat([stempAdaptMaster, tempAdapt], ignore_index=True)
    ##
    ##
    tempAdapt = pd.read_sql_query("SELECT t, delay, vfrequency, bfrequency, run_id \
                                FROM btemporal_adaptation \
                                WHERE run_id in ({})".format(', '.join(map(str, runIDs['run_id'].values))), conHTemp)\
                                .merge(runIDs,on=['run_id']).drop_duplicates()
    tempAdapt = tempAdapt[tempAdapt['t'] > 10]
    tempAdapt['TA'] = (1 - np.array(tempAdapt['vfrequency'])) * np.array(tempAdapt['bfrequency'])
    tempAdapt = tempAdapt.groupby(['t','delay','run_id', 'combo_id'])\
        .agg(TA=('TA', 'sum')).reset_index()
    tempAdapt = tempAdapt.groupby(['delay','run_id', 'combo_id'])\
        .agg(mean=('TA', 'mean'),std=('TA','std')).reset_index()\
        .dropna()
    tempAdapt = tempAdapt.groupby(['combo_id', 'delay'])\
        .agg(mean=('mean', 'mean'),std=('mean','std'),n=('run_id','size')).reset_index()\
        .dropna()
    stempAdaptMaster = pd.concat([stempAdaptMaster, tempAdapt], ignore_index=True)
    ##















scale = 10
combo = [(2,1)]
combos = [(1,0),(1,1),(2,0),(2,1)]
combos = [(1,0)]
figsize = (14,7)
resolve = 500
c1 = 'darkorange'
a = 0.3
nThresh = 10 
localxmin = 0
localxmax = 0
tempxmin = 0
tempxmax = 0

for combo in combos:
    # (bcomm,micMutSpacer) = combo[i]
    fig, axes = plt.subplots(1, 2, figsize=figsize, sharex='col')
    plt.subplots_adjust(wspace=0.3)
    (bcomm,micMutSpacer) = combo
    if bcomm == 1:
        if micMutSpacer == 0:
            c2 = 'darkblue'
        if micMutSpacer == 1:
            c2 = 'deepskyblue'
        c1 = 'darkorange'
        divtype = 'Monomorphic'
    if bcomm == 2:
        if micMutSpacer == 0:
            c2 = 'darkred'
        if micMutSpacer == 1:
            c2 = 'm'
        c1 = 'darkorange'
        divtype = 'Polymorphic'
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
    localAdaptTrunc = localAdaptMaster[(localAdaptMaster.n >= nThresh) & (localAdaptMaster.combo_id == cID)].drop_duplicates()
    maxTime = 600
    localAdaptTrunc = localAdaptTrunc[localAdaptTrunc['t'] <= maxTime]
    tempAdaptTrunc = tempAdaptMaster[(tempAdaptMaster.n >= nThresh) & (tempAdaptMaster.combo_id == cID)].drop_duplicates()
    maxDelay = abs(min(tempAdaptTrunc['delay']))
    maxDelay = 300
    tempAdaptTrunc = tempAdaptTrunc[tempAdaptTrunc['delay'] <= maxDelay]
    tempAdaptTrunc = tempAdaptTrunc[tempAdaptTrunc['delay'] >= -1*maxDelay]
    localAdaptTrunc0 = localAdaptMaster[(localAdaptMaster.n >= nThresh) & (localAdaptMaster.combo_id == cID0)].drop_duplicates()
    localAdaptTrunc0 = localAdaptTrunc0[localAdaptTrunc0['t'] <= maxTime]
    tempAdaptTrunc0 = tempAdaptMaster[(tempAdaptMaster.n >= nThresh) & (tempAdaptMaster.combo_id == cID0)].drop_duplicates()
    tempAdaptTrunc0 = tempAdaptTrunc0[tempAdaptTrunc0['delay'] <= maxDelay]
    tempAdaptTrunc0 = tempAdaptTrunc0[tempAdaptTrunc0['delay'] >= -1*maxDelay]
    axes[0].fill_between(tempAdaptTrunc0['delay'], tempAdaptTrunc0['mean']-tempAdaptTrunc0['std'],
                            tempAdaptTrunc0['mean'] + tempAdaptTrunc0['std'], color=c1, alpha=a)
    axes[0].plot(tempAdaptTrunc0['delay'], tempAdaptTrunc0['mean'], 
                    linewidth=2, color=c1, label = r'$\sigma = 0$')
    axes[0].fill_between(tempAdaptTrunc['delay'], tempAdaptTrunc['mean']-tempAdaptTrunc['std'], 
                        tempAdaptTrunc['mean'] + tempAdaptTrunc['std'], color=c2, alpha=a)
    axes[0].plot(tempAdaptTrunc['delay'],tempAdaptTrunc['mean'],
                    linewidth=2, color=c2, label=''.join([r'$\sigma =$','{}'.format(scale)]))
    axes[1].fill_between(localAdaptTrunc0['t'], localAdaptTrunc0['mean']-localAdaptTrunc0['std'],
                            localAdaptTrunc0['mean'] + localAdaptTrunc0['std'], color=c1, alpha=a)
    axes[1].plot(localAdaptTrunc0['t'], localAdaptTrunc0['mean'],
                    linewidth=2, color=c1, label= r'$\sigma = 0$')
    axes[1].fill_between(localAdaptTrunc['t'], localAdaptTrunc['mean']-localAdaptTrunc['std'],
                        localAdaptTrunc['mean'] + localAdaptTrunc['std'], color=c2, alpha=a)
    axes[1].plot(localAdaptTrunc['t'], localAdaptTrunc['mean'],
                    linewidth=2, color=c2, label = ''.join([r'$\sigma =$','{}'.format(scale)]))
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
    axes[0].set_xlim(-1*maxDelay, maxDelay)
    lim = axes[1].get_xlim()
    axes[1].set_xlim(0, maxTime)
    axes[0].margins(x=0)
    axes[1].margins(x=0)
    fig.savefig(os.path.join('/Volumes/Yadgah/sylvain-martin-collab',simDir,'adaptation-{0}-{1}_{2}-trunc.pdf'.format(divtype,int(micMutSpacer),nThresh)),dpi=resolve)
    plt.close('all')

plt.show()





####
# Local Adaptation Decomposition
####
nThresh = 20
micMutSpacer = 0 
bcomm =  1
scale = 10

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
#
localAdapt0 = pd.DataFrame({'t': [], 'run_id': [], 'intraFitness': [], 'combo_id': []})
interFitness = pd.DataFrame({'t': [], 'run_id': [], 'interFitness': [], 'combo_id': []})
for runID in runIDs0['run_id'].values:
    print(runID)
    sumFitness = pd.read_sql_query("SELECT t, vfrequency, bfrequency, run_id \
                FROM vlocal_adaptation WHERE run_id = {0} AND brun_id = {0}".format(runID), conLoc)\
                .merge(runIDs0,on=['run_id'])
    sumFitness['intraFitness'] = sumFitness['vfrequency'] * sumFitness['bfrequency']
    sumFitness = sumFitness.groupby(['t','run_id', 'combo_id']).agg(intraFitness=('intraFitness','sum')).reset_index()
    localAdapt0 = pd.concat([localAdapt0, sumFitness], ignore_index=True)
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
localAdapt0 = interFitness.merge(localAdapt0,on=['t','run_id', 'combo_id'])
localAdapt0['LA'] = localAdapt0['intraFitness'] - localAdapt0['interFitness']
localAdapt0 = localAdapt0.groupby(['t'])\
                .agg(mean=('LA', 'mean'),std=('LA','std'),
                        meanIntra=('intraFitness', 'mean'),stdIntra=('intraFitness','std'), 
                        meanInter=('interFitness', 'mean'),stdInter=('interFitness','std'),
                        n=('LA','size')).reset_index()
##
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
localAdapt = interFitness.merge(localAdapt,on=['t','run_id', 'combo_id'])
localAdapt['LA'] = localAdapt['intraFitness'] - localAdapt['interFitness']
localAdapt = localAdapt.groupby(['t'])\
                .agg(mean=('LA', 'mean'),std=('LA','std'),
                        meanIntra=('intraFitness', 'mean'),stdIntra=('intraFitness','std'), 
                        meanInter=('interFitness', 'mean'),stdInter=('interFitness','std'),
                        n=('LA','size')).reset_index()
#
localAdaptTrunc0 = localAdapt0[localAdapt0.n >= nThresh]
localAdaptTrunc = localAdapt[localAdapt.n >= nThresh]
######
######
virus_total = pd.read_sql_query("SELECT run_id, t, viral_abundance FROM summary \
                            WHERE run_id in ({}) AND viral_abundance > 0"
                            .format(', '.join(map(str, runIDs0['run_id']))), conSim)\
                            .rename(columns={"viral_abundance": "vtotal"})
vstats0 = virus_total.groupby(['t'])\
.agg(exp_vtotal=('vtotal', 'mean'), std_vtotal=('vtotal', 'std'), n=('vtotal', 'size'))\
.reset_index()
vstatsTrunc0 = vstats0[vstats0.n >= nThresh]
vstatsTrunc0 = vstatsTrunc0[vstatsTrunc0.t <= max(localAdaptTrunc0['t'])]
###
virus_total = pd.read_sql_query("SELECT run_id, t, viral_abundance FROM summary \
                        WHERE run_id in ({}) AND viral_abundance > 0"
                        .format(', '.join(map(str, runIDs['run_id']))), conSim)\
                        .rename(columns={"viral_abundance": "vtotal"})
vstats = virus_total.groupby(['t'])\
.agg(exp_vtotal=('vtotal', 'mean'), std_vtotal=('vtotal', 'std'), n=('vtotal', 'size'))\
.reset_index()
vstatsTrunc = vstats[vstats.n >= nThresh]
vstatsTrunc = vstatsTrunc[vstatsTrunc.t <= max(localAdaptTrunc['t'])]
######
######
######
fig, ax = plt.subplots(2, 2, sharex='col',
                    figsize=(14, 14))
axes = [ax[0,0], ax[0,1], ax[1,0], ax[1,1], ax[0,0].twinx(), ax[0,1].twinx()]
axes[1].sharey(axes[0])
axes[3].sharey(axes[2])
axes[5].sharey(axes[4])
plt.subplots_adjust(hspace=0.3, wspace=0.4)
a = 0.1
axes[0].fill_between(vstatsTrunc0['t'],
                vstatsTrunc0['exp_vtotal'] -
                vstatsTrunc0['std_vtotal'],
                vstatsTrunc0['exp_vtotal'] +
                vstatsTrunc0['std_vtotal'], color='grey', alpha=0.1)
axes[0].plot(vstatsTrunc0['t'],
        vstatsTrunc0['exp_vtotal'],
        linewidth=2, color='grey', label='Viral Abund.', linestyle='solid', alpha=0.75)
axes[1].fill_between(vstatsTrunc['t'],
                vstatsTrunc['exp_vtotal'] -
                vstatsTrunc['std_vtotal'],
                vstatsTrunc['exp_vtotal'] +
                vstatsTrunc['std_vtotal'], color='grey', alpha=0.1)
axes[1].plot(vstatsTrunc['t'],
        vstatsTrunc['exp_vtotal'],
        linewidth=2, color='grey', linestyle='solid', alpha=0.75)
for i in [0,1]:
    axes[i].yaxis.set_label_position("right")
    axes[i].yaxis.tick_right()
    lim = axes[i].get_ylim()
    axes[i].set_ylim(0, lim[1])
axes[0].set_ylabel('')
# axes[0].tick_params(axis='y', labelcolor='w', top=False,
#                 bottom=False, left=False, right=False)
axes[1].set_ylabel(ylabel='Mean Viral Abundance', rotation=270, labelpad=35,fontsize=20)
####
###
axes[4].fill_between(localAdaptTrunc0['t'], localAdaptTrunc0['meanIntra']-localAdaptTrunc0['stdIntra'],
                    localAdaptTrunc0['meanIntra'] + localAdaptTrunc0['stdIntra'], color='darkorange', alpha=a)
axes[4].plot(localAdaptTrunc0['t'], localAdaptTrunc0['meanIntra'],
                linewidth=2, color='darkorange', label = r'Sym. $\sigma = 0$')
axes[5].fill_between(localAdaptTrunc['t'], localAdaptTrunc['meanIntra']-localAdaptTrunc['stdIntra'],
                    localAdaptTrunc['meanIntra'] + localAdaptTrunc['stdIntra'], color='darkblue', alpha=a)
axes[5].plot(localAdaptTrunc['t'], localAdaptTrunc['meanIntra'],
                linewidth=2, color='darkblue', label = r'Sym. $\sigma = 10$')
axes[4].yaxis.tick_left()
axes[5].yaxis.tick_left()
axes[4].set_ylabel(ylabel='Mean Viral Fitness',labelpad=15,fontsize=20)
axes[4].yaxis.set_label_position("left")
axes[5].yaxis.set_label_position("left")

###
###
axes[2].fill_between(localAdaptTrunc0['t'], localAdaptTrunc0['meanIntra']-localAdaptTrunc0['stdIntra'],
                    localAdaptTrunc0['meanIntra'] + localAdaptTrunc0['stdIntra'], color='darkorange', alpha=a)
axes[2].plot(localAdaptTrunc0['t'], localAdaptTrunc0['meanIntra'],
                linewidth=2, color='darkorange', label = r'Sym. $\sigma = 0$')
###
axes[3].fill_between(localAdaptTrunc['t'], localAdaptTrunc['meanIntra']-localAdaptTrunc['stdIntra'],
                    localAdaptTrunc['meanIntra'] + localAdaptTrunc['stdIntra'], color='darkblue', alpha=a)
axes[3].plot(localAdaptTrunc['t'], localAdaptTrunc['meanIntra'],
                linewidth=2, color='darkblue', label = r'Sym. $\sigma = 10$')
# axes[1].set_yticklabels([])
axes[2].fill_between(localAdaptTrunc0['t'], localAdaptTrunc0['meanInter']-localAdaptTrunc0['stdInter'],
                localAdaptTrunc0['meanInter'] + localAdaptTrunc0['stdInter'], color='darkred', alpha=a)
axes[2].plot(localAdaptTrunc0['t'], localAdaptTrunc0['meanInter'],
                linewidth=2, color='darkred', label = r'Allo. $\sigma = 0$')
###
axes[3].fill_between(localAdaptTrunc['t'], localAdaptTrunc['meanInter']-localAdaptTrunc['stdInter'],
                    localAdaptTrunc['meanInter'] + localAdaptTrunc['stdInter'], color='darkgreen', alpha=a)
axes[3].plot(localAdaptTrunc['t'], localAdaptTrunc['meanInter'],
                linewidth=2, color='darkgreen', label = r'Allo. $\sigma = 10$')
###
axes[2].set_ylabel(ylabel='Mean Viral Fitness',labelpad=15,fontsize=20)
axes[2].yaxis.tick_left()
axes[2].yaxis.set_label_position("left")
###
handles = []
labels = []
handle, label = axes[0].get_legend_handles_labels()
handles.extend(handle)
labels.extend(label)
handle, label = axes[4].get_legend_handles_labels()
handles.extend(handle)
labels.extend(label)
axes[0].legend(handles, labels, loc='upper right', fontsize=20)
###
handles = []
labels = []
handle, label = axes[5].get_legend_handles_labels()
handles.extend(handle)
labels.extend(label)
axes[1].legend(handles, labels, loc='upper right', fontsize=20)
###
axes[2].legend(loc='upper right', fontsize=20)
axes[3].legend(loc='upper right', fontsize=20)
#
axes[0].yaxis.get_offset_text().set_fontsize(20)
axes[1].yaxis.get_offset_text().set_fontsize(20)
axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(50)) 
axes[1].xaxis.set_minor_locator(ticker.MultipleLocator(50)) 
axes[2].set_xlabel(xlabel = 'Time t',fontsize=20,labelpad=15)
axes[3].set_xlabel(xlabel = 'Time t',fontsize=20,labelpad=15)
axes[0].tick_params(axis='x', labelsize=20)
axes[0].tick_params(axis='y', labelsize=20)
axes[1].tick_params(axis='x', labelsize=20)
axes[1].tick_params(axis='y', labelsize=20)
axes[4].tick_params(axis='x', labelsize=20)
axes[4].tick_params(axis='y', labelsize=20)
axes[5].tick_params(axis='x', labelsize=20)
axes[5].tick_params(axis='y', labelsize=20)
axes[2].tick_params(axis='x', labelsize=20)
axes[2].tick_params(axis='y', labelsize=20)
axes[3].tick_params(axis='x', labelsize=20)
axes[3].tick_params(axis='y', labelsize=20)
#
lim = axes[0].get_xlim()
axes[0].set_xlim(0, 1000)
lim = axes[1].get_xlim()
axes[1].set_xlim(0, 1000)
axes[0].margins(x=0)
axes[1].margins(x=0)
fig.savefig(os.path.join('/Volumes/Yadgah/sylvain-martin-collab',simDir,'supp-LA-decomposition_{0}.pdf'.format(nThresh)),dpi=resolve)
plt.close('all')
plt.show()


#####
#####



scale = 10
resolve = 500
c0 = 'gold'
cScale = 'darkblue'
cMut = 'darkred'
cHost0 = 'darkgreen'
cHostScale = 'darkblue'
cHostMut = 'darkred'
a = 0.3
nThresh = 10 
tempxmin = 0
tempxmax = 0
figsize = (17,22)
fig = plt.figure(figsize=figsize)
outer_grid = fig.add_gridspec(2, 1, wspace=1.2, hspace=.3)
for bcomm in [1,2]:
    if bcomm == 1:
        divtype = 'Monomorphic'
    if bcomm == 2:
        divtype = 'Polymorphic'
    comboSpace = pd.read_sql_query(
        "SELECT combo_id \
            FROM param_combos WHERE \
            microbe_mutation_prob_per_spacer = {0} \
            AND init_bcomm_function = {1} \
            AND n_particles_per_vstrain > 0 \
            AND evofunctionScale = {2} \
            ORDER BY combo_id"
        .format(0, bcomm, scale),
        conSim)
    cID = comboSpace['combo_id'].values[0]
    comboSpaceMut = pd.read_sql_query(
        "SELECT combo_id \
            FROM param_combos WHERE \
            microbe_mutation_prob_per_spacer = {0} \
            AND init_bcomm_function = {1} \
            AND n_particles_per_vstrain > 0 \
            AND evofunctionScale = {2} \
            ORDER BY combo_id"
        .format(1, bcomm, scale),
        conSim)
    cIDMut = comboSpaceMut['combo_id'].values[0]
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
    cID0 = comboSpace0['combo_id'].values[0]
    ##
    ## 
    tempAdaptTrunc = tempAdaptMaster[(tempAdaptMaster.n >= nThresh) & (tempAdaptMaster.combo_id == cID)].drop_duplicates()
    maxDelay = abs(min(tempAdaptTrunc['delay']))
    tempAdaptTrunc = tempAdaptTrunc[tempAdaptTrunc['delay'] <= maxDelay]
    htempAdaptTrunc = htempAdaptMaster[(htempAdaptMaster.n >= nThresh) & (htempAdaptMaster.combo_id == cID)].drop_duplicates()
    # maxDelay = 500
    # htempAdaptTrunc = htempAdaptTrunc[htempAdaptTrunc['delay'] <= maxDelay]
    # htempAdaptTrunc = htempAdaptTrunc[htempAdaptTrunc['delay'] >= -1*maxDelay]
    tempAdaptMutTrunc = tempAdaptMaster[(tempAdaptMaster.n >= nThresh) & (tempAdaptMaster.combo_id == cIDMut)].drop_duplicates()
    maxDelay = abs(min(tempAdaptMutTrunc['delay']))
    tempAdaptMutTrunc = tempAdaptMutTrunc[tempAdaptMutTrunc['delay'] <= maxDelay]
    htempAdaptMutTrunc = htempAdaptMaster[(htempAdaptMaster.n >= nThresh) & (htempAdaptMaster.combo_id == cIDMut)].drop_duplicates()
    # maxDelay = 500
    # htempAdaptMutTrunc = htempAdaptMutTrunc[htempAdaptMutTrunc['delay'] <= maxDelay]
    # htempAdaptMutTrunc = htempAdaptMutTrunc[htempAdaptMutTrunc['delay'] >= -1*maxDelay]
    tempAdaptTrunc0 = tempAdaptMaster[(tempAdaptMaster.n >= nThresh) & (tempAdaptMaster.combo_id == cID0)].drop_duplicates()
    maxDelay = abs(min(tempAdaptTrunc0['delay']))
    tempAdaptTrunc0 = tempAdaptTrunc0[tempAdaptTrunc0['delay'] <= maxDelay]
    htempAdaptTrunc0 = htempAdaptMaster[(htempAdaptMaster.n >= nThresh) & (htempAdaptMaster.combo_id == cID0)].drop_duplicates()
    htempAdaptTrunc0 = htempAdaptTrunc0[htempAdaptTrunc0['delay'] <= maxDelay]
    htempAdaptTrunc0 = htempAdaptTrunc0[htempAdaptTrunc0['delay'] >= -1*maxDelay]
    ##
    ##
    inner_grid = outer_grid[bcomm-1].subgridspec(2, 6, hspace=0.4, wspace=.6)
    ax1 = fig.add_subplot(inner_grid[0,1:3]) # row 0 with axes spanning 2 cols on odds
    ax2 = fig.add_subplot(inner_grid[0,3:5])
    ax3 = fig.add_subplot(inner_grid[1,0:2]) # row 0 with axes spanning 2 cols on evens
    ax4 = fig.add_subplot(inner_grid[1,2:4])
    ax5 = fig.add_subplot(inner_grid[1,4:])
    # axes = [ax1, ax2, ax3, ax4, ax5,
    #     ax3.twinx(), ax4.twinx(), ax5.twinx()]
    axes = [ax1, ax2, ax3, ax4, ax5]

    axes[0].fill_between(tempAdaptTrunc0['delay'], tempAdaptTrunc0['mean']-tempAdaptTrunc0['std'],
                            tempAdaptTrunc0['mean'] + tempAdaptTrunc0['std'], color=c0, alpha=a)
    axes[0].plot(tempAdaptTrunc0['delay'], tempAdaptTrunc0['mean'], 
                    linewidth=1.5, color=c0, label = r'virus, $\sigma = 0$')
    axes[1].fill_between(tempAdaptTrunc0['delay'], tempAdaptTrunc0['mean']-tempAdaptTrunc0['std'],
                            tempAdaptTrunc0['mean'] + tempAdaptTrunc0['std'], color=c0, alpha=a)
    axes[1].plot(tempAdaptTrunc0['delay'], tempAdaptTrunc0['mean'], 
                    linewidth=1.5, color=c0)
    axes[2].fill_between(tempAdaptTrunc0['delay'], tempAdaptTrunc0['mean']-tempAdaptTrunc0['std'],
                        tempAdaptTrunc0['mean'] + tempAdaptTrunc0['std'], color=c0, alpha=a)
    axes[2].plot(tempAdaptTrunc0['delay'], tempAdaptTrunc0['mean'], 
                    linewidth=1.5, color=c0, label = r'virus, $\sigma = 0$')
    axes[2].fill_between(htempAdaptTrunc0['delay'], htempAdaptTrunc0['mean']-htempAdaptTrunc0['std'],
                        htempAdaptTrunc0['mean'] + htempAdaptTrunc0['std'], color=cHost0, alpha=a)
    axes[2].plot(htempAdaptTrunc0['delay'], htempAdaptTrunc0['mean'], 
                    linewidth=1.5, color=cHost0, label = r'host, $\sigma = 0$')
    #
    axes[0].fill_between(tempAdaptTrunc['delay'], tempAdaptTrunc['mean']-tempAdaptTrunc['std'], 
                        tempAdaptTrunc['mean'] + tempAdaptTrunc['std'], color=cScale, alpha=a)
    axes[0].plot(tempAdaptTrunc['delay'],tempAdaptTrunc['mean'],
                    linewidth=1.5, color=cScale, label=''.join(['virus, ', r'$\sigma =$','{}, '.format(scale), r'$\mu_s = 0$']))
    axes[3].fill_between(tempAdaptTrunc['delay'], tempAdaptTrunc['mean']-tempAdaptTrunc['std'], 
                        tempAdaptTrunc['mean'] + tempAdaptTrunc['std'], color=cScale, alpha=a)
    axes[3].plot(tempAdaptTrunc['delay'],tempAdaptTrunc['mean'],
                    linewidth=1.5, color=cScale, label=''.join(['virus, ', r'$\sigma =$','{}, '.format(scale), r'$\mu_s = 0$']))
    axes[4].fill_between(htempAdaptTrunc['delay'], htempAdaptTrunc['mean']-htempAdaptTrunc['std'], 
                        htempAdaptTrunc['mean'] + htempAdaptTrunc['std'], color=cHostScale, alpha=a)
    axes[4].plot(htempAdaptTrunc['delay'],htempAdaptTrunc['mean'],
                    linewidth=1.5, color=cHostScale, label=''.join(['host, ', r'$\sigma =$','{}, '.format(scale), r'$\mu_s = 0$']))
    #
    axes[1].fill_between(tempAdaptMutTrunc['delay'], tempAdaptMutTrunc['mean']-tempAdaptMutTrunc['std'], 
                        tempAdaptMutTrunc['mean'] + tempAdaptMutTrunc['std'], color=cMut, alpha=a)
    axes[1].plot(tempAdaptMutTrunc['delay'],tempAdaptMutTrunc['mean'],
                    linewidth=1.5, color=cMut, label=''.join(['virus, ', r'$\sigma =$','{}, '.format(scale), r'$\mu_s = 1$']))
    axes[3].fill_between(tempAdaptMutTrunc['delay'], tempAdaptMutTrunc['mean']-tempAdaptMutTrunc['std'], 
                        tempAdaptMutTrunc['mean'] + tempAdaptMutTrunc['std'], color=cMut, alpha=a)
    axes[3].plot(tempAdaptMutTrunc['delay'],tempAdaptMutTrunc['mean'],
                    linewidth=1.5, color=cMut, label=''.join(['virus, ', r'$\sigma =$','{}, '.format(scale), r'$\mu_s = 1$']))
    axes[4].fill_between(htempAdaptMutTrunc['delay'], htempAdaptMutTrunc['mean']-htempAdaptMutTrunc['std'], 
                        htempAdaptMutTrunc['mean'] + htempAdaptMutTrunc['std'], color=cHostMut, alpha=a)
    axes[4].plot(htempAdaptMutTrunc['delay'],htempAdaptMutTrunc['mean'],
                    linewidth=1.5, color=cHostMut, label=''.join(['host, ', r'$\sigma =$','{}, '.format(scale), r'$\mu_s = 1$']))
    #
    axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(100))
    axes[1].xaxis.set_minor_locator(ticker.MultipleLocator(100))
    axes[2].xaxis.set_minor_locator(ticker.MultipleLocator(100))
    axes[3].xaxis.set_minor_locator(ticker.MultipleLocator(100))
    axes[4].xaxis.set_minor_locator(ticker.MultipleLocator(100))

    axes[0].tick_params(axis='x', labelsize=15)
    axes[0].tick_params(axis='y', labelsize=15)
    axes[1].tick_params(axis='x', labelsize=15)
    axes[1].tick_params(axis='y', labelsize=15)
    axes[2].tick_params(axis='x', labelsize=15)
    axes[2].tick_params(axis='y', labelsize=15)
    axes[3].tick_params(axis='x', labelsize=15)
    axes[3].tick_params(axis='y', labelsize=15)
    axes[4].tick_params(axis='x', labelsize=15)
    axes[4].tick_params(axis='y', labelsize=15)

    axes[0].legend(loc='upper right', fontsize=15)
    axes[1].legend(loc='upper right', fontsize=15)
    axes[2].legend(loc='lower left', fontsize=15)
    axes[3].legend(loc='upper right', fontsize=15)
    axes[4].legend(loc='lower left', fontsize=15)
    axes[0].set_xlabel(xlabel = r'Delay $\tau$',fontsize=15,labelpad=15)
    axes[0].set_ylabel(ylabel ='Temporal\nAdaptation (TA)',labelpad=15,fontsize=15)
    axes[2].set_ylabel(ylabel ='Temporal\nAdaptation (TA)',labelpad=15,fontsize=15)
    

fig.savefig(os.path.join('/Volumes/Yadgah/sylvain-martin-collab',simDir,'supp-TA-all_{}.pdf'.format(nThresh)),dpi=resolve)
plt.close('all')



scale = 10
resolve = 500
c0 = 'darkorange'
cScale = 'darkblue'
cMut = 'darkred'
cHost0 = 'darkorange'
cHostScale = 'darkgreen'
cHostMut = 'indigo'
a = 0.3
nThresh = 20 
tempxmin = 0
tempxmax = 0
figsize = (17,17)
fig = plt.figure(figsize=figsize)
outer_grid = fig.add_gridspec(2, 1, hspace=.3)
for bcomm in [1,2]:
    if bcomm == 1:
        divtype = 'Monomorphic'
    if bcomm == 2:
        divtype = 'Polymorphic'
    comboSpace = pd.read_sql_query(
        "SELECT combo_id \
            FROM param_combos WHERE \
            microbe_mutation_prob_per_spacer = {0} \
            AND init_bcomm_function = {1} \
            AND n_particles_per_vstrain > 0 \
            AND evofunctionScale = {2} \
            ORDER BY combo_id"
        .format(0, bcomm, scale),
        conSim)
    cID = comboSpace['combo_id'].values[0]
    comboSpaceMut = pd.read_sql_query(
        "SELECT combo_id \
            FROM param_combos WHERE \
            microbe_mutation_prob_per_spacer = {0} \
            AND init_bcomm_function = {1} \
            AND n_particles_per_vstrain > 0 \
            AND evofunctionScale = {2} \
            ORDER BY combo_id"
        .format(1, bcomm, scale),
        conSim)
    cIDMut = comboSpaceMut['combo_id'].values[0]
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
    cID0 = comboSpace0['combo_id'].values[0]
    ##
    ## 
    localAdaptTrunc = localAdaptMaster[(localAdaptMaster.n >= nThresh) & (localAdaptMaster.combo_id == cID)].drop_duplicates()
    localAdaptMutTrunc = localAdaptMaster[(localAdaptMaster.n >= nThresh) & (localAdaptMaster.combo_id == cIDMut)].drop_duplicates()
    localAdaptTrunc0 = localAdaptMaster[(localAdaptMaster.n >= nThresh) & (localAdaptMaster.combo_id == cID0)].drop_duplicates()
    ##
    ##
    inner_grid = outer_grid[bcomm-1].subgridspec(1, 2, wspace=.3)
    ax1 = fig.add_subplot(inner_grid[0,0:1]) # row 0 with axes spanning 2 cols on odds
    ax2 = fig.add_subplot(inner_grid[0,1:2])
    axes = [ax1, ax2]
    ##
    axes[0].fill_between(localAdaptTrunc0['t'], localAdaptTrunc0['mean']-localAdaptTrunc0['std'],
                            localAdaptTrunc0['mean'] + localAdaptTrunc0['std'], color=c0, alpha=a)
    axes[0].plot(localAdaptTrunc0['t'], localAdaptTrunc0['mean'], 
                    linewidth=1.5, color=c0, label = r'virus, $\sigma = 0$')
    axes[1].fill_between(localAdaptTrunc0['t'], localAdaptTrunc0['mean']-localAdaptTrunc0['std'],
                            localAdaptTrunc0['mean'] + localAdaptTrunc0['std'], color=c0, alpha=a)
    axes[1].plot(localAdaptTrunc0['t'], localAdaptTrunc0['mean'], 
                    linewidth=1.5, color=c0)
    #
    axes[0].fill_between(localAdaptTrunc['t'], localAdaptTrunc['mean']-localAdaptTrunc['std'], 
                        localAdaptTrunc['mean'] + localAdaptTrunc['std'], color=cScale, alpha=a)
    axes[0].plot(localAdaptTrunc['t'],localAdaptTrunc['mean'],
                    linewidth=1.5, color=cScale, label=''.join(['virus, ', r'$\sigma =$','{}, '.format(scale), r'$\mu_s = 0$']))
    #
    axes[1].fill_between(localAdaptMutTrunc['t'], localAdaptMutTrunc['mean']-localAdaptMutTrunc['std'], 
                        localAdaptMutTrunc['mean'] + localAdaptMutTrunc['std'], color=cMut, alpha=a)
    axes[1].plot(localAdaptMutTrunc['t'],localAdaptMutTrunc['mean'],
                    linewidth=1.5, color=cMut, label=''.join(['virus, ', r'$\sigma =$','{}, '.format(scale), r'$\mu_s = 1$']))
    #
    axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(100))
    axes[1].xaxis.set_minor_locator(ticker.MultipleLocator(100))
    ##
    axes[0].tick_params(axis='x', labelsize=15)
    axes[0].tick_params(axis='y', labelsize=15)
    axes[1].tick_params(axis='x', labelsize=15)
    axes[1].tick_params(axis='y', labelsize=15)
    ##
    axes[0].legend(loc='upper right', fontsize=15)
    axes[1].legend(loc='upper right', fontsize=15)
    axes[0].set_xlabel(xlabel = r'Time $t$',fontsize=15,labelpad=15)
    axes[0].set_ylabel(ylabel ='Local\nAdaptation (LA)',labelpad=15,fontsize=15)


fig.savefig(os.path.join('/Volumes/Yadgah/sylvain-martin-collab',simDir,'supp-LA-all_{0}.pdf'.format(nThresh)),dpi=resolve)
plt.close('all')





scale = 10
resolve = 500
c0 = 'gold'
cScale = 'darkblue'
cMut = 'darkred'
cHost0 = 'darkgreen'
cHostScale = 'darkblue'
cHostMut = 'darkred'
a = 0.3
nThresh = 10 
tempxmin = 0
tempxmax = 0
figsize = (17,22)
fig = plt.figure(figsize=figsize)
outer_grid = fig.add_gridspec(2, 1, wspace=1.2, hspace=.3)
for bcomm in [1,2]:
    if bcomm == 1:
        divtype = 'Monomorphic'
    if bcomm == 2:
        divtype = 'Polymorphic'
    comboSpace = pd.read_sql_query(
        "SELECT combo_id \
            FROM param_combos WHERE \
            microbe_mutation_prob_per_spacer = {0} \
            AND init_bcomm_function = {1} \
            AND n_particles_per_vstrain > 0 \
            AND evofunctionScale = {2} \
            ORDER BY combo_id"
        .format(0, bcomm, scale),
        conSim)
    cID = comboSpace['combo_id'].values[0]
    comboSpaceMut = pd.read_sql_query(
        "SELECT combo_id \
            FROM param_combos WHERE \
            microbe_mutation_prob_per_spacer = {0} \
            AND init_bcomm_function = {1} \
            AND n_particles_per_vstrain > 0 \
            AND evofunctionScale = {2} \
            ORDER BY combo_id"
        .format(1, bcomm, scale),
        conSim)
    cIDMut = comboSpaceMut['combo_id'].values[0]
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
    cID0 = comboSpace0['combo_id'].values[0]
    ##
    ## 
    tempAdaptTrunc = tempAdaptMaster[(tempAdaptMaster.n >= nThresh) & (tempAdaptMaster.combo_id == cID)].drop_duplicates()
    maxDelay = abs(min(tempAdaptTrunc['delay']))
    tempAdaptTrunc = tempAdaptTrunc[tempAdaptTrunc['delay'] <= maxDelay]
    stempAdaptTrunc = stempAdaptMaster[(stempAdaptMaster.n >= nThresh) & (stempAdaptMaster.combo_id == cID)].drop_duplicates()
    # maxDelay = 500
    # stempAdaptTrunc = stempAdaptTrunc[stempAdaptTrunc['delay'] <= maxDelay]
    # stempAdaptTrunc = stempAdaptTrunc[stempAdaptTrunc['delay'] >= -1*maxDelay]
    tempAdaptMutTrunc = tempAdaptMaster[(tempAdaptMaster.n >= nThresh) & (tempAdaptMaster.combo_id == cIDMut)].drop_duplicates()
    maxDelay = abs(min(tempAdaptMutTrunc['delay']))
    tempAdaptMutTrunc = tempAdaptMutTrunc[tempAdaptMutTrunc['delay'] <= maxDelay]
    stempAdaptMutTrunc = stempAdaptMaster[(stempAdaptMaster.n >= nThresh) & (stempAdaptMaster.combo_id == cIDMut)].drop_duplicates()
    # maxDelay = 500
    # stempAdaptMutTrunc = stempAdaptMutTrunc[stempAdaptMutTrunc['delay'] <= maxDelay]
    # stempAdaptMutTrunc = stempAdaptMutTrunc[stempAdaptMutTrunc['delay'] >= -1*maxDelay]
    tempAdaptTrunc0 = tempAdaptMaster[(tempAdaptMaster.n >= nThresh) & (tempAdaptMaster.combo_id == cID0)].drop_duplicates()
    maxDelay = abs(min(tempAdaptTrunc0['delay']))
    tempAdaptTrunc0 = tempAdaptTrunc0[tempAdaptTrunc0['delay'] <= maxDelay]
    stempAdaptTrunc0 = stempAdaptMaster[(stempAdaptMaster.n >= nThresh) & (stempAdaptMaster.combo_id == cID0)].drop_duplicates()
    stempAdaptTrunc0 = stempAdaptTrunc0[stempAdaptTrunc0['delay'] <= maxDelay]
    stempAdaptTrunc0 = stempAdaptTrunc0[stempAdaptTrunc0['delay'] >= -1*maxDelay]
    ##
    ##
    inner_grid = outer_grid[bcomm-1].subgridspec(2, 6, hspace=0.4, wspace=.6)
    ax1 = fig.add_subplot(inner_grid[0,1:3]) # row 0 with axes spanning 2 cols on odds
    ax2 = fig.add_subplot(inner_grid[0,3:5])
    ax3 = fig.add_subplot(inner_grid[1,0:2]) # row 0 with axes spanning 2 cols on evens
    ax4 = fig.add_subplot(inner_grid[1,2:4])
    ax5 = fig.add_subplot(inner_grid[1,4:])
    # axes = [ax1, ax2, ax3, ax4, ax5,
    #     ax3.twinx(), ax4.twinx(), ax5.twinx()]
    axes = [ax1, ax2, ax3, ax4, ax5]

    axes[0].fill_between(tempAdaptTrunc0['delay'], tempAdaptTrunc0['mean']-tempAdaptTrunc0['std'],
                            tempAdaptTrunc0['mean'] + tempAdaptTrunc0['std'], color=c0, alpha=a)
    axes[0].plot(tempAdaptTrunc0['delay'], tempAdaptTrunc0['mean'], 
                    linewidth=1.5, color=c0, label = r'virus, $\sigma = 0$')
    axes[1].fill_between(tempAdaptTrunc0['delay'], tempAdaptTrunc0['mean']-tempAdaptTrunc0['std'],
                            tempAdaptTrunc0['mean'] + tempAdaptTrunc0['std'], color=c0, alpha=a)
    axes[1].plot(tempAdaptTrunc0['delay'], tempAdaptTrunc0['mean'], 
                    linewidth=1.5, color=c0)
    axes[2].fill_between(tempAdaptTrunc0['delay'], tempAdaptTrunc0['mean']-tempAdaptTrunc0['std'],
                        tempAdaptTrunc0['mean'] + tempAdaptTrunc0['std'], color=c0, alpha=a)
    axes[2].plot(tempAdaptTrunc0['delay'], tempAdaptTrunc0['mean'], 
                    linewidth=1.5, color=c0, label = r'virus, $\sigma = 0$')
    axes[2].fill_between(stempAdaptTrunc0['delay'], stempAdaptTrunc0['mean']-stempAdaptTrunc0['std'],
                        stempAdaptTrunc0['mean'] + stempAdaptTrunc0['std'], color=cHost0, alpha=a)
    axes[2].plot(stempAdaptTrunc0['delay'], stempAdaptTrunc0['mean'], 
                    linewidth=1.5, color=cHost0, label = r'host, $\sigma = 0$')
    #
    axes[0].fill_between(tempAdaptTrunc['delay'], tempAdaptTrunc['mean']-tempAdaptTrunc['std'], 
                        tempAdaptTrunc['mean'] + tempAdaptTrunc['std'], color=cScale, alpha=a)
    axes[0].plot(tempAdaptTrunc['delay'],tempAdaptTrunc['mean'],
                    linewidth=1.5, color=cScale, label=''.join(['virus, ', r'$\sigma =$','{}, '.format(scale), r'$\mu_s = 0$']))
    axes[3].fill_between(tempAdaptTrunc['delay'], tempAdaptTrunc['mean']-tempAdaptTrunc['std'], 
                        tempAdaptTrunc['mean'] + tempAdaptTrunc['std'], color=cScale, alpha=a)
    axes[3].plot(tempAdaptTrunc['delay'],tempAdaptTrunc['mean'],
                    linewidth=1.5, color=cScale, label=''.join(['virus, ', r'$\sigma =$','{}, '.format(scale), r'$\mu_s = 0$']))
    axes[4].fill_between(stempAdaptTrunc['delay'], stempAdaptTrunc['mean']-stempAdaptTrunc['std'], 
                        stempAdaptTrunc['mean'] + stempAdaptTrunc['std'], color=cHostScale, alpha=a)
    axes[4].plot(stempAdaptTrunc['delay'],stempAdaptTrunc['mean'],
                    linewidth=1.5, color=cHostScale, label=''.join(['host, ', r'$\sigma =$','{}, '.format(scale), r'$\mu_s = 0$']))
    #
    axes[1].fill_between(tempAdaptMutTrunc['delay'], tempAdaptMutTrunc['mean']-tempAdaptMutTrunc['std'], 
                        tempAdaptMutTrunc['mean'] + tempAdaptMutTrunc['std'], color=cMut, alpha=a)
    axes[1].plot(tempAdaptMutTrunc['delay'],tempAdaptMutTrunc['mean'],
                    linewidth=1.5, color=cMut, label=''.join(['virus, ', r'$\sigma =$','{}, '.format(scale), r'$\mu_s = 1$']))
    axes[3].fill_between(tempAdaptMutTrunc['delay'], tempAdaptMutTrunc['mean']-tempAdaptMutTrunc['std'], 
                        tempAdaptMutTrunc['mean'] + tempAdaptMutTrunc['std'], color=cMut, alpha=a)
    axes[3].plot(tempAdaptMutTrunc['delay'],tempAdaptMutTrunc['mean'],
                    linewidth=1.5, color=cMut, label=''.join(['virus, ', r'$\sigma =$','{}, '.format(scale), r'$\mu_s = 1$']))
    axes[4].fill_between(stempAdaptMutTrunc['delay'], stempAdaptMutTrunc['mean']-stempAdaptMutTrunc['std'], 
                        stempAdaptMutTrunc['mean'] + stempAdaptMutTrunc['std'], color=cHostMut, alpha=a)
    axes[4].plot(stempAdaptMutTrunc['delay'],stempAdaptMutTrunc['mean'],
                    linewidth=1.5, color=cHostMut, label=''.join(['host, ', r'$\sigma =$','{}, '.format(scale), r'$\mu_s = 1$']))
    #
    axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(100))
    axes[1].xaxis.set_minor_locator(ticker.MultipleLocator(100))
    axes[2].xaxis.set_minor_locator(ticker.MultipleLocator(100))
    axes[3].xaxis.set_minor_locator(ticker.MultipleLocator(100))
    axes[4].xaxis.set_minor_locator(ticker.MultipleLocator(100))

    axes[0].tick_params(axis='x', labelsize=15)
    axes[0].tick_params(axis='y', labelsize=15)
    axes[1].tick_params(axis='x', labelsize=15)
    axes[1].tick_params(axis='y', labelsize=15)
    axes[2].tick_params(axis='x', labelsize=15)
    axes[2].tick_params(axis='y', labelsize=15)
    axes[3].tick_params(axis='x', labelsize=15)
    axes[3].tick_params(axis='y', labelsize=15)
    axes[4].tick_params(axis='x', labelsize=15)
    axes[4].tick_params(axis='y', labelsize=15)

    axes[0].legend(loc='upper right', fontsize=15)
    axes[1].legend(loc='upper right', fontsize=15)
    axes[2].legend(loc='upper right', fontsize=15)
    axes[3].legend(loc='upper right', fontsize=15)
    axes[4].legend(loc='upper left', fontsize=15)
    axes[0].set_xlabel(xlabel = r'Delay $\tau$',fontsize=15,labelpad=15)
    axes[0].set_ylabel(ylabel ='Temporal\nAdaptation (TA)',labelpad=15,fontsize=15)
    axes[2].set_ylabel(ylabel ='Temporal\nAdaptation (TA)',labelpad=15,fontsize=15)
    
fig.savefig(os.path.join('/Volumes/Yadgah/sylvain-martin-collab',simDir,'supp-SUS-TA-all_{}.pdf'.format(nThresh)),dpi=resolve)
plt.close('all')


####################
####################
#################### Normalized By Contact Probabili
####################
####################
####################


tempAdaptInfMaster = pd.DataFrame({'combo_id': [], 'delay': [], 'mean': [], 'std': [], 'n' : []})
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
    if cID0 not in tempAdaptInfMaster['combo_id'].values:
        tempAdaptInf = pd.read_sql_query("SELECT t, delay, vfrequency, bfrequency, run_id \
                                    FROM vtemporal_adaptation \
                                    WHERE run_id in ({})".format(', '.join(map(str, runIDs0['run_id'].values))), conTemp)\
                                    .merge(runIDs0,on=['run_id'])
        tempAdaptInf = tempAdaptInf[tempAdaptInf['t'] > 10]
        tempAdaptInf['TA'] = tempAdaptInf['vfrequency'] * tempAdaptInf['bfrequency']
        tempAdaptInf = tempAdaptInf.groupby(['t', 'delay', 'run_id'])\
            .agg(TA=('TA', 'sum')).reset_index().merge(tempAdaptInf.drop(columns=['TA']),on=['t', 'delay', 'run_id'])
        tempAdaptInf['TA'] =  (np.array(tempAdaptInf['vfrequency']) * np.array(tempAdaptInf['bfrequency'])  / np.array(tempAdaptInf['TA'])) * np.array(tempAdaptInf['bfrequency'])
        tempAdaptInf = tempAdaptInf.groupby(['t','delay','run_id', 'combo_id'])\
            .agg(TA=('TA', 'sum')).reset_index()
        tempAdaptInf = tempAdaptInf.groupby(['delay','run_id','combo_id'])\
            .agg(mean=('TA', 'mean'),std=('TA','std')).reset_index()\
            .dropna()
        tempAdaptInf = tempAdaptInf.groupby(['combo_id', 'delay'])\
            .agg(mean=('mean', 'mean'),std=('mean','std'),n=('run_id','size')).reset_index()\
            .dropna()
        tempAdaptInfMaster = pd.concat([tempAdaptInfMaster, tempAdaptInf], ignore_index=True)
    ##
    ##
    tempAdaptInf = pd.read_sql_query("SELECT t, delay, vfrequency, bfrequency, run_id \
                                FROM vtemporal_adaptation \
                                WHERE run_id in ({})".format(', '.join(map(str, runIDs['run_id'].values))), conTemp)\
                                .merge(runIDs,on=['run_id'])
    tempAdaptInf = tempAdaptInf[tempAdaptInf['t'] > 10]
    tempAdaptInf['TA'] = tempAdaptInf['vfrequency'] * tempAdaptInf['bfrequency']
    tempAdaptInf = tempAdaptInf.groupby(['t', 'delay', 'run_id'])\
        .agg(TA=('TA', 'sum')).reset_index().merge(tempAdaptInf.drop(columns=['TA']),on=['t', 'delay', 'run_id'])
    tempAdaptInf['TA'] =  (np.array(tempAdaptInf['vfrequency']) * np.array(tempAdaptInf['bfrequency'])  / np.array(tempAdaptInf['TA'])) * np.array(tempAdaptInf['bfrequency'])
    tempAdaptInf = tempAdaptInf.groupby(['t','delay','run_id', 'combo_id'])\
        .agg(TA=('TA', 'sum')).reset_index()
    tempAdaptInf = tempAdaptInf.groupby(['delay','run_id', 'combo_id'])\
        .agg(mean=('TA', 'mean'),std=('TA','std')).reset_index()\
        .dropna()
    tempAdaptInf = tempAdaptInf.groupby(['combo_id', 'delay'])\
        .agg(mean=('mean', 'mean'),std=('mean','std'),n=('run_id','size')).reset_index()\
        .dropna()
    tempAdaptInfMaster = pd.concat([tempAdaptInfMaster, tempAdaptInf], ignore_index=True)
    ## 
    ##

tempAdaptInfMaster = tempAdaptInfMaster.drop_duplicates()

#######
## HOST TEMPORAL ADAPTATION (WEIGHTED BY CONTACT PROBABILITIES)
#######

htempAdaptInfMaster = pd.DataFrame({'combo_id': [], 'delay': [], 'mean': [], 'std': [], 'n' : []})
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
    if cID0 not in htempAdaptInfMaster['combo_id'].values:
        tempAdaptInf = pd.read_sql_query("SELECT t, delay, vfrequency, bfrequency, run_id \
                                    FROM btemporal_adaptation \
                                    WHERE run_id in ({})".format(', '.join(map(str, runIDs0['run_id'].values))), conHTemp)\
                                    .merge(runIDs0,on=['run_id']).drop_duplicates()
        tempAdaptInf = tempAdaptInf[tempAdaptInf['t'] > 10]
        tempAdaptInf['TA'] = tempAdaptInf['vfrequency'] * tempAdaptInf['bfrequency']
        tempAdaptInf = tempAdaptInf.groupby(['t', 'delay', 'run_id'])\
            .agg(TA=('TA', 'sum')).reset_index().merge(tempAdaptInf.drop(columns=['TA']),on=['t', 'delay', 'run_id'])
        tempAdaptInf['TA'] =  (np.array(tempAdaptInf['vfrequency']) * np.array(tempAdaptInf['bfrequency'])  / np.array(tempAdaptInf['TA'])) * np.array(tempAdaptInf['vfrequency'])
        tempAdaptInf = tempAdaptInf.groupby(['t','delay','run_id', 'combo_id'])\
            .agg(TA=('TA', 'sum')).reset_index()
        tempAdaptInf = tempAdaptInf.groupby(['delay','run_id','combo_id'])\
            .agg(mean=('TA', 'mean'),std=('TA','std')).reset_index()\
            .dropna()
        tempAdaptInf = tempAdaptInf.groupby(['combo_id', 'delay'])\
            .agg(mean=('mean', 'mean'),std=('mean','std'),n=('run_id','size')).reset_index()\
            .dropna()
        htempAdaptInfMaster = pd.concat([htempAdaptInfMaster, tempAdaptInf], ignore_index=True)
    ##
    ##
    tempAdaptInf = pd.read_sql_query("SELECT t, delay, vfrequency, bfrequency, run_id \
                                FROM btemporal_adaptation \
                                WHERE run_id in ({})".format(', '.join(map(str, runIDs['run_id'].values))), conHTemp)\
                                .merge(runIDs,on=['run_id']).drop_duplicates()
    tempAdaptInf = tempAdaptInf[tempAdaptInf['t'] > 10]
    tempAdaptInf['TA'] = tempAdaptInf['vfrequency'] * tempAdaptInf['bfrequency']
    tempAdaptInf = tempAdaptInf.groupby(['t', 'delay', 'run_id'])\
        .agg(TA=('TA', 'sum')).reset_index().merge(tempAdaptInf.drop(columns=['TA']),on=['t', 'delay', 'run_id'])
    tempAdaptInf['TA'] =  (np.array(tempAdaptInf['vfrequency']) * np.array(tempAdaptInf['bfrequency'])  / np.array(tempAdaptInf['TA'])) * np.array(tempAdaptInf['vfrequency'])
    tempAdaptInf = tempAdaptInf.groupby(['t','delay','run_id', 'combo_id'])\
        .agg(TA=('TA', 'sum')).reset_index()
    tempAdaptInf = tempAdaptInf.groupby(['delay','run_id', 'combo_id'])\
        .agg(mean=('TA', 'mean'),std=('TA','std')).reset_index()\
        .dropna()
    tempAdaptInf = tempAdaptInf.groupby(['combo_id', 'delay'])\
        .agg(mean=('mean', 'mean'),std=('mean','std'),n=('run_id','size')).reset_index()\
        .dropna()
    htempAdaptInfMaster = pd.concat([htempAdaptInfMaster, tempAdaptInf], ignore_index=True)
    ##

htempAdaptInfMaster = htempAdaptInfMaster.drop_duplicates()

scale = 10
resolve = 500
c0 = 'gold'
cScale = 'darkblue'
cMut = 'darkred'
cHost0 = 'darkgreen'
cHostScale = 'darkblue'
cHostMut = 'darkred'
a = 0.3
nThresh = 10 
tempxmin = 0
tempxmax = 0
figsize = (17,22)
fig = plt.figure(figsize=figsize)
outer_grid = fig.add_gridspec(2, 1, wspace=1.2, hspace=.3)
for bcomm in [1,2]:
    if bcomm == 1:
        divtype = 'Monomorphic'
    if bcomm == 2:
        divtype = 'Polymorphic'
    comboSpace = pd.read_sql_query(
        "SELECT combo_id \
            FROM param_combos WHERE \
            microbe_mutation_prob_per_spacer = {0} \
            AND init_bcomm_function = {1} \
            AND n_particles_per_vstrain > 0 \
            AND evofunctionScale = {2} \
            ORDER BY combo_id"
        .format(0, bcomm, scale),
        conSim)
    cID = comboSpace['combo_id'].values[0]
    comboSpaceMut = pd.read_sql_query(
        "SELECT combo_id \
            FROM param_combos WHERE \
            microbe_mutation_prob_per_spacer = {0} \
            AND init_bcomm_function = {1} \
            AND n_particles_per_vstrain > 0 \
            AND evofunctionScale = {2} \
            ORDER BY combo_id"
        .format(1, bcomm, scale),
        conSim)
    cIDMut = comboSpaceMut['combo_id'].values[0]
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
    cID0 = comboSpace0['combo_id'].values[0]
    ##
    ## 
    tempAdaptInfTrunc = tempAdaptInfMaster[(tempAdaptInfMaster.n >= nThresh) & (tempAdaptInfMaster.combo_id == cID)].drop_duplicates()
    maxDelay = abs(min(tempAdaptInfTrunc['delay']))
    tempAdaptInfTrunc = tempAdaptInfTrunc[tempAdaptInfTrunc['delay'] <= maxDelay]
    htempAdaptInfTrunc = htempAdaptInfMaster[(htempAdaptInfMaster.n >= nThresh) & (htempAdaptInfMaster.combo_id == cID)].drop_duplicates()
    # maxDelay = 500
    # htempAdaptInfTrunc = htempAdaptInfTrunc[htempAdaptInfTrunc['delay'] <= maxDelay]
    # htempAdaptInfTrunc = htempAdaptInfTrunc[htempAdaptInfTrunc['delay'] >= -1*maxDelay]
    tempAdaptInfMutTrunc = tempAdaptInfMaster[(tempAdaptInfMaster.n >= nThresh) & (tempAdaptInfMaster.combo_id == cIDMut)].drop_duplicates()
    maxDelay = abs(min(tempAdaptInfMutTrunc['delay']))
    tempAdaptInfMutTrunc = tempAdaptInfMutTrunc[tempAdaptInfMutTrunc['delay'] <= maxDelay]
    htempAdaptInfMutTrunc = htempAdaptInfMaster[(htempAdaptInfMaster.n >= nThresh) & (htempAdaptInfMaster.combo_id == cIDMut)].drop_duplicates()
    # maxDelay = 500
    # htempAdaptInfMutTrunc = htempAdaptInfMutTrunc[htempAdaptInfMutTrunc['delay'] <= maxDelay]
    # htempAdaptInfMutTrunc = htempAdaptInfMutTrunc[htempAdaptInfMutTrunc['delay'] >= -1*maxDelay]
    tempAdaptInfTrunc0 = tempAdaptInfMaster[(tempAdaptInfMaster.n >= nThresh) & (tempAdaptInfMaster.combo_id == cID0)].drop_duplicates()
    maxDelay = abs(min(tempAdaptInfTrunc0['delay']))
    tempAdaptInfTrunc0 = tempAdaptInfTrunc0[tempAdaptInfTrunc0['delay'] <= maxDelay]
    htempAdaptInfTrunc0 = htempAdaptInfMaster[(htempAdaptInfMaster.n >= nThresh) & (htempAdaptInfMaster.combo_id == cID0)].drop_duplicates()
    htempAdaptInfTrunc0 = htempAdaptInfTrunc0[htempAdaptInfTrunc0['delay'] <= maxDelay]
    htempAdaptInfTrunc0 = htempAdaptInfTrunc0[htempAdaptInfTrunc0['delay'] >= -1*maxDelay]
    ##
    ##
    inner_grid = outer_grid[bcomm-1].subgridspec(2, 6, hspace=0.4, wspace=.6)
    ax1 = fig.add_subplot(inner_grid[0,1:3]) # row 0 with axes spanning 2 cols on odds
    ax2 = fig.add_subplot(inner_grid[0,3:5])
    ax3 = fig.add_subplot(inner_grid[1,0:2]) # row 0 with axes spanning 2 cols on evens
    ax4 = fig.add_subplot(inner_grid[1,2:4])
    ax5 = fig.add_subplot(inner_grid[1,4:])
    # axes = [ax1, ax2, ax3, ax4, ax5,
    #     ax3.twinx(), ax4.twinx(), ax5.twinx()]
    axes = [ax1, ax2, ax3, ax4, ax5]

    axes[0].fill_between(tempAdaptInfTrunc0['delay'], tempAdaptInfTrunc0['mean']-tempAdaptInfTrunc0['std'],
                            tempAdaptInfTrunc0['mean'] + tempAdaptInfTrunc0['std'], color=c0, alpha=a)
    axes[0].plot(tempAdaptInfTrunc0['delay'], tempAdaptInfTrunc0['mean'], 
                    linewidth=1.5, color=c0, label = r'virus, $\sigma = 0$')
    axes[1].fill_between(tempAdaptInfTrunc0['delay'], tempAdaptInfTrunc0['mean']-tempAdaptInfTrunc0['std'],
                            tempAdaptInfTrunc0['mean'] + tempAdaptInfTrunc0['std'], color=c0, alpha=a)
    axes[1].plot(tempAdaptInfTrunc0['delay'], tempAdaptInfTrunc0['mean'], 
                    linewidth=1.5, color=c0)
    axes[2].fill_between(tempAdaptInfTrunc0['delay'], tempAdaptInfTrunc0['mean']-tempAdaptInfTrunc0['std'],
                        tempAdaptInfTrunc0['mean'] + tempAdaptInfTrunc0['std'], color=c0, alpha=a)
    axes[2].plot(tempAdaptInfTrunc0['delay'], tempAdaptInfTrunc0['mean'], 
                    linewidth=1.5, color=c0, label = r'virus, $\sigma = 0$')
    axes[2].fill_between(htempAdaptInfTrunc0['delay'], htempAdaptInfTrunc0['mean']-htempAdaptInfTrunc0['std'],
                        htempAdaptInfTrunc0['mean'] + htempAdaptInfTrunc0['std'], color=cHost0, alpha=a)
    axes[2].plot(htempAdaptInfTrunc0['delay'], htempAdaptInfTrunc0['mean'], 
                    linewidth=1.5, color=cHost0, label = r'host, $\sigma = 0$')
    #
    axes[0].fill_between(tempAdaptInfTrunc['delay'], tempAdaptInfTrunc['mean']-tempAdaptInfTrunc['std'], 
                        tempAdaptInfTrunc['mean'] + tempAdaptInfTrunc['std'], color=cScale, alpha=a)
    axes[0].plot(tempAdaptInfTrunc['delay'],tempAdaptInfTrunc['mean'],
                    linewidth=1.5, color=cScale, label=''.join(['virus, ', r'$\sigma =$','{}, '.format(scale), r'$\mu_s = 0$']))
    axes[3].fill_between(tempAdaptInfTrunc['delay'], tempAdaptInfTrunc['mean']-tempAdaptInfTrunc['std'], 
                        tempAdaptInfTrunc['mean'] + tempAdaptInfTrunc['std'], color=cScale, alpha=a)
    axes[3].plot(tempAdaptInfTrunc['delay'],tempAdaptInfTrunc['mean'],
                    linewidth=1.5, color=cScale, label=''.join(['virus, ', r'$\sigma =$','{}, '.format(scale), r'$\mu_s = 0$']))
    axes[4].fill_between(htempAdaptInfTrunc['delay'], htempAdaptInfTrunc['mean']-htempAdaptInfTrunc['std'], 
                        htempAdaptInfTrunc['mean'] + htempAdaptInfTrunc['std'], color=cHostScale, alpha=a)
    axes[4].plot(htempAdaptInfTrunc['delay'],htempAdaptInfTrunc['mean'],
                    linewidth=1.5, color=cHostScale, label=''.join(['host, ', r'$\sigma =$','{}, '.format(scale), r'$\mu_s = 0$']))
    #
    axes[1].fill_between(tempAdaptInfMutTrunc['delay'], tempAdaptInfMutTrunc['mean']-tempAdaptInfMutTrunc['std'], 
                        tempAdaptInfMutTrunc['mean'] + tempAdaptInfMutTrunc['std'], color=cMut, alpha=a)
    axes[1].plot(tempAdaptInfMutTrunc['delay'],tempAdaptInfMutTrunc['mean'],
                    linewidth=1.5, color=cMut, label=''.join(['virus, ', r'$\sigma =$','{}, '.format(scale), r'$\mu_s = 1$']))
    axes[3].fill_between(tempAdaptInfMutTrunc['delay'], tempAdaptInfMutTrunc['mean']-tempAdaptInfMutTrunc['std'], 
                        tempAdaptInfMutTrunc['mean'] + tempAdaptInfMutTrunc['std'], color=cMut, alpha=a)
    axes[3].plot(tempAdaptInfMutTrunc['delay'],tempAdaptInfMutTrunc['mean'],
                    linewidth=1.5, color=cMut, label=''.join(['virus, ', r'$\sigma =$','{}, '.format(scale), r'$\mu_s = 1$']))
    axes[4].fill_between(htempAdaptInfMutTrunc['delay'], htempAdaptInfMutTrunc['mean']-htempAdaptInfMutTrunc['std'], 
                        htempAdaptInfMutTrunc['mean'] + htempAdaptInfMutTrunc['std'], color=cHostMut, alpha=a)
    axes[4].plot(htempAdaptInfMutTrunc['delay'],htempAdaptInfMutTrunc['mean'],
                    linewidth=1.5, color=cHostMut, label=''.join(['host, ', r'$\sigma =$','{}, '.format(scale), r'$\mu_s = 1$']))
    #
    axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(100))
    axes[1].xaxis.set_minor_locator(ticker.MultipleLocator(100))
    axes[2].xaxis.set_minor_locator(ticker.MultipleLocator(100))
    axes[3].xaxis.set_minor_locator(ticker.MultipleLocator(100))
    axes[4].xaxis.set_minor_locator(ticker.MultipleLocator(100))

    axes[0].tick_params(axis='x', labelsize=15)
    axes[0].tick_params(axis='y', labelsize=15)
    axes[1].tick_params(axis='x', labelsize=15)
    axes[1].tick_params(axis='y', labelsize=15)
    axes[2].tick_params(axis='x', labelsize=15)
    axes[2].tick_params(axis='y', labelsize=15)
    axes[3].tick_params(axis='x', labelsize=15)
    axes[3].tick_params(axis='y', labelsize=15)
    axes[4].tick_params(axis='x', labelsize=15)
    axes[4].tick_params(axis='y', labelsize=15)

    axes[0].legend(loc='upper right', fontsize=15)
    axes[1].legend(loc='upper right', fontsize=15)
    axes[2].legend(loc='lower left', fontsize=15)
    axes[3].legend(loc='upper right', fontsize=15)
    axes[4].legend(loc='lower left', fontsize=15)
    axes[0].set_xlabel(xlabel = r'Delay $\tau$',fontsize=15,labelpad=15)
    axes[0].set_ylabel(ylabel ='Temporal\nAdaptation (TA)',labelpad=15,fontsize=15)
    axes[2].set_ylabel(ylabel ='Temporal\nAdaptation (TA)',labelpad=15,fontsize=15)
    

fig.savefig(os.path.join('/Volumes/Yadgah/sylvain-martin-collab',simDir,'supp-CONTACT-PROBS-TA-all_{}.pdf'.format(nThresh)),dpi=resolve)
plt.close('all')
