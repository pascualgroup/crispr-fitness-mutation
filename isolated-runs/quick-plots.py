#!/usr/bin/env python3
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
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.cm as cm
from matplotlib.colors import Normalize
# import mpl_scatter_density

SCRIPT_PATH = os.path.abspath(os.path.dirname(__file__))  # cluster

# dir = 'sylvain-martin-collab/trial-3/gathered-analyses'
# DBWALLS_PATH = os.path.join('/Volumes/Yadgah/{0}/walls-shannon/mean-walls-shannon.sqlite'.format(dir))
# DBEXT_PATH = os.path.join('/Volumes/Yadgah/{0}/extinctions/mean-extinctions.sqlite'.format(dir))


# DBSIM_PATH = os.path.join('/Volumes/Yadgah/sylvain-martin-collab/8_MOI3/sweep_db.sqlite')
# DBEXT_PATH = os.path.join('/Volumes/Yadgah/sylvain-martin-collab/8_MOI3/gathered-analyses/extinctions/extinctions.sqlite')
# DBESC_PATH = os.path.join('/Volumes/Yadgah/sylvain-martin-collab/8_MOI3/gathered-analyses/escape-count/escape-count.sqlite')
DBSIM_PATH = os.path.join('/Volumes/Yadgah/sylvain-martin-collab/12_MOI3/sweep_db_gathered.sqlite')
# DBSIM_PATH = os.path.join('/Volumes/Yadgah/sylvain-martin-collab/12_MOI3/sweep_db.sqlite')
DBTRAIT_PATH = os.path.join('/Volumes/Yadgah/sylvain-martin-collab/12_MOI3/gathered-analyses/trait-evolution/trait-evolution.sqlite')
DBSIM_PATHnull = os.path.join('/Volumes/Yadgah/sylvain-martin-collab/14_MOI3/sweep_db_gathered.sqlite')
DBSIM_PATHnoV = os.path.join('/Volumes/Yadgah/sylvain-martin-collab/13_MOI3/sweep_db_gathered.sqlite')
# DBSIM_PATHnull = os.path.join('/Volumes/Yadgah/sylvain-martin-collab/14_MOI3/sweep_db.sqlite')
DBTRAIT_PATHnull = os.path.join('/Volumes/Yadgah/sylvain-martin-collab/14_MOI3/gathered-analyses/trait-evolution/trait-evolution.sqlite')
DBTRAIT_PATHnoV = os.path.join('/Volumes/Yadgah/sylvain-martin-collab/13_MOI3/gathered-analyses/trait-evolution/trait-evolution.sqlite')


conSim = sqlite3.connect(DBSIM_PATH)
curSim = conSim.cursor()
conTrait = sqlite3.connect(DBTRAIT_PATH)
curTrait = conTrait.cursor()
conSimNull = sqlite3.connect(DBSIM_PATHnull)
curSimNull = conSimNull.cursor()
conSimNoV = sqlite3.connect(DBSIM_PATHnoV)
curSimNoV = conSimNoV.cursor()
conTraitNull = sqlite3.connect(DBTRAIT_PATHnull)
curTraitNull = conTraitNull.cursor()
conTraitNoV = sqlite3.connect(DBTRAIT_PATHnoV)
curTraitNoV = conTraitNoV.cursor()

######
######
###### Distribution of Fitness
######
######
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
        ORDER BY combo_id"
    .format(micMutProb, viralMutProb, viralDecay, failureProb, carryCap, micDeath),
    conSim)
comboSpaceNoV = pd.read_sql_query(
    "SELECT combo_id, evofunctionScale \
FROM param_combos WHERE microbe_mutation_prob = {0} AND viral_mutation_rate = {1} \
        AND viral_decay_rate = {2} AND crispr_failure_prob = {3} AND microbe_carrying_capacity = {4} \
        AND microbe_death_rate = {5} \
        ORDER BY combo_id"
    .format(micMutProb, viralMutProb, viralDecay, failureProb, carryCap, micDeath),
    conTraitNoV)
comboSpaceNull = pd.read_sql_query(
    "SELECT combo_id, evofunctionScale \
FROM param_combos WHERE microbe_mutation_prob = {0} AND viral_mutation_rate = {1} \
        AND viral_decay_rate = {2} AND crispr_failure_prob = {3} AND microbe_carrying_capacity = {4} \
        AND microbe_death_rate = {5} \
        ORDER BY combo_id"
    .format(micMutProb, viralMutProb, viralDecay, failureProb, carryCap, micDeath),
    conTraitNull)
cIDNull = 3
cIDNoV = 7
cID = 27
runIDs = pd.read_sql_query(
    "SELECT combo_id, run_id FROM runs \
        WHERE combo_id in ({})"
    .format(', '.join(map(str, comboSpace['combo_id']))), conSim)
runIDsNoV = pd.read_sql_query(
    "SELECT combo_id, run_id FROM runs \
        WHERE combo_id in ({})"
    .format(cIDNoV), conTraitNoV)
runIDsNull = pd.read_sql_query(
    "SELECT combo_id, run_id FROM runs \
        WHERE combo_id in ({})"
    .format(cIDNull), conTraitNull)
##########
noVDist = pd.read_sql_query("SELECT t, expectation, run_id FROM b_intrinsic_fitness_expectation \
                           WHERE run_id in ({0})"
                              .format(', '.join(map(str, runIDsNoV['run_id']))), conTraitNoV)\
    .merge(pd.read_sql_query("SELECT t, variance, run_id FROM b_intrinsic_fitness_variance \
                           WHERE run_id in ({0})"
                             .format(', '.join(map(str, runIDsNoV['run_id']))), conTraitNoV), on=['t', 'run_id'])
noVDist['variance'] = np.array(noVDist['variance']**(1/2))
noVDist = noVDist.rename(columns={'variance': 'std'})
times = list(range(0,2001,1))
for runID in np.unique(noVDist['run_id']):
    # print('runID: {}'.format(runID))
    mTimes = list(set(times) - set(noVDist[noVDist.run_id == runID]['t']))
    mTimes.sort()
    timeInt = []
    for i in range(0,len(mTimes),1):
        # print('time_idx: {}'.format(i))
        if i == 0:
            timeInt.append(mTimes[i])
        elif i == len(mTimes)-1:
            if mTimes[i] - mTimes[i-1] == 1:
                timeInt.append(mTimes[i])
                t = min(timeInt)
                missing = noVDist[(noVDist.run_id == runID)
                                   & (noVDist.t <= t)].copy()
                missing = noVDist[(noVDist.run_id == runID) & (
                    noVDist.t == max(missing['t']))].copy()
                missing = pd.DataFrame(
                    {'t': timeInt, 'expectation': len(timeInt)*list(missing['expectation'].values),
                     'run_id': len(timeInt)*list(missing['run_id'].values),
                     'std': len(timeInt)*list(missing['std'].values)})
                noVDist = pd.concat([noVDist, missing], ignore_index=True)
                timeInt = [mTimes[i]]
            else: 
                t = min(mTimes[i])
                missing = noVDist[(noVDist.run_id==runID)&(noVDist.t <= t)].copy()
                missing = noVDist[(noVDist.run_id==runID)&(noVDist.t == max(missing['t']))].copy()
                missing = pd.DataFrame(
                    {'t': timeInt, 'expectation': len(timeInt)*list(missing['expectation'].values),
                     'run_id': len(timeInt)*list(missing['run_id'].values),
                     'std': len(timeInt)*list(missing['std'].values)})
                noVDist = pd.concat([noVDist,missing],ignore_index=True)
                timeInt = [mTimes[i]]
        else:
            if mTimes[i] - mTimes[i-1] == 1:
                timeInt.append(mTimes[i])
            else: 
                t = min(timeInt)
                missing = noVDist[(noVDist.run_id==runID)&(noVDist.t <= t)].copy()
                missing = noVDist[(noVDist.run_id==runID)&(noVDist.t == max(missing['t']))].copy()
                missing = pd.DataFrame(
                    {'t': timeInt, 'expectation': len(timeInt)*list(missing['expectation'].values),
                     'run_id': len(timeInt)*list(missing['run_id'].values),
                     'std': len(timeInt)*list(missing['std'].values)})
                noVDist = pd.concat([noVDist,missing],ignore_index=True)
                timeInt = [mTimes[i]]
    
noVDist = noVDist.sort_values(by=['run_id','t'])
noVDist = noVDist.merge(runIDsNoV, on=['run_id']).merge(comboSpaceNoV, on=['combo_id'])\
    .sort_values(by=['evofunctionScale', 'run_id'])
meanNoVDist = noVDist.groupby(['combo_id', 'evofunctionScale','t'])\
    .agg(meanExp=('expectation', 'mean'), stdExp=('expectation', 'std'), \
         meanStd=('std', 'mean'), stdStd=('std', 'std')).reset_index()
numRuns = noVDist.groupby(['t','combo_id', 'evofunctionScale'])\
    .agg(n=('expectation', 'size')).reset_index()
meanNoVDist = meanNoVDist.merge(numRuns, on=['combo_id', 'evofunctionScale','t'])
#
dist = pd.read_sql_query("SELECT t, expectation, run_id FROM b_intrinsic_fitness_expectation \
                           WHERE run_id in ({0})"
                         .format(', '.join(map(str, runIDs[runIDs.combo_id==cID]['run_id']))), conTrait)\
    .merge(pd.read_sql_query("SELECT t, variance, run_id FROM b_intrinsic_fitness_variance \
                           WHERE run_id in ({0})"
                             .format(', '.join(map(str, runIDs[runIDs.combo_id == cID]['run_id']))), conTrait), on=['t', 'run_id'])
dist['variance'] = np.array(dist['variance']**(1/2))
dist = dist.rename(columns={'variance': 'std'})
times = list(range(0, 2001, 1))
for runID in np.unique(dist['run_id']):
    # print('runID: {}'.format(runID))
    mTimes = list(set(times) - set(dist[dist.run_id == runID]['t']))
    mTimes.sort()
    timeInt = []
    for i in range(0, len(mTimes), 1):
        # print('time_idx: {}'.format(i))
        if i == 0:
            timeInt.append(mTimes[i])
        elif i == len(mTimes)-1:
            if mTimes[i] - mTimes[i-1] == 1:
                timeInt.append(mTimes[i])
                t = min(timeInt)
                missing = dist[(dist.run_id == runID)
                               & (dist.t <= t)].copy()
                missing = dist[(dist.run_id == runID) & (
                    dist.t == max(missing['t']))].copy()
                missing = pd.DataFrame(
                    {'t': timeInt, 'expectation': len(timeInt)*list(missing['expectation'].values),
                     'run_id': len(timeInt)*list(missing['run_id'].values),
                     'std': len(timeInt)*list(missing['std'].values)})
                dist = pd.concat([dist, missing], ignore_index=True)
                timeInt = [mTimes[i]]
            else:
                t = min(mTimes[i])
                missing = dist[(dist.run_id == runID) & (dist.t <= t)].copy()
                missing = dist[(dist.run_id == runID) & (
                    dist.t == max(missing['t']))].copy()
                missing = pd.DataFrame(
                    {'t': timeInt, 'expectation': len(timeInt)*list(missing['expectation'].values),
                     'run_id': len(timeInt)*list(missing['run_id'].values),
                     'std': len(timeInt)*list(missing['std'].values)})
                dist = pd.concat([dist, missing], ignore_index=True)
                timeInt = [mTimes[i]]
        else:
            if mTimes[i] - mTimes[i-1] == 1:
                timeInt.append(mTimes[i])
            else:
                t = min(timeInt)
                missing = dist[(dist.run_id == runID) & (dist.t <= t)].copy()
                missing = dist[(dist.run_id == runID) & (
                    dist.t == max(missing['t']))].copy()
                missing = pd.DataFrame(
                    {'t': timeInt, 'expectation': len(timeInt)*list(missing['expectation'].values),
                     'run_id': len(timeInt)*list(missing['run_id'].values),
                     'std': len(timeInt)*list(missing['std'].values)})
                dist = pd.concat([dist, missing], ignore_index=True)
                timeInt = [mTimes[i]]

dist = dist.sort_values(by=['run_id', 't'])
dist = dist.merge(runIDs, on=['run_id']).merge(comboSpace, on=['combo_id'])\
    .sort_values(by=['evofunctionScale', 'run_id'])
meanDist = dist.groupby(['combo_id', 'evofunctionScale', 't'])\
    .agg(meanExp=('expectation', 'mean'), stdExp=('expectation', 'std'),
         meanStd=('std', 'mean'), stdStd=('std', 'std')).reset_index()
numRuns = dist.groupby(['t', 'combo_id', 'evofunctionScale'])\
    .agg(n=('expectation', 'size')).reset_index()
meanDist = meanDist.merge(numRuns, on=['combo_id', 'evofunctionScale', 't'])
######
######
###### Clade Immune Proportions
######
######
cladeImmune = pd.read_sql_query("SELECT bstrain_id, growth_rate, run_id FROM bgrowthrates \
                           WHERE run_id in ({0})"
                             .format(', '.join(map(str, runIDs[runIDs.combo_id == cID]['run_id']))), conSim)\
                             .merge(\
                               pd.read_sql_query("SELECT t, bstrain_id, abundance, run_id FROM babundance \
                                                    WHERE run_id in ({0})"
                                                .format(', '.join(\
                                                    map(str, runIDs[runIDs.combo_id == cID]['run_id'])\
                                                        )), conSim), on = ['run_id','bstrain_id'] 
                             )\
                             .merge(\
                               pd.read_sql_query("SELECT t, microbial_abundance, run_id FROM summary \
                                                    WHERE run_id in ({0})"
                                                .format(', '.join(\
                                                    map(str, runIDs[runIDs.combo_id == cID]['run_id'])\
                                                        )), conSim), on = ['run_id','t'] 
                             )\
                             .rename(columns={'microbial_abundance':'btotal'})                            
cladeImmune['bfreq'] = cladeImmune['abundance']/cladeImmune['btotal']
cladeImmune = cladeImmune.drop(columns=['abundance','btotal'])
cladeImmune['shannon'] = np.exp(-1*np.array(cladeImmune['bfreq'])*np.log(cladeImmune['bfreq']))
cladeImmune['simpson'] = np.array(cladeImmune['bfreq'])**2
shannon = cladeImmune.groupby(['t','run_id']).agg(stotal=('shannon','prod')).reset_index()
simpson = cladeImmune.groupby(['t','run_id']).agg(ltotal=('simpson','sum')).reset_index()
simpson['ltotal'] = 1/simpson['ltotal']
richness = cladeImmune[['run_id','t','bstrain_id']].drop_duplicates().groupby(['t', 'run_id']).agg(
    rtotal=('bstrain_id', 'size')).reset_index()
shannonStats = shannon.groupby(['t']).agg(
    shanMean=('stotal', 'mean'), shanStd=('stotal', 'std')).reset_index()
simpsonStats = simpson.groupby(['t']).agg(
    simpMean=('ltotal', 'mean'), simpStd=('ltotal', 'std')).reset_index()
richnessStats = richness.groupby(['t']).agg(
    richMean=('rtotal', 'mean'), richStd=('rtotal', 'std')).reset_index()

divNoV = pd.read_sql_query("SELECT t, bstrain_id, abundance, run_id FROM babundance \
                        WHERE run_id in ({0})"
                    .format(', '.join(
                        map(str,
                            runIDsNoV[runIDsNoV.combo_id == cIDNoV]['run_id'])
                            )), conSimNoV)\
    .merge(
    pd.read_sql_query("SELECT t, microbial_abundance, run_id FROM summary \
                                                    WHERE run_id in ({0})"
                      .format(', '.join(
                          map(str, runIDsNoV[runIDsNoV.combo_id == cIDNoV]['run_id'])
                      )), conSimNoV), on=['run_id', 't'])\
    .rename(columns={'microbial_abundance': 'btotal'})
divNoV['bfreq'] = divNoV['abundance']/divNoV['btotal']
divNoV['shannon'] = np.exp(-1*np.array(divNoV['bfreq'])*np.log(divNoV['bfreq']))  
divNoV['simpson'] = np.array(divNoV['bfreq'])**2
divNoV = divNoV.groupby(['t','run_id'])\
            .agg(shannon=('shannon','prod'),
                 simpson=('simpson','sum'),
                 richness=('bstrain_id','size')).reset_index()
divNoV['simpson'] = 1/divNoV['simpson']
divNoV = divNoV.groupby(['t'])\
            .agg(shanMean=('shannon','mean'),
                shanStd=('shannon','std'),
                simpMean=('simpson','mean'),
                simpStd=('simpson','std'),
                richMean=('richness','mean'),
                 richStd=('richness', 'std')).reset_index()
#
cladeImmune = cladeImmune.groupby(['run_id','t','growth_rate'])\
    .agg(shannon=('shannon', 'prod'), richness=('bstrain_id', 'size')).reset_index()\
                        .merge(shannon,on=['run_id','t'])\
                        .merge(richness, on=['run_id', 't'])
cladeImmune['shanP'] = cladeImmune['shannon']/cladeImmune['stotal']
cladeImmune['richP'] = cladeImmune['richness']/cladeImmune['rtotal']

for runID in np.unique(cladeImmune['run_id']):
    # print('runID: {}'.format(runID))
    for gRate in np.unique(cladeImmune['growth_rate']):
        maxT = int(max(cladeImmune[(cladeImmune.run_id == runID)
                                      & (cladeImmune.growth_rate == gRate)]['t']))
        times = list(range(0, maxT+1, 1))
        mTimes = list(set(times) - \
                      set(cladeImmune[(cladeImmune.run_id == runID)
                                      & (cladeImmune.growth_rate == gRate)]['t']))
        mTimes.sort()
        timeInt = []
        for i in range(0, len(mTimes), 1):
            # print('time_idx: {}'.format(i))
            if i == 0:
                timeInt.append(mTimes[i])
            elif i == len(mTimes)-1:
                if mTimes[i] - mTimes[i-1] == 1:
                    timeInt.append(mTimes[i])
                    t = min(timeInt)
                    missing = cladeImmune[(cladeImmune.run_id == runID)
                                & (cladeImmune.growth_rate == gRate)
                                & (cladeImmune.t <= t)].copy()
                    missing = cladeImmune[(cladeImmune.run_id == runID)
                                & (cladeImmune.growth_rate == gRate) 
                                & (cladeImmune.t == max(missing['t']))].copy()
                    missing = pd.DataFrame(
                        {'run_id': len(timeInt)*runID,
                         't': timeInt,
                         'growth_rate': len(timeInt)*[gRate],
                         'shannon': len(timeInt)*missing['shannon'].values,
                         'stotal' : len(timeInt)*missing['stotal'].values,
                         'shanP' : len(timeInt)*missing['shanP'].values,
                         'richness': len(timeInt)*missing['richness'].values,
                         'rtotal': len(timeInt)*missing['rtotal'].values,
                         'richP': len(timeInt)*missing['richP'].values
                         })
                    cladeImmune = pd.concat([cladeImmune, missing], ignore_index=True)
                    timeInt = [mTimes[i]]
                else:
                    t = min(mTimes[i])
                    missing = cladeImmune[(cladeImmune.run_id == runID)
                                & (cladeImmune.growth_rate == gRate)
                                & (cladeImmune.t <= t)].copy()
                    missing = cladeImmune[(cladeImmune.run_id == runID)
                                & (cladeImmune.growth_rate == gRate) 
                                & (cladeImmune.t == max(missing['t']))].copy()
                    missing = pd.DataFrame(
                        {'run_id': len(timeInt)*runID,
                         't': timeInt,
                         'growth_rate': len(timeInt)*[gRate],
                         'shannon': len(timeInt)*missing['shannon'].values,
                         'stotal' : len(timeInt)*missing['stotal'].values,
                         'shanP': len(timeInt)*missing['shanP'].values,
                         'richness': len(timeInt)*missing['richness'].values,
                         'rtotal': len(timeInt)*missing['rtotal'].values,
                         'richP': len(timeInt)*missing['richP'].values
                         })
                    cladeImmune = pd.concat([cladeImmune, missing], ignore_index=True)
                    timeInt = [mTimes[i]]
            else:
                if mTimes[i] - mTimes[i-1] == 1:
                    timeInt.append(mTimes[i])
                else:
                    t = min(timeInt)
                    missing = cladeImmune[(cladeImmune.run_id == runID)
                                & (cladeImmune.growth_rate == gRate)
                                & (cladeImmune.t <= t)].copy()
                    missing = cladeImmune[(cladeImmune.run_id == runID)
                                & (cladeImmune.growth_rate == gRate) 
                                & (cladeImmune.t == max(missing['t']))].copy()
                    missing = pd.DataFrame(
                        {'run_id': len(timeInt)*runID,
                         't': timeInt,
                         'growth_rate': len(timeInt)*[gRate],
                         'shannon': len(timeInt)*missing['shannon'].values,
                         'stotal' : len(timeInt)*missing['stotal'].values,
                         'shanP': len(timeInt)*missing['shanP'].values,
                         'richness': len(timeInt)*missing['richness'].values,
                         'rtotal': len(timeInt)*missing['rtotal'].values,
                         'richP': len(timeInt)*missing['richP'].values
                         })
                    cladeImmune = pd.concat([cladeImmune, missing], ignore_index=True)
                    timeInt = [mTimes[i]]

cladeImmune = cladeImmune.sort_values(by=['run_id', 't','growth_rate'])
cladeIstats = cladeImmune.groupby(['t','growth_rate'])\
                .agg(shanPmean=('shanP','mean'),n=('shanP','size'),
                     shanPstd=('shanP','std'),
                     shanPmedian=('shanP','median'),
                     richPmean=('richP', 'mean'),
                     richPstd=('richP', 'std'),
                     richPmedian=('richP', 'median'))\
                .reset_index()
cladeIstats['growth_rate'] = cladeIstats['growth_rate'] - .025
cladeIstatsTrunc = cladeIstats[(cladeIstats.n >= 150) & (cladeIstats.t != 0)]
# cladeIstatsTrunc['t'] = list(np.log(cladeIstatsTrunc['t']))
cladeIstatsTrunc['growth_rate'] = np.round(cladeIstatsTrunc['growth_rate'],3)
###
###
###
virus_total = pd.read_sql_query("SELECT run_id, t, viral_abundance FROM summary \
                                WHERE run_id in ({}) AND viral_abundance > 0"
                                .format(', '.join(map(str, runIDs['run_id']))), conSim)\
    .rename(columns={"viral_abundance": "vtotal"})
vstats = virus_total.groupby(['t'])\
    .agg(exp_vtotal=('vtotal', 'mean'), std_vtotal=('vtotal', 'std'), n=('vtotal', 'size'))\
    .reset_index()
vstats = vstats[vstats.n >= 50]
#
fig = plt.figure()
gs = plt.GridSpec(4, 2, figure=fig, hspace=.065, wspace=.6)
ax1 = fig.add_subplot(gs[0:2, 0])
ax2 = fig.add_subplot(gs[0:1, 1])
ax3 = fig.add_subplot(gs[1:2, 1], sharex=ax2)
ax4 = fig.add_subplot(gs[2:4, 0], sharex=ax1)
ax5 = fig.add_subplot(gs[2:4, 1],sharex=ax2)
axes = [ax1, ax2, ax3, ax4, ax5]
axes = [ax1, ax2, ax3, ax4, ax5,
        ax1.twinx(),ax2.twinx(),ax3.twinx(),
        ax4.twinx()]
axes[0].set_xscale('log', base=10)
# axes[1].set_xscale('log', base=10)
axes[2].set_xscale('log', base=10)
for i in range(0,4):
    # i = 0
    axes[i].fill_between(vstats['t'],
                        vstats['exp_vtotal'] -
                        vstats['std_vtotal'],
                        vstats['exp_vtotal'] +
                        vstats['std_vtotal'], color='grey', alpha=0.1)
    axes[i].plot(vstats['t'],
                vstats['exp_vtotal'],
                linewidth=2, color='grey', label='Viral Abund.', linestyle='solid', alpha=0.75)
    lim = axes[i].get_ylim()
    axes[i].set_ylim(0, lim[1])
    axes[i].yaxis.tick_right()
    # axes[i].set_xlim(0, 2000)
    axes[i].tick_params(axis='x', labelsize=15)
    # axes[i].set_xscale('log', base=10)
    if i == 0:
        axes[i].set_ylabel(ylabel='Viral Abundance',
                           rotation=270, labelpad=25, fontsize=15)
        axes[i].tick_params(axis='y', labelsize=15)
        axes[i].yaxis.get_offset_text().set_fontsize(15)
        axes[i].yaxis.set_label_position("right")
    if i != 0:
        print(i)
        axes[i].set_yticks([])
#
axes[0].tick_params(axis='x', labelcolor='w', top=False,
                    bottom=False, left=False, right=False)
axes[1].tick_params(axis='x', labelcolor='w', top=False,
                    bottom=False, left=False, right=False)
axes[2].tick_params(axis='x', labelcolor='w', top=False,
                    bottom=False, left=False, right=False)
axes[5].fill_between(meanDist[meanDist['combo_id'] == cID]['t'],
                     meanDist[meanDist['combo_id'] == cID]['meanExp'] -
                     meanDist[meanDist['combo_id'] == cID]['stdExp'],
                     meanDist[meanDist['combo_id'] == cID]['meanExp'] +
                     meanDist[meanDist['combo_id'] == cID]['stdExp'], color='mediumblue', alpha=0.3)
axes[5].plot(meanDist[meanDist['combo_id'] == cID]['t'],
             meanDist[meanDist['combo_id'] == cID]['meanExp'],
             linewidth=2, color='mediumblue', label=r'$\sigma = 3$')
axes[5].fill_between(meanNoVDist[meanNoVDist['combo_id'] == cIDNoV]['t'],
                     meanNoVDist[meanNoVDist['combo_id'] == cIDNoV]['meanExp'] -
                     meanNoVDist[meanNoVDist['combo_id']
                                 == cIDNoV]['stdExp'],
                     meanNoVDist[meanNoVDist['combo_id'] == cIDNoV]['meanExp'] +
                     meanNoVDist[meanNoVDist['combo_id'] == cIDNoV]['stdExp'], color='lime', alpha=0.3)
axes[5].plot(meanNoVDist[meanNoVDist['combo_id'] == cIDNoV]['t'],
             meanNoVDist[meanNoVDist['combo_id'] == cIDNoV]['meanExp'],
             linewidth=2, color='lime', label=r'$\sigma = 3$, no virus', linestyle='dashed')
axes[5].tick_params(axis='x', labelsize=15)
axes[5].tick_params(axis='y', labelsize=15)
axes[5].set_ylabel(
    ylabel=r'Expected Fitness Mean $\mathbb{E}[\bar{f}]$', labelpad=10, fontsize=15)
axes[5].yaxis.set_label_position("left")
axes[5].yaxis.tick_left()
lim = axes[5].get_ylim()
axes[5].set_ylim(0.4, lim[1])
axes[8].fill_between(meanDist[meanDist['combo_id'] == cID]['t'],
                     meanDist[meanDist['combo_id'] == cID]['meanStd'] -
                     meanDist[meanDist['combo_id'] == cID]['stdStd'],
                     meanDist[meanDist['combo_id'] == cID]['meanStd'] +
                     meanDist[meanDist['combo_id'] == cID]['stdStd'], color='mediumblue', alpha=0.3)
axes[8].plot(meanDist[meanDist['combo_id'] == cID]['t'],
             meanDist[meanDist['combo_id'] == cID]['meanStd'],
             linewidth=2, color='mediumblue', label=r'$\sigma=3$')
axes[8].fill_between(meanNoVDist[meanNoVDist['combo_id'] == cIDNoV]['t'],
                     meanNoVDist[meanNoVDist['combo_id'] == cIDNoV]['meanStd'] -
                     meanNoVDist[meanNoVDist['combo_id']
                                 == cIDNoV]['stdStd'],
                     meanNoVDist[meanNoVDist['combo_id'] == cIDNoV]['meanStd'] +
                     meanNoVDist[meanNoVDist['combo_id'] == cIDNoV]['stdStd'], color='lime', alpha=0.3)
axes[8].plot(meanNoVDist[meanNoVDist['combo_id'] == cIDNoV]['t'],
             meanNoVDist[meanNoVDist['combo_id'] == cIDNoV]['meanStd'],
             linewidth=2, color='lime', label=r'$\sigma = 3$, no virus', linestyle='dashed')
axes[3].set_xlabel(xlabel='Time t', fontsize=15, labelpad=15)
axes[6].tick_params(axis='x', labelsize=15)
axes[6].tick_params(axis='y', labelsize=15)
axes[7].tick_params(axis='x', labelsize=15)
axes[7].tick_params(axis='y', labelsize=15)
axes[8].tick_params(axis='x', labelsize=15)
axes[8].tick_params(axis='y', labelsize=15)
axes[4].tick_params(axis='x', labelsize=15)
axes[4].tick_params(axis='y', labelsize=15)
axes[4].set_xlabel(xlabel='Time t', fontsize=15, labelpad=15)
s = ['Expected Fitness SD ', r'$\mathbb{E}[\sigma_f]$']
axes[8].set_ylabel(ylabel=r''.join(s), labelpad=10, fontsize=15)
axes[8].set_xlabel(xlabel = 'Time t',fontsize=15,labelpad=15)
# axes[8].legend(loc='upper right',fontsize=15)
gRates = sorted(np.unique(cladeIstatsTrunc['growth_rate']))
gRateMap = cm.get_cmap('hot').reversed()
norm = Normalize(
    vmin=-1*float(len(gRates))/8, vmax=float(len(gRates)))
for i in range(len(gRates)-1, 0-1, -1):
    gRate = gRates[i]
    axes[4].fill_between(cladeIstatsTrunc[cladeIstatsTrunc.growth_rate == gRate]['t'],
                    cladeIstatsTrunc[cladeIstatsTrunc.growth_rate == gRate]['richPmean'] -
                    cladeIstatsTrunc[cladeIstatsTrunc.growth_rate ==
                                     gRate]['richPstd'],
                    cladeIstatsTrunc[cladeIstatsTrunc.growth_rate == gRate]['richPmean'] +
                    cladeIstatsTrunc[cladeIstatsTrunc.growth_rate == gRate]['richPstd'], 
                    color=gRateMap(norm(float(i))), alpha=0.3)
for i in range(len(gRates)-1, 0-1, -1):
    print(i)
    gRate = gRates[i]
    axes[4].plot(cladeIstatsTrunc[cladeIstatsTrunc.growth_rate == gRate]['t'],
            cladeIstatsTrunc[cladeIstatsTrunc.growth_rate ==
                             gRate]['richPmean'],
            color=gRateMap(norm(float(i))), label='{}'.format(gRate), linewidth=2, linestyle='solid',alpha=0.75)
axes[7].fill_between(richnessStats['t'],
                     richnessStats['richMean'] -
                     richnessStats['richStd'],
                     richnessStats['richMean'] +
                     richnessStats['richStd'], color='mediumblue', alpha=0.3)
axes[7].plot(richnessStats['t'],
             richnessStats['richMean'],
             linewidth=2, color='mediumblue')
axes[6].fill_between(simpsonStats['t'],
                     simpsonStats['simpMean'] -
                     simpsonStats['simpStd'],
                     simpsonStats['simpMean'] +
                     simpsonStats['simpStd'], color='mediumblue', alpha=0.3)
axes[6].plot(simpsonStats['t'],
             simpsonStats['simpMean'],
             linewidth=2, color='mediumblue')
axes[6].fill_between(divNoV['t'],
                     divNoV['simpMean'] -
                     divNoV['simpStd'],
                     divNoV['simpMean'] +
                     divNoV['simpStd'], color='lime', alpha=0.3)
axes[6].plot(divNoV['t'],
             divNoV['simpMean'],
             linewidth=2, color='lime', linestyle='dashed')
axes[7].fill_between(divNoV['t'],
                     divNoV['richMean'] -
                     divNoV['richStd'],
                     divNoV['richMean'] +
                     divNoV['richStd'], color='lime', alpha=0.3)
axes[7].plot(divNoV['t'],
             divNoV['richMean'],
             linewidth=2, color='lime', linestyle='dashed')
axes[6].set_ylabel(ylabel ='Simpson Diversity',labelpad=10,fontsize=15)
axes[7].set_ylabel(ylabel='Richness', labelpad=10, fontsize=15)
axes[4].set_ylabel(ylabel='Proportion of Richness', labelpad=10, fontsize=15)
lim = axes[5].get_xlim()
axes[5].set_xlim(lim[0], 2000)
lim = axes[6].get_xlim()
axes[6].set_xlim(lim[0], 2000)
lim = axes[7].get_xlim()
axes[7].set_xlim(lim[0], 2000)
lim = axes[8].get_xlim()
axes[8].set_xlim(lim[0], 2000)
lim = axes[4].get_xlim()
axes[4].set_xlim(lim[0], 2000)
axes[6].yaxis.set_label_position("left")
axes[6].yaxis.tick_left()
axes[7].yaxis.set_label_position("left")
axes[7].yaxis.tick_left()
axes[8].yaxis.set_label_position("left")
axes[8].yaxis.tick_left()
axes[4].yaxis.set_label_position("left")
axes[4].yaxis.tick_left()
handles = []
labels = []
handle, label = axes[0].get_legend_handles_labels()
handles.extend(handle)
labels.extend(label)
handle, label = axes[5].get_legend_handles_labels()
handles.extend(handle)
labels.extend(label)
axes[5].legend(handles, labels, loc='lower left', fontsize=12)
axes[4].legend(loc='upper left', fontsize=10, title='Host Fitness', title_fontsize=10)
plt.show()


virus_total = pd.read_sql_query("SELECT run_id, t, viral_abundance FROM summary \
                                WHERE run_id in ({}) AND viral_abundance > 0"
                                  .format(', '.join(map(str, runIDs['run_id']))), conSim)\
    .rename(columns={"viral_abundance": "vtotal"})
vstats = virus_total.groupby(['t'])\
    .agg(exp_vtotal=('vtotal', 'mean'), std_vtotal=('vtotal', 'std'), n = ('vtotal', 'size'))\
        .reset_index()
vstats = vstats[vstats.n>=50]


fig = plt.figure()
gs = plt.GridSpec(2, 1, figure=fig, hspace=.075)
ax1 = fig.add_subplot(gs[0:1, 0])
ax2 = fig.add_subplot(gs[1:2, 0], sharex=ax1)
axes = [ax1, ax2, ax1.twinx(), ax2.twinx()]
axes[2].fill_between(vstats['t'],
                     vstats['exp_vtotal'] -
                     vstats['std_vtotal'],
                     vstats['exp_vtotal'] +
                     vstats['std_vtotal'], color='grey', alpha=0.1)
axes[2].plot(vstats['t'],
             vstats['exp_vtotal'],
             linewidth=2, color='grey', label='Viral Abund.', linestyle='solid', alpha = 0.75)
axes[3].fill_between(vstats['t'],
                     vstats['exp_vtotal'] -
                     vstats['std_vtotal'],
                     vstats['exp_vtotal'] +
                     vstats['std_vtotal'], color='grey', alpha=0.1)
axes[3].plot(vstats['t'],
             vstats['exp_vtotal'],
             linewidth=2, color='grey', label='Viral Abund.', linestyle='solid', alpha=0.75)
axes[0].fill_between(meanNoVDist[meanNoVDist['combo_id'] == cIDNoV]['t'],
                     meanNoVDist[meanNoVDist['combo_id'] == cIDNoV]['meanExp'] -
                     meanNoVDist[meanNoVDist['combo_id']
                                 == cIDNoV]['stdExp'],
                     meanNoVDist[meanNoVDist['combo_id'] == cIDNoV]['meanExp'] +
                     meanNoVDist[meanNoVDist['combo_id'] == cIDNoV]['stdExp'], color='lime', alpha=0.3)
axes[0].fill_between(meanDist[meanDist['combo_id'] == cID]['t'],
                     meanDist[meanDist['combo_id'] == cID]['meanExp'] -
                     meanDist[meanDist['combo_id'] == cID]['stdExp'],
                     meanDist[meanDist['combo_id'] == cID]['meanExp'] +
                     meanDist[meanDist['combo_id'] == cID]['stdExp'], color='mediumblue', alpha=0.3)
axes[0].plot(meanDist[meanDist['combo_id'] == cID]['t'],
             meanDist[meanDist['combo_id'] == cID]['meanExp'],
             linewidth=2, color='mediumblue', label=r'$\sigma = 3$')
axes[0].plot(meanNoVDist[meanNoVDist['combo_id'] == cIDNoV]['t'],
             meanNoVDist[meanNoVDist['combo_id'] == cIDNoV]['meanExp'],
             linewidth=2, color='lime', label=r'$\sigma = 3$, No Virus', linestyle='dashed')
axes[0].set_xscale('log', base=10)
axes[1].set_xscale('log', base=10)
axes[0].tick_params(axis='x', labelsize=15)
axes[0].tick_params(axis='y', labelsize=15)
axes[1].tick_params(axis='x', labelsize=15)
axes[1].tick_params(axis='y', labelsize=15)
axes[0].tick_params(axis='x', labelcolor='w', top=False,
                    bottom=False, left=False, right=False)
# axes[1].tick_params(axis='x', labelcolor='w', top=False,
#                     bottom=False, left=False, right=False)
axes[0].set_ylabel(
    ylabel=r'Expected Fitness Mean $\mathbb{E}[\bar{f}]$', labelpad=10, fontsize=15)
axes[1].fill_between(simpsonStats['t'],
                     simpsonStats['simpMean'] -
                     simpsonStats['simpStd'],
                     simpsonStats['simpMean'] +
                     simpsonStats['simpStd'], color='mediumblue', alpha=0.3)
axes[1].plot(simpsonStats['t'],
             simpsonStats['simpMean'],
             linewidth=2, color='mediumblue')
axes[1].fill_between(divNoV['t'],
                     divNoV['simpMean'] -
                     divNoV['simpStd'],
                     divNoV['simpMean'] +
                     divNoV['simpStd'], color='lime', alpha=0.3)
axes[1].plot(divNoV['t'],
             divNoV['simpMean'],
             linewidth=2, color='lime', linestyle='dashed')
axes[1].set_ylabel(ylabel='Simpson\nImmune Strain Diversity', labelpad=10, fontsize=15)
axes[1].set_xlabel(ylabel='Time', labelpad=10, fontsize=15)
axes[0].set_xlim(0, 2000)
axes[1].set_xlim(0, 2000)
axes[2].set_xlim(0, 2000)
axes[3].set_xlim(0, 2000)
lim = axes[2].get_ylim()
axes[2].set_ylim(0,lim[1])
axes[2].set_ylabel(ylabel='Viral Abundance',
                   rotation=270, labelpad=25, fontsize=15)
lim = axes[3].get_ylim()
axes[3].set_ylim(0, lim[1])
axes[3].set_ylabel(ylabel='Viral Abundance',
                   rotation=270, labelpad=25, fontsize=15)
# axes[2].set_yticks([])
# axes[2].set_yticklabels([])
# axes[0].legend(loc='lower left', fontsize=15)
# axes[3].legend(loc='upper left', fontsize=15)
axes[2].tick_params(axis='x', labelsize=15)
axes[2].tick_params(axis='y', labelsize=15)
axes[3].tick_params(axis='x', labelsize=15)
axes[3].tick_params(axis='y', labelsize=15)
handles = []
labels = []
handle, label = axes[0].get_legend_handles_labels()
handles.extend(handle)
labels.extend(label)
handle, label = axes[1].get_legend_handles_labels()
handles.extend(handle)
labels.extend(label)
handle, label = axes[2].get_legend_handles_labels()
handles.extend(handle)
labels.extend(label)
axes[0].legend(handles, labels, loc='upper left', fontsize=15)
axes[2].yaxis.get_offset_text().set_fontsize(15)
axes[3].yaxis.get_offset_text().set_fontsize(15)
plt.show()


fig = plt.figure()
gs = plt.GridSpec(2, 1, figure=fig, hspace=.075)
ax1 = fig.add_subplot(gs[0:1, 0])
ax2 = fig.add_subplot(gs[1:2, 0], sharex=ax1)
axes = [ax1, ax2, ax1.twinx()]
axes[0].set_xscale('log', base=10)
axes[1].set_xscale('log', base=10)
axes[2].fill_between(vstats['t'],
                     vstats['exp_vtotal'] -
                     vstats['std_vtotal'],
                     vstats['exp_vtotal'] +
                     vstats['std_vtotal'], color='grey', alpha=0.1)
axes[2].plot(vstats['t'],
             vstats['exp_vtotal'],
             linewidth=2, color='grey', label='Viral Abund.', linestyle='solid', alpha=0.75)
# axes[3].fill_between(vstats['t'],
#                      vstats['exp_vtotal'] -
#                      vstats['std_vtotal'],
#                      vstats['exp_vtotal'] +
#                      vstats['std_vtotal'], color='grey', alpha=0.1)
# axes[3].plot(vstats['t'],
#              vstats['exp_vtotal'],
#              linewidth=2, color='grey', label='Viral Abund.', linestyle='solid', alpha=0.75)
axes[0].tick_params(axis='x', labelsize=15)
axes[0].tick_params(axis='y', labelsize=15)
axes[1].tick_params(axis='x', labelsize=15)
axes[1].tick_params(axis='y', labelsize=15)
axes[2].tick_params(axis='x', labelsize=15)
axes[2].tick_params(axis='y', labelsize=15)
axes[0].tick_params(axis='x', labelcolor='w', top=False,
                    bottom=False, left=False, right=False)
# axes[1].set_xlabel(xlabel='Time', fontsize=15, labelpad=15)
axes[1].set_xlabel(ylabel='Time', labelpad=10, fontsize=15)
gRates = sorted(np.unique(cladeIstatsTrunc['growth_rate']))
gRateMap = cm.get_cmap('inferno').reversed()
norm = Normalize(
    vmin=float(0), vmax=float(len(gRates)))
for i in range(len(gRates)-1, 0-1, -1):
    gRate = gRates[i]
    axes[1].fill_between(cladeIstatsTrunc[cladeIstatsTrunc.growth_rate == gRate]['t'],
                         cladeIstatsTrunc[cladeIstatsTrunc.growth_rate == gRate]['richPmean'] -
                         cladeIstatsTrunc[cladeIstatsTrunc.growth_rate ==
                                          gRate]['richPstd'],
                         cladeIstatsTrunc[cladeIstatsTrunc.growth_rate == gRate]['richPmean'] +
                         cladeIstatsTrunc[cladeIstatsTrunc.growth_rate == gRate]['richPstd'], color=gRateMap(norm(float(i))), alpha=0.3)
for i in range(len(gRates)-1, 0-1, -1):
    print(i)
    gRate = gRates[i]
    axes[1].plot(cladeIstatsTrunc[cladeIstatsTrunc.growth_rate == gRate]['t'],
                 cladeIstatsTrunc[cladeIstatsTrunc.growth_rate ==
                                  gRate]['richPmean'],
                 color=gRateMap(norm(float(i))), label='{}'.format(gRate), linewidth=2, linestyle='solid', alpha=0.75)
axes[0].fill_between(divNoV['t'],
                     divNoV['richMean'] -
                     divNoV['richStd'],
                     divNoV['richMean'] +
                     divNoV['richStd'], color='lime', alpha=0.3)
axes[0].fill_between(richnessStats['t'],
                     richnessStats['richMean'] -
                     richnessStats['richStd'],
                     richnessStats['richMean'] +
                     richnessStats['richStd'], color='mediumblue', alpha=0.3)
axes[0].plot(richnessStats['t'],
             richnessStats['richMean'],
             linewidth=2, color='mediumblue', label=r'$\sigma = 3$')
axes[0].plot(divNoV['t'],
             divNoV['richMean'],
             linewidth=2, color='lime', linestyle='dashed', label=r'$\sigma = 3$, No Virus')
axes[0].set_ylabel(ylabel='Richness', labelpad=10, fontsize=15)
axes[1].set_ylabel(ylabel='Proportion of Richness', labelpad=10, fontsize=15)
axes[0].set_xlim(0, 2000)
axes[1].set_xlim(0, 2000)
axes[2].set_xlim(0, 2000)
# axes[3].set_xlim(0, 2000)
lim = axes[2].get_ylim()
axes[2].set_ylim(0, lim[1])
axes[2].set_ylabel(ylabel='Viral Abundance',
                   rotation=270, labelpad=25, fontsize=15)
lim = axes[3].get_ylim()
# axes[3].set_ylim(0, lim[1])
# axes[3].set_ylabel(ylabel='Viral Abundance',
#                    rotation=270, labelpad=25, fontsize=15)
# axes[2].set_yticks([])
# axes[2].set_yticklabels([])
# axes[0].legend(loc='lower left', fontsize=15)
axes[1].legend(loc='upper left', fontsize=10,
               title='Fitness', title_fontsize=10)
handles = []
labels = []
handle, label = axes[0].get_legend_handles_labels()
handles.extend(handle)
labels.extend(label)
handle, label = axes[2].get_legend_handles_labels()
handles.extend(handle)
labels.extend(label)
axes[0].legend(handles, labels, loc='upper left', fontsize=15)
axes[2].yaxis.get_offset_text().set_fontsize(15)
# axes[3].yaxis.get_offset_text().set_fontsize(15)

plt.show()

def fgm(s,x):
    return np.exp(-s*x**2)

domain = list(np.linspace(-2, 2, 100))
selection = [0, .5, 1, 1.5, 2, 2.5, 3]
alleles = [0.0, -0.5, 0.54, 0.12, 0.72, -0.92, 0.9, -0.64]
ymax = 1.25
xmax = 1.5
for s in selection:
    fig, ax = plt.subplots(1,figsize=(8,8/1.5))
    ax.plot(domain,list(map(lambda x: fgm(s,x), domain)),label=r'$\sigma =$ {0}'.format(s),\
            color='black',linewidth=2)
    ax.set_ylim(0,ymax)
    ax.set_xlim(-xmax, xmax)
    ax.legend()
    ax.set_xlabel(xlabel='Trait Value',labelpad=10,fontsize=15)
    ax.set_ylabel(ylabel='Fitness', labelpad=10, fontsize=15)
    ax.tick_params(axis='x', labelsize=15)
    ax.tick_params(axis='y', labelsize=15)
    for a in alleles:
        ax.axvline(x=a, ymax=fgm(s,a)/ymax, linestyle='dotted', color='darkred', linewidth=1)
    fig.tight_layout()
    fig.savefig(os.path.join('/Users/armun/Desktop/fgmSigma{0}.png'.format(s)),dpi=500)



