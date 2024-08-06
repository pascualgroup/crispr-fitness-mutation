import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import sqlite3
import matplotlib.cm as cm
from matplotlib.colors import Normalize
import matplotlib.colors as mc
from matplotlib.collections import LineCollection
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.ticker as ticker
import seaborn as sns

DBSIM_PATH = os.path.join('/Volumes/Yadgah/sylvain-martin-collab/12_MOI3/sweep_db_gathered.sqlite')
DBTRAIT_PATH = os.path.join('/Volumes/Yadgah/sylvain-martin-collab/12_MOI3/gathered-analyses/trait-evolution/trait-evolution.sqlite')
DBSIM_PATHnull = os.path.join('/Volumes/Yadgah/sylvain-martin-collab/14_MOI3/sweep_db_gathered.sqlite')
DBSIM_PATHnoV = os.path.join('/Volumes/Yadgah/sylvain-martin-collab/13_MOI3/sweep_db_gathered.sqlite')
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
cID2 = 3
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

cladeImmune2 = pd.read_sql_query("SELECT bstrain_id, growth_rate, run_id FROM bgrowthrates \
                                    WHERE run_id in ({0})"
                             .format(', '.join(map(str, runIDs[runIDs.combo_id == cID2]['run_id']))), conSim)\
                             .merge(\
                               pd.read_sql_query("SELECT t, bstrain_id, abundance, run_id FROM babundance \
                                                    WHERE run_id in ({0})"
                                                .format(', '.join(\
                                                    map(str, runIDs[runIDs.combo_id == cID2]['run_id'])\
                                                        )), conSim), on = ['run_id','bstrain_id'] 
                             )\
                             .merge(\
                               pd.read_sql_query("SELECT t, microbial_abundance, run_id FROM summary \
                                                    WHERE run_id in ({0})"
                                                .format(', '.join(\
                                                    map(str, runIDs[runIDs.combo_id == cID2]['run_id'])\
                                                        )), conSim), on = ['run_id','t'] 
                             )\
                             .rename(columns={'microbial_abundance':'btotal'})                            
cladeImmune2['bfreq'] = cladeImmune2['abundance']/cladeImmune2['btotal']
cladeImmune2 = cladeImmune2.drop(columns=['abundance','btotal'])
cladeImmune2['shannon'] = np.exp(-1*np.array(cladeImmune2['bfreq'])*np.log(cladeImmune2['bfreq']))
cladeImmune2['simpson'] = np.array(cladeImmune2['bfreq'])**2
shannon2 = cladeImmune2.groupby(['t','run_id']).agg(stotal=('shannon','prod')).reset_index()
simpson2 = cladeImmune2.groupby(['t','run_id']).agg(ltotal=('simpson','sum')).reset_index()
simpson2['ltotal'] = 1/simpson2['ltotal']
richness2 = cladeImmune2[['run_id','t','bstrain_id']].drop_duplicates().groupby(['t', 'run_id']).agg(
    rtotal=('bstrain_id', 'size')).reset_index()
shannonStats2 = shannon2.groupby(['t']).agg(
    shanMean=('stotal', 'mean'), shanStd=('stotal', 'std')).reset_index()
simpsonStats2 = simpson2.groupby(['t']).agg(
    simpMean=('ltotal', 'mean'), simpStd=('ltotal', 'std')).reset_index()
richnessStats2 = richness2.groupby(['t']).agg(
    richMean=('rtotal', 'mean'), richStd=('rtotal', 'std')).reset_index()

cladeImmune2 = cladeImmune2.groupby(['run_id','t','growth_rate'])\
    .agg(shannon=('shannon', 'prod'), richness=('bstrain_id', 'size')).reset_index()\
                        .merge(shannon2,on=['run_id','t'])\
                        .merge(richness2, on=['run_id', 't'])
cladeImmune2['shanP'] = cladeImmune2['shannon']/cladeImmune2['stotal']
cladeImmune2['richP'] = cladeImmune2['richness']/cladeImmune2['rtotal']
for runID in np.unique(cladeImmune2['run_id']):
    # print('runID: {}'.format(runID))
    for gRate in np.unique(cladeImmune2['growth_rate']):
        maxT = int(max(cladeImmune2[(cladeImmune2.run_id == runID)
                                      & (cladeImmune2.growth_rate == gRate)]['t']))
        times = list(range(0, maxT+1, 1))
        mTimes = list(set(times) - \
                      set(cladeImmune2[(cladeImmune2.run_id == runID)
                                      & (cladeImmune2.growth_rate == gRate)]['t']))
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
                    missing = cladeImmune2[(cladeImmune2.run_id == runID)
                                & (cladeImmune2.growth_rate == gRate)
                                & (cladeImmune2.t <= t)].copy()
                    missing = cladeImmune2[(cladeImmune2.run_id == runID)
                                & (cladeImmune2.growth_rate == gRate) 
                                & (cladeImmune2.t == max(missing['t']))].copy()
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
                    cladeImmune2 = pd.concat([cladeImmune2, missing], ignore_index=True)
                    timeInt = [mTimes[i]]
                else:
                    t = min(mTimes[i])
                    missing = cladeImmune2[(cladeImmune2.run_id == runID)
                                & (cladeImmune2.growth_rate == gRate)
                                & (cladeImmune2.t <= t)].copy()
                    missing = cladeImmune2[(cladeImmune2.run_id == runID)
                                & (cladeImmune2.growth_rate == gRate) 
                                & (cladeImmune2.t == max(missing['t']))].copy()
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
                    cladeImmune2 = pd.concat([cladeImmune2, missing], ignore_index=True)
                    timeInt = [mTimes[i]]
            else:
                if mTimes[i] - mTimes[i-1] == 1:
                    timeInt.append(mTimes[i])
                else:
                    t = min(timeInt)
                    missing = cladeImmune2[(cladeImmune2.run_id == runID)
                                & (cladeImmune2.growth_rate == gRate)
                                & (cladeImmune2.t <= t)].copy()
                    missing = cladeImmune2[(cladeImmune2.run_id == runID)
                                & (cladeImmune2.growth_rate == gRate) 
                                & (cladeImmune2.t == max(missing['t']))].copy()
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
                    cladeImmune2 = pd.concat([cladeImmune2, missing], ignore_index=True)
                    timeInt = [mTimes[i]]

cladeImmune2 = cladeImmune2.sort_values(by=['run_id', 't','growth_rate'])
cladeIstats2 = cladeImmune2.groupby(['t','growth_rate'])\
                .agg(shanPmean=('shanP','mean'),n=('shanP','size'),
                     shanPstd=('shanP','std'),
                     shanPmedian=('shanP','median'),
                     richPmean=('richP', 'mean'),
                     richPstd=('richP', 'std'),
                     richPmedian=('richP', 'median'))\
                .reset_index()
cladeIstats2['growth_rate'] = cladeIstats2['growth_rate'] - .025
cladeIstatsTrunc2 = cladeIstats2[(cladeIstats2.n >= 150) & (cladeIstats2.t != 0)]
# cladeIstatsTrunc['t'] = list(np.log(cladeIstatsTrunc['t']))
cladeIstatsTrunc2['growth_rate'] = np.round(cladeIstatsTrunc2['growth_rate'],3)





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
virus_total = pd.read_sql_query("SELECT run_id, t, viral_abundance FROM summary \
                                WHERE run_id in ({}) AND viral_abundance > 0"
                                .format(', '.join(map(str, runIDs[runIDs.combo_id==cID]['run_id']))), conSim)\
    .rename(columns={"viral_abundance": "vtotal"})
vstats = virus_total.groupby(['t'])\
    .agg(exp_vtotal=('vtotal', 'mean'), std_vtotal=('vtotal', 'std'), n=('vtotal', 'size'))\
    .reset_index()
vstats = vstats[vstats.n >= 2]



fig = plt.figure(figsize=(20, 8))
gs = fig.add_gridspec(1, 2, hspace=0, wspace=.35)
(ax1, ax2) = gs.subplots()
axes = [ax1, ax2]
gRates = sorted(np.unique(cladeIstatsTrunc['growth_rate']))
gRateMap = cm.get_cmap('inferno').reversed()
norm = Normalize(
    vmin=float(0), vmax=float(len(gRates)))
# speciesIDX = strainTimes[strainTimes[tree_parent_strain_id] == 0].sort_values(by='growth_rate')[tree_strain_id]
for i in range(len(gRates)-1, 0-1, -1):
    print(i)
    gRate = gRates[i]
    axes[0].fill_between(cladeIstatsTrunc[cladeIstatsTrunc.growth_rate == gRate]['t'],
                    cladeIstatsTrunc[cladeIstatsTrunc.growth_rate == gRate]['richPmean'] -
                    cladeIstatsTrunc[cladeIstatsTrunc.growth_rate ==
                                    gRate]['richPstd'],
                    cladeIstatsTrunc[cladeIstatsTrunc.growth_rate == gRate]['richPmean'] +
                    cladeIstatsTrunc[cladeIstatsTrunc.growth_rate == gRate]['richPstd'], 
                    color=gRateMap(norm(float(i))), alpha=0.3)
for i in range(len(gRates)-1, 0-1, -1):
    print(i)
    gRate = gRates[i]
    axes[0].plot(cladeIstatsTrunc[cladeIstatsTrunc.growth_rate == gRate]['t'],
            cladeIstatsTrunc[cladeIstatsTrunc.growth_rate ==
                            gRate]['richPmean'],
            color=gRateMap(norm(float(i))), label='{}'.format(gRate), linewidth=2, linestyle='solid',alpha=0.75)
    
gRates = sorted(np.unique(cladeIstatsTrunc2['growth_rate']))
norm = Normalize(
    vmin=float(0), vmax=float(len(gRates)))
# speciesIDX = strainTimes[strainTimes[tree_parent_strain_id] == 0].sort_values(by='growth_rate')[tree_strain_id]
for i in range(len(gRates)-1, 0-1, -1):
    print(i)
    gRate = gRates[i]
    axes[1].fill_between(cladeIstatsTrunc2[cladeIstatsTrunc2.growth_rate == gRate]['t'],
                    cladeIstatsTrunc2[cladeIstatsTrunc2.growth_rate == gRate]['richPmean'] -
                    cladeIstatsTrunc2[cladeIstatsTrunc2.growth_rate ==
                                    gRate]['richPstd'],
                    cladeIstatsTrunc2[cladeIstatsTrunc2.growth_rate == gRate]['richPmean'] +
                    cladeIstatsTrunc2[cladeIstatsTrunc2.growth_rate == gRate]['richPstd'], 
                    color=gRateMap(norm(float(i))), alpha=0.3)
for i in range(len(gRates)-1, 0-1, -1):
    print(i)
    gRate = gRates[i]
    axes[1].plot(cladeIstatsTrunc2[cladeIstatsTrunc2.growth_rate == gRate]['t'],
            cladeIstatsTrunc2[cladeIstatsTrunc2.growth_rate ==
                            gRate]['richPmean'],
            color=gRateMap(norm(float(i))), label='{}'.format(gRate), linewidth=2, linestyle='solid',alpha=0.75)
axes[0].set_xscale('log', base=10)
axes[1].set_xscale('log', base=10)
plt.show()


