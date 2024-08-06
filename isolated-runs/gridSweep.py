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
from sci_analysis import analyze




DBCOMBO_PATH = os.path.join('/Volumes/Yadgah/sylvain-martin-collab/25_MOI3/sweep_db.sqlite')
DBSIM_PATH = os.path.join('/Volumes/Yadgah/sylvain-martin-collab/25_MOI3/sweep_db_gathered.sqlite')
DBEXT_PATH = os.path.join('/Volumes/Yadgah/sylvain-martin-collab/25_MOI3/vextinctions-local.sqlite')
DBWALLS_PATH = os.path.join('/Volumes/Yadgah/sylvain-martin-collab/25_MOI3/walls-shannon.sqlite')
simDir = '25_MOI3'
####
conCombo = sqlite3.connect(DBCOMBO_PATH)
curCombo = conCombo.cursor()
conRSim = sqlite3.connect(DBSIM_PATH)
curRSim = conRSim.cursor()
conSim = sqlite3.connect(DBEXT_PATH)
curSim = conSim.cursor()
conWalls = sqlite3.connect(DBWALLS_PATH)
curWalls = conWalls.cursor()
#
viraldiversities = [comm for (comm,) in \
             curCombo.execute("SELECT DISTINCT init_bcomm_function\
                              FROM param_combos \
                              ORDER BY init_bcomm_function").fetchall()]
vparticles = [ind for (ind,) in \
             curCombo.execute("SELECT DISTINCT n_particles_per_vstrain\
                              FROM param_combos \
                              ORDER BY init_bcomm_function").fetchall()]
mutations = [mut for (mut,) in \
             curCombo.execute("SELECT DISTINCT microbe_mutation_prob_per_spacer \
                              FROM param_combos \
                              ORDER BY microbe_mutation_prob_per_spacer").fetchall()]
ccapacities = [cc for (cc,) in \
             curCombo.execute("SELECT DISTINCT microbe_carrying_capacity \
                              FROM param_combos \
                              WHERE microbe_carrying_capacity <= 400000 \
                              ORDER BY microbe_carrying_capacity").fetchall()]
selectionfactors = [cc for (cc,) in \
             curCombo.execute("SELECT DISTINCT evofunctionScale \
                              FROM param_combos \
                              ORDER BY microbe_carrying_capacity").fetchall()]
maxT = max([t for (t,) in \
             curCombo.execute("SELECT DISTINCT t_final \
                              FROM param_combos").fetchall()])
#

treatment = 'Monomorphic'
# thresholds = [lp for (lp,) in curWalls.execute("SELECT DISTINCT lower_percent FROM threshold_values").fetchall()]
thresholdValNumwalls = 0.55
thresholdValDurations = 0.55
maxRuns = 400
resolve = 500
#
for treatment in ['Monomorphic','Polymorphic']:
    if treatment == 'Monomorphic':
        vdiv = 1
    if treatment == 'Polymorphic':
        vdiv = 2
    comboSpace0 = pd.read_sql_query(
        "SELECT combo_id, evofunctionScale \
        FROM param_combos WHERE init_bcomm_function = {0} \
        AND evofunctionScale = 0 \
        AND n_particles_per_vstrain = 100 \
        ORDER BY combo_id" \
            .format(vdiv),conCombo)
    runIDs0 = pd.read_sql_query(
        "SELECT combo_id, run_id FROM runs \
            WHERE combo_id in ({})"\
            .format(', '.join(map(str, comboSpace0['combo_id']))),conCombo)
    timeExt0 = pd.read_sql_query("SELECT t_extinction, run_id FROM vextinctions WHERE run_id in ({})"
                            .format(', '.join(map(str, runIDs0['run_id']))), conSim).merge(runIDs0,on=['run_id'])\
                            .merge(comboSpace0,on=['combo_id']).drop_duplicates()
    vstrains0 = pd.read_sql_query("SELECT vstrain_id, run_id FROM vstrains WHERE run_id in ({})"
                            .format(', '.join(map(str, runIDs0['run_id']))), conRSim).merge(runIDs0,on=['run_id'])\
                            .merge(comboSpace0,on=['combo_id']).drop_duplicates()
    vstrains0 = vstrains0.groupby(['combo_id', 'run_id', 'evofunctionScale']).agg(vstrains=('vstrain_id', 'size')).reset_index() 
    bstrains0 = pd.read_sql_query("SELECT bstrain_id, run_id FROM bstrains WHERE run_id in ({})"
                            .format(', '.join(map(str, runIDs0['run_id']))), conRSim).merge(runIDs0,on=['run_id'])\
                            .merge(comboSpace0,on=['combo_id']).drop_duplicates()
    bstrains0 = bstrains0.groupby(['combo_id', 'run_id', 'evofunctionScale']).agg(bstrains=('bstrain_id', 'size')).reset_index()
    # 
    vstrainsTree0 = pd.read_sql_query("SELECT vstrain_id, parent_vstrain_id, run_id FROM vstrains WHERE run_id in ({})"
                            .format(', '.join(map(str, runIDs0['run_id']))), conRSim).merge(runIDs0,on=['run_id'])\
                            .merge(comboSpace0,on=['combo_id']).drop_duplicates()
    vstrainsTree0 = vstrainsTree0.groupby(['combo_id', 'run_id', 'evofunctionScale', 'parent_vstrain_id']).agg(vstrains=('vstrain_id', 'size')).reset_index() 
    vstrainsTree0 = vstrainsTree0.groupby(['combo_id', 'run_id', 'evofunctionScale']).agg(vstrains=('vstrains', 'mean')).reset_index() 
    #
    bstrainsTree0 = pd.read_sql_query("SELECT bstrain_id, parent_bstrain_id, run_id FROM bstrains WHERE run_id in ({})"
                            .format(', '.join(map(str, runIDs0['run_id']))), conRSim).merge(runIDs0,on=['run_id'])\
                            .merge(comboSpace0,on=['combo_id']).drop_duplicates()
    bstrainsTree0 = bstrainsTree0.groupby(['combo_id', 'run_id', 'evofunctionScale', 'parent_bstrain_id']).agg(bstrains=('bstrain_id', 'size')).reset_index() 
    bstrainsTree0 = bstrainsTree0.groupby(['combo_id', 'run_id', 'evofunctionScale']).agg(bstrains=('bstrains', 'mean')).reset_index() 
    #
    durations0 = pd.read_sql_query("SELECT duration, run_id FROM microbial_peakwall_durations \
                                WHERE run_id in ({0}) AND lower_percent = {1}"
                            .format(', '.join(map(str, runIDs0['run_id'])),thresholdValDurations), conWalls).merge(runIDs0,on=['run_id'])\
                            .merge(comboSpace0,on=['combo_id'])
    #
    mve0 = pd.read_sql_query("SELECT begin_t, run_id FROM microbial_peakwall_durations \
                                WHERE run_id in ({0}) AND upper_percent = 0.98 AND lower_percent = {1}"
                            .format(', '.join(map(str, runIDs0['run_id'])),thresholdValDurations), conWalls).merge(runIDs0,on=['run_id'])\
                            .merge(comboSpace0,on=['combo_id'])
    mveMax = pd.read_sql_query("SELECT end_t, run_id FROM microbial_peakwall_durations \
                                WHERE run_id in ({0}) AND upper_percent = 0.98 AND lower_percent = {1}"
                            .format(', '.join(map(str, runIDs0['run_id'])),thresholdValDurations), conWalls).merge(runIDs0,on=['run_id'])\
                            .merge(comboSpace0,on=['combo_id'])\
                            .groupby(['combo_id', 'run_id', 'evofunctionScale']).agg(begin_t=('end_t', 'max')).reset_index()
    mve0 = pd.concat([mve0,mveMax]).drop_duplicates()
    # duration_0length = list(set(runIDs0['run_id']) - set(np.unique(durations0['run_id'].values)))
    # duration_0length = runIDs0[runIDs0['run_id'].isin(duration_0length)]
    # duration_0length.insert(0, "duration", len(duration_0length)*[0])
    # durations0 = pd.concat([durations0,duration_0length.merge(comboSpace0,on=['combo_id'])]).drop_duplicates() 
    numwalls0 = pd.read_sql_query("SELECT num_peaks, run_id FROM microbial_peakwall_count \
                                WHERE run_id in ({0}) AND upper_percent = 0.98 AND lower_percent = {1}"
                            .format(', '.join(map(str, runIDs0['run_id'])),thresholdValNumwalls), conWalls).merge(runIDs0,on=['run_id'])\
                            .merge(comboSpace0,on=['combo_id']).drop_duplicates()
    numoutbreaks0 = pd.read_sql_query("SELECT num_peaks, run_id FROM microbial_peakwall_count \
                            WHERE run_id in ({0}) AND upper_percent = 0.98 AND lower_percent = {1}"
                        .format(', '.join(map(str, runIDs0['run_id'])),0.85), conWalls).merge(runIDs0,on=['run_id'])\
                        .merge(comboSpace0,on=['combo_id']).drop_duplicates()
    #
    fig, ax = plt.subplots(7, len(mutations),layout="constrained",
                        figsize=(20, 22),sharex=True, sharey='row')
    for j in range(0,len(mutations)):
        mut = mutations[j]
        print(j)
        comboSpace = pd.read_sql_query(
        "SELECT combo_id, evofunctionScale \
        FROM param_combos WHERE init_bcomm_function = {0} \
                AND microbe_mutation_prob_per_spacer = {1} \
                AND evofunctionScale in (0,1,3,6,10) \
                AND n_particles_per_vstrain = 100 \
                ORDER BY combo_id" \
                    .format(vdiv,mut),conCombo)
        runIDs = pd.read_sql_query(
            "SELECT combo_id, run_id FROM runs \
                WHERE combo_id in ({})"\
                .format(', '.join(map(str, comboSpace['combo_id']))),conCombo)
        timeExt = pd.read_sql_query("SELECT t_extinction, run_id FROM vextinctions WHERE run_id in ({})"
                                .format(', '.join(map(str, runIDs['run_id']))), conSim).merge(runIDs,on=['run_id'])\
                                .merge(comboSpace,on=['combo_id']).drop_duplicates()
        timeExt = pd.concat([timeExt,timeExt0]).drop_duplicates()
        timeExt['run_id'] = list(range(1,len(timeExt['run_id'])+1))
        timeExt = timeExt.pivot(index='run_id',columns='evofunctionScale',values='t_extinction')
        boxprops = dict(linewidth=1.5,color='darkblue')
        meanpointprops = dict(marker='D', markeredgecolor='black',
                        markerfacecolor='black') 
        medianprops = dict(linewidth=1.5,color='blue')
        whiskerprops = dict(linewidth=1,color='black')
        boxplot = timeExt.boxplot(ax=ax[0,j],showfliers=False,grid=False, 
                                boxprops=boxprops, medianprops=medianprops, whiskerprops=whiskerprops, whis=[20,80],
                                meanprops=meanpointprops, meanline=False, showmeans=True)
        ax[0,j].tick_params(axis='x', labelsize=15)
        ax[0,j].tick_params(axis='y', labelsize=15)
        ax[0,j].set_title(r'$\mu_s =$ {0}'.format(mutations[j]), fontsize=12)
        ###
        ##
        numoutbreaks = pd.read_sql_query("SELECT num_peaks, run_id FROM microbial_peakwall_count \
                                WHERE run_id in ({0}) AND upper_percent = 0.98 AND lower_percent = {1}"
                            .format(', '.join(map(str, runIDs['run_id'])),0.85), conWalls).merge(runIDs,on=['run_id'])\
                            .merge(comboSpace,on=['combo_id']).drop_duplicates()
        numoutbreaks = pd.concat([numoutbreaks,numoutbreaks0]).drop_duplicates()
        numoutbreaks['run_id'] = list(range(1,len(numoutbreaks['run_id'])+1))
        numoutbreaks = numoutbreaks.pivot(index='run_id',columns='evofunctionScale',values='num_peaks')
        ##
        boxprops = dict(linewidth=1.5,color='purple') 
        medianprops = dict(linewidth=1.5,color='mediumpurple')
        whiskerprops = dict(linewidth=1,color='black')
        boxplotWalls = numoutbreaks.boxplot(ax=ax[1,j],showfliers=False,grid=False, 
                                boxprops=boxprops, medianprops=medianprops, whiskerprops=whiskerprops, whis=[20,80],
                                meanprops=meanpointprops, meanline=False, showmeans=True)
        ax[1,j].tick_params(axis='x', labelsize=15)
        ax[1,j].tick_params(axis='y', labelsize=15)
        ####
        ##
        numwalls = pd.read_sql_query("SELECT num_peaks, run_id FROM microbial_peakwall_count \
                                WHERE run_id in ({0}) AND upper_percent = 0.98 AND lower_percent = {1}"
                            .format(', '.join(map(str, runIDs['run_id'])),thresholdValNumwalls), conWalls).merge(runIDs,on=['run_id'])\
                            .merge(comboSpace,on=['combo_id']).drop_duplicates()
        numwalls = pd.concat([numwalls,numwalls0]).drop_duplicates()
        numwalls['run_id'] = list(range(1,len(numwalls['run_id'])+1))
        numwalls = numwalls.pivot(index='run_id',columns='evofunctionScale',values='num_peaks')
        ##
        boxprops = dict(linewidth=1.5,color='darkgreen') 
        medianprops = dict(linewidth=1.5,color='limegreen')
        whiskerprops = dict(linewidth=1,color='black')
        boxplotWalls = numwalls.boxplot(ax=ax[2,j],showfliers=False,grid=False, 
                                boxprops=boxprops, medianprops=medianprops, whiskerprops=whiskerprops, whis=[20,80],
                                meanprops=meanpointprops, meanline=False, showmeans=True)
        ax[2,j].tick_params(axis='x', labelsize=15)
        ax[2,j].tick_params(axis='y', labelsize=15)
        ###
        ###
        durations = pd.read_sql_query("SELECT duration, run_id FROM microbial_peakwall_durations \
                                WHERE run_id in ({0}) AND upper_percent = 0.98 AND lower_percent = {1}"
                            .format(', '.join(map(str, runIDs['run_id'])),thresholdValDurations), conWalls).merge(runIDs,on=['run_id'])\
                            .merge(comboSpace,on=['combo_id']).drop_duplicates()
        durations = pd.concat([durations,durations0]).drop_duplicates()
        # duration_0length = list(set(runIDs['run_id'].values) - set(np.unique(durations['run_id'].values)))
        # duration_0length = runIDs[runIDs['run_id'].isin(duration_0length)]
        # duration_0length.insert(0, "duration", len(duration_0length)*[0])
        # durations = pd.concat([durations,duration_0length.merge(comboSpace,on=['combo_id'])]).drop_duplicates()
        # .groupby(['combo_id', 't'])\
        durations = durations.groupby(['combo_id', 'run_id', 'evofunctionScale']).agg(duration=('duration', 'mean')).reset_index() 
        durations['run_id'] = list(range(1,len(durations['run_id'])+1))
        durations = durations.pivot(index='run_id',columns='evofunctionScale',values='duration')
        ##
        boxprops = dict(linewidth=1.5,color='darkorange') 
        medianprops = dict(linewidth=1.5,color='orange')
        whiskerprops = dict(linewidth=1,color='black')
        boxplotWalls = durations.boxplot(ax=ax[3,j],showfliers=False,grid=False, 
                                boxprops=boxprops, medianprops=medianprops, whiskerprops=whiskerprops, whis=[20,80],
                                meanprops=meanpointprops, meanline=False, showmeans=True)
        ax[3,j].tick_params(axis='x', labelsize=15)
        ax[3,j].tick_params(axis='y', labelsize=15)
        ###
        ###
        ###
        vstrainsTree = pd.read_sql_query("SELECT vstrain_id, parent_vstrain_id, run_id FROM vstrains WHERE run_id in ({})"
                            .format(', '.join(map(str, runIDs['run_id']))), conRSim).merge(runIDs,on=['run_id'])\
                            .merge(comboSpace,on=['combo_id']).drop_duplicates()
        vstrainsTree = vstrainsTree.groupby(['combo_id', 'run_id', 'evofunctionScale', 'parent_vstrain_id']).agg(vstrains=('vstrain_id', 'size')).reset_index() 
        vstrainsTree = vstrainsTree.groupby(['combo_id', 'run_id', 'evofunctionScale']).agg(vstrains=('vstrains', 'mean')).reset_index()    
        vstrainsTree = pd.concat([vstrainsTree,vstrainsTree0]).drop_duplicates()
        vstrainsTree = vstrainsTree.pivot(index='run_id',columns='evofunctionScale',values='vstrains')
        ##
        boxprops = dict(linewidth=1.5,color='darkcyan') 
        medianprops = dict(linewidth=1.5,color='cyan')
        whiskerprops = dict(linewidth=1,color='black')
        boxplotStrains = vstrainsTree.boxplot(ax=ax[5,j],showfliers=False,grid=False, 
                                boxprops=boxprops, medianprops=medianprops, whiskerprops=whiskerprops, whis=[20,80],
                                meanprops=meanpointprops, meanline=False, showmeans=True)
        ax[5,j].tick_params(axis='x', labelsize=15)
        ax[5,j].tick_params(axis='y', labelsize=15)
        ###
        ###
        bstrainsTree = pd.read_sql_query("SELECT bstrain_id, parent_bstrain_id, run_id FROM bstrains WHERE run_id in ({})"
                            .format(', '.join(map(str, runIDs['run_id']))), conRSim).merge(runIDs,on=['run_id'])\
                            .merge(comboSpace,on=['combo_id']).drop_duplicates()
        bstrainsTree = bstrainsTree.groupby(['combo_id', 'run_id', 'evofunctionScale', 'parent_bstrain_id']).agg(bstrains=('bstrain_id', 'size')).reset_index() 
        bstrainsTree = bstrainsTree.groupby(['combo_id', 'run_id', 'evofunctionScale']).agg(bstrains=('bstrains', 'mean')).reset_index()    
        bstrainsTree = pd.concat([bstrainsTree,bstrainsTree0]).drop_duplicates()
        bstrainsTree = bstrainsTree.pivot(index='run_id',columns='evofunctionScale',values='bstrains')
        ##
        boxprops = dict(linewidth=1.5,color='mediumvioletred') 
        medianprops = dict(linewidth=1.5,color='deeppink')
        whiskerprops = dict(linewidth=1,color='black')
        boxplotStrains = bstrainsTree.boxplot(ax=ax[6,j],showfliers=False,grid=False, 
                                boxprops=boxprops, medianprops=medianprops, whiskerprops=whiskerprops, whis=[20,80],
                                meanprops=meanpointprops, meanline=False, showmeans=True)
        ax[6,j].tick_params(axis='x', labelsize=15)
        ax[6,j].tick_params(axis='y', labelsize=15)
        ###
        ###
        mve = pd.read_sql_query("SELECT begin_t, run_id FROM microbial_peakwall_durations \
                                WHERE run_id in ({0}) AND upper_percent = 0.98 AND lower_percent = {1}"
                            .format(', '.join(map(str, runIDs['run_id'])),thresholdValDurations), conWalls).merge(runIDs,on=['run_id'])\
                            .merge(comboSpace,on=['combo_id'])
        mveMax = pd.read_sql_query("SELECT end_t, run_id FROM microbial_peakwall_durations \
                                    WHERE run_id in ({0}) AND upper_percent = 0.98 AND lower_percent = {1}"
                                .format(', '.join(map(str, runIDs['run_id'])),thresholdValDurations), conWalls).merge(runIDs,on=['run_id'])\
                                .merge(comboSpace,on=['combo_id'])\
                                .groupby(['combo_id', 'run_id', 'evofunctionScale']).agg(begin_t=('end_t', 'max')).reset_index()
        mve = pd.concat([mve,mveMax]).drop_duplicates()
        mve = pd.concat([mve,mve0]).drop_duplicates()
        mve['run_id'] = list(range(1,len(mve['run_id'])+1))
        mve = mve.pivot(index='run_id',columns='evofunctionScale',values='begin_t')
        ###
        boxprops = dict(linewidth=1.5,color='darkred') 
        medianprops = dict(linewidth=1.5,color='red')
        whiskerprops = dict(linewidth=1,color='black')
        boxplotStrains = mve.boxplot(ax=ax[4,j],showfliers=False,grid=False, 
                                boxprops=boxprops, medianprops=medianprops, whiskerprops=whiskerprops, whis=[20,80],
                                meanprops=meanpointprops, meanline=False, showmeans=True)
        ax[4,j].tick_params(axis='x', labelsize=15)
        ax[4,j].tick_params(axis='y', labelsize=15)

    ax[6,0].set_xlabel(xlabel=r'Selection Intensity $\sigma$', labelpad=15, fontsize=15)
    ax[0,0].set_ylabel(ylabel='Ultimate Viral\nExtinction Time', labelpad=15, fontsize=15)
    ax[1,0].set_ylabel(ylabel='Number of\nViral Outbreaks', labelpad=15, fontsize=15)
    ax[2,0].set_ylabel(ylabel='Number of\nMajor Viral Epidemics', labelpad=15, fontsize=15)
    ax[3,0].set_ylabel(ylabel='Duration of\nHost Control', labelpad=15, fontsize=15)
    # ax[3,0].set_ylabel(ylabel=r'Viral Diversity Generated', labelpad=15, fontsize=15)
    # ax[4,0].set_ylabel(ylabel=r'Host Diversity Generated', labelpad=15, fontsize=15)
    ax[5,0].set_ylabel(ylabel='Mean Number of\nViral Mutants\nper Outbreak', labelpad=15, fontsize=15)
    ax[6,0].set_ylabel(ylabel='Mean Number of\nSpacers Acquired\nper Outbreak', labelpad=15, fontsize=15)
    ax[4,0].set_ylabel(ylabel='Time of\nMajor Viral Epidemics', labelpad=15, fontsize=15)
    fig.suptitle('{} Treatment'.format(treatment), fontsize=16)
    fig.savefig(os.path.join('/Volumes/Yadgah/sylvain-martin-collab',simDir,'{0}-extinction_55_quantile_mean_ALL.png'.format(treatment)),dpi=resolve)
    plt.close('all')


plt.show()  

## NEW PAPER RESULT OFFSET 5 SUBPLOTS
plt.close('all')
fig = plt.figure(figsize=(15,8)) 
ax1 = plt.subplot2grid(shape=(2,6), loc=(0,0), colspan=2)
ax2 = plt.subplot2grid((2,6), (0,2), colspan=2)
ax3 = plt.subplot2grid((2,6), (0,4), colspan=2)
ax4 = plt.subplot2grid((2,6), (1,1), colspan=2)
ax5 = plt.subplot2grid((2,6), (1,3), colspan=2)
plt.subplots_adjust(hspace=0.5, wspace=1.8)
thresholdValNumwalls = 0.55
thresholdValDurations = 0.55
maxRuns = 400
resolve = 500
treatment = 'Monomorphic'
vdiv = 1
comboSpace0 = pd.read_sql_query(
    "SELECT combo_id, evofunctionScale \
    FROM param_combos WHERE init_bcomm_function = {0} \
    AND evofunctionScale = 0 \
    AND n_particles_per_vstrain = 100 \
    ORDER BY combo_id" \
        .format(vdiv),conCombo)
runIDs0 = pd.read_sql_query(
    "SELECT combo_id, run_id FROM runs \
        WHERE combo_id in ({})"\
        .format(', '.join(map(str, comboSpace0['combo_id']))),conCombo)
timeExt0 = pd.read_sql_query("SELECT t_extinction, run_id FROM vextinctions WHERE run_id in ({})"
                        .format(', '.join(map(str, runIDs0['run_id']))), conSim).merge(runIDs0,on=['run_id'])\
                        .merge(comboSpace0,on=['combo_id']).drop_duplicates()
vstrains0 = pd.read_sql_query("SELECT vstrain_id, run_id FROM vstrains WHERE run_id in ({})"
                        .format(', '.join(map(str, runIDs0['run_id']))), conRSim).merge(runIDs0,on=['run_id'])\
                        .merge(comboSpace0,on=['combo_id']).drop_duplicates()
vstrains0 = vstrains0.groupby(['combo_id', 'run_id', 'evofunctionScale']).agg(vstrains=('vstrain_id', 'size')).reset_index() 
bstrains0 = pd.read_sql_query("SELECT bstrain_id, run_id FROM bstrains WHERE run_id in ({})"
                        .format(', '.join(map(str, runIDs0['run_id']))), conRSim).merge(runIDs0,on=['run_id'])\
                        .merge(comboSpace0,on=['combo_id']).drop_duplicates()
bstrains0 = bstrains0.groupby(['combo_id', 'run_id', 'evofunctionScale']).agg(bstrains=('bstrain_id', 'size')).reset_index()
# 
vstrainsTree0 = pd.read_sql_query("SELECT vstrain_id, parent_vstrain_id, run_id FROM vstrains WHERE run_id in ({})"
                        .format(', '.join(map(str, runIDs0['run_id']))), conRSim).merge(runIDs0,on=['run_id'])\
                        .merge(comboSpace0,on=['combo_id']).drop_duplicates()
vstrainsTree0 = vstrainsTree0.groupby(['combo_id', 'run_id', 'evofunctionScale', 'parent_vstrain_id']).agg(vstrains=('vstrain_id', 'size')).reset_index() 
vstrainsTree0 = vstrainsTree0.groupby(['combo_id', 'run_id', 'evofunctionScale']).agg(vstrains=('vstrains', 'mean')).reset_index() 
#
bstrainsTree0 = pd.read_sql_query("SELECT bstrain_id, parent_bstrain_id, run_id FROM bstrains WHERE run_id in ({})"
                        .format(', '.join(map(str, runIDs0['run_id']))), conRSim).merge(runIDs0,on=['run_id'])\
                        .merge(comboSpace0,on=['combo_id']).drop_duplicates()
bstrainsTree0 = bstrainsTree0.groupby(['combo_id', 'run_id', 'evofunctionScale', 'parent_bstrain_id']).agg(bstrains=('bstrain_id', 'size')).reset_index() 
bstrainsTree0 = bstrainsTree0.groupby(['combo_id', 'run_id', 'evofunctionScale']).agg(bstrains=('bstrains', 'mean')).reset_index() 
#
durations0 = pd.read_sql_query("SELECT duration, run_id FROM microbial_peakwall_durations \
                            WHERE run_id in ({0}) AND lower_percent = {1}"
                        .format(', '.join(map(str, runIDs0['run_id'])),thresholdValDurations), conWalls).merge(runIDs0,on=['run_id'])\
                        .merge(comboSpace0,on=['combo_id'])
#
numwalls0 = pd.read_sql_query("SELECT num_peaks, run_id FROM microbial_peakwall_count \
                            WHERE run_id in ({0}) AND upper_percent = 0.98 AND lower_percent = {1}"
                        .format(', '.join(map(str, runIDs0['run_id'])),thresholdValNumwalls), conWalls).merge(runIDs0,on=['run_id'])\
                        .merge(comboSpace0,on=['combo_id']).drop_duplicates()
##
mut = 0
comboSpace = pd.read_sql_query(
"SELECT combo_id, evofunctionScale \
FROM param_combos WHERE init_bcomm_function = {0} \
        AND microbe_mutation_prob_per_spacer = {1} \
        AND evofunctionScale in (0,1,3,6,10) \
        AND n_particles_per_vstrain = 100 \
        ORDER BY combo_id" \
            .format(vdiv,mut),conCombo)
runIDs = pd.read_sql_query(
    "SELECT combo_id, run_id FROM runs \
        WHERE combo_id in ({})"\
        .format(', '.join(map(str, comboSpace['combo_id']))),conCombo)
timeExt = pd.read_sql_query("SELECT t_extinction, run_id FROM vextinctions WHERE run_id in ({})"
                        .format(', '.join(map(str, runIDs['run_id']))), conSim).merge(runIDs,on=['run_id'])\
                        .merge(comboSpace,on=['combo_id']).drop_duplicates()
timeExt = pd.concat([timeExt,timeExt0]).drop_duplicates()
timeExt['run_id'] = list(range(1,len(timeExt['run_id'])+1))
timeExt = timeExt.pivot(index='run_id',columns='evofunctionScale',values='t_extinction')
boxprops = dict(linewidth=1.5,color='darkblue')
meanpointprops = dict(marker='D', markeredgecolor='black',
                markerfacecolor='black') 
medianprops = dict(linewidth=1.5,color='blue')
whiskerprops = dict(linewidth=1,color='black')
boxplot = timeExt.boxplot(ax=ax1,showfliers=False,grid=False, 
                        boxprops=boxprops, medianprops=medianprops, whiskerprops=whiskerprops, whis=[20,80],
                        meanprops=meanpointprops, meanline=False, showmeans=True)
ax1.tick_params(axis='x', labelsize=15)
ax1.tick_params(axis='y', labelsize=15)
####
##
numwalls = pd.read_sql_query("SELECT num_peaks, run_id FROM microbial_peakwall_count \
                        WHERE run_id in ({0}) AND upper_percent = 0.98 AND lower_percent = {1}"
                    .format(', '.join(map(str, runIDs['run_id'])),thresholdValNumwalls), conWalls).merge(runIDs,on=['run_id'])\
                    .merge(comboSpace,on=['combo_id']).drop_duplicates()
numwalls = pd.concat([numwalls,numwalls0]).drop_duplicates()
numwalls['run_id'] = list(range(1,len(numwalls['run_id'])+1))
numwalls = numwalls.pivot(index='run_id',columns='evofunctionScale',values='num_peaks')
##
boxprops = dict(linewidth=1.5,color='darkgreen') 
medianprops = dict(linewidth=1.5,color='limegreen')
whiskerprops = dict(linewidth=1,color='black')
boxplotWalls = numwalls.boxplot(ax=ax2,showfliers=False,grid=False, 
                        boxprops=boxprops, medianprops=medianprops, whiskerprops=whiskerprops, whis=[20,80],
                        meanprops=meanpointprops, meanline=False, showmeans=True)
ax2.tick_params(axis='x', labelsize=15)
ax2.tick_params(axis='y', labelsize=15)
###
###
durations = pd.read_sql_query("SELECT duration, run_id FROM microbial_peakwall_durations \
                        WHERE run_id in ({0}) AND upper_percent = 0.98 AND lower_percent = {1}"
                    .format(', '.join(map(str, runIDs['run_id'])),thresholdValDurations), conWalls).merge(runIDs,on=['run_id'])\
                    .merge(comboSpace,on=['combo_id']).drop_duplicates()
durations = pd.concat([durations,durations0]).drop_duplicates()
# duration_0length = list(set(runIDs['run_id'].values) - set(np.unique(durations['run_id'].values)))
# duration_0length = runIDs[runIDs['run_id'].isin(duration_0length)]
# duration_0length.insert(0, "duration", len(duration_0length)*[0])
# durations = pd.concat([durations,duration_0length.merge(comboSpace,on=['combo_id'])]).drop_duplicates()
# .groupby(['combo_id', 't'])\
durations = durations.groupby(['combo_id', 'run_id', 'evofunctionScale']).agg(duration=('duration', 'mean')).reset_index() 
durations['run_id'] = list(range(1,len(durations['run_id'])+1))
durations = durations.pivot(index='run_id',columns='evofunctionScale',values='duration')
##
boxprops = dict(linewidth=1.5,color='darkorange') 
medianprops = dict(linewidth=1.5,color='orange')
whiskerprops = dict(linewidth=1,color='black')
boxplotWalls = durations.boxplot(ax=ax3,showfliers=False,grid=False, 
                        boxprops=boxprops, medianprops=medianprops, whiskerprops=whiskerprops, whis=[20,80],
                        meanprops=meanpointprops, meanline=False, showmeans=True)
ax3.tick_params(axis='x', labelsize=15)
ax3.tick_params(axis='y', labelsize=15)
###
###
###
vstrainsTree = pd.read_sql_query("SELECT vstrain_id, parent_vstrain_id, run_id FROM vstrains WHERE run_id in ({})"
                    .format(', '.join(map(str, runIDs['run_id']))), conRSim).merge(runIDs,on=['run_id'])\
                    .merge(comboSpace,on=['combo_id']).drop_duplicates()
vstrainsTree = vstrainsTree.groupby(['combo_id', 'run_id', 'evofunctionScale', 'parent_vstrain_id']).agg(vstrains=('vstrain_id', 'size')).reset_index() 
vstrainsTree = vstrainsTree.groupby(['combo_id', 'run_id', 'evofunctionScale']).agg(vstrains=('vstrains', 'mean')).reset_index()    
vstrainsTree = pd.concat([vstrainsTree,vstrainsTree0]).drop_duplicates()
ratio = vstrainsTree.copy()
vstrainsTree = vstrainsTree.pivot(index='run_id',columns='evofunctionScale',values='vstrains')
##
boxprops = dict(linewidth=1.5,color='darkcyan') 
medianprops = dict(linewidth=1.5,color='cyan')
whiskerprops = dict(linewidth=1,color='black')
boxplotStrains = vstrainsTree.boxplot(ax=ax4,showfliers=False,grid=False, 
                        boxprops=boxprops, medianprops=medianprops, whiskerprops=whiskerprops, whis=[20,80],
                        meanprops=meanpointprops, meanline=False, showmeans=True)
ax4.tick_params(axis='x', labelsize=15)
ax4.tick_params(axis='y', labelsize=15)
###
###
bstrainsTree = pd.read_sql_query("SELECT bstrain_id, parent_bstrain_id, run_id FROM bstrains WHERE run_id in ({})"
                    .format(', '.join(map(str, runIDs['run_id']))), conRSim).merge(runIDs,on=['run_id'])\
                    .merge(comboSpace,on=['combo_id']).drop_duplicates()
bstrainsTree = bstrainsTree.groupby(['combo_id', 'run_id', 'evofunctionScale', 'parent_bstrain_id']).agg(bstrains=('bstrain_id', 'size')).reset_index() 
bstrainsTree = bstrainsTree.groupby(['combo_id', 'run_id', 'evofunctionScale']).agg(bstrains=('bstrains', 'mean')).reset_index()    
bstrainsTree = pd.concat([bstrainsTree,bstrainsTree0]).drop_duplicates()
ratio = ratio.merge(bstrainsTree,on=['combo_id','run_id','evofunctionScale'])
ratio['ratio'] = ratio['vstrains']/ratio['bstrains']
ratio = ratio.pivot(index='run_id',columns='evofunctionScale',values='ratio')
bstrainsTree = bstrainsTree.pivot(index='run_id',columns='evofunctionScale',values='bstrains')
##
boxprops = dict(linewidth=1.5,color='mediumvioletred') 
medianprops = dict(linewidth=1.5,color='deeppink')
whiskerprops = dict(linewidth=1,color='black')
boxplotStrains = bstrainsTree.boxplot(ax=ax5,showfliers=False,grid=False, 
                        boxprops=boxprops, medianprops=medianprops, whiskerprops=whiskerprops, whis=[20,80],
                        meanprops=meanpointprops, meanline=False, showmeans=True)
ax5.tick_params(axis='x', labelsize=15)
ax5.tick_params(axis='y', labelsize=15)
###


ax1.set_xlabel(xlabel=r'Selection Intensity $\sigma$', labelpad=15, fontsize=15)
ax1.set_ylabel(ylabel='Ultimate Viral\nExtinction Time', labelpad=15, fontsize=15)
ax2.set_ylabel(ylabel='Number of\nMajor Viral Epidemics', labelpad=15, fontsize=15)
ax3.set_ylabel(ylabel='Duration of\nHost Control', labelpad=15, fontsize=15)
ax4.set_ylabel(ylabel='Mean Number of\nViral Mutations\nper Outbreak', labelpad=15, fontsize=15)
ax5.set_ylabel(ylabel='Mean Number of\nSpacers Acquired\nper Outbreak', labelpad=15, fontsize=15)
fig.savefig(os.path.join('/Volumes/Yadgah/sylvain-martin-collab',simDir,'{0}-extinction_55_quantile_mean_ALL_result2_5offset.png'.format(treatment)),dpi=resolve)





## PAPER RESULT!!!!!
thresholdValNumwalls = 0.85
thresholdValDurations = 0.35
maxRuns = 400
resolve = 500
fig, ax = plt.subplots(3,2,layout="constrained",
                    figsize=(10, 12),sharex=True, sharey='row')
for j in [0,1]:
    vdiv = viraldiversities[j]
    mut = 0
    comboSpace0 = pd.read_sql_query(
        "SELECT combo_id, evofunctionScale \
        FROM param_combos WHERE init_bcomm_function = {0} \
        AND evofunctionScale = 0 \
        AND n_particles_per_vstrain = 100 \
        ORDER BY combo_id" \
            .format(vdiv),conCombo)
    runIDs0 = pd.read_sql_query(
        "SELECT combo_id, run_id FROM runs \
            WHERE combo_id in ({})"\
            .format(', '.join(map(str, comboSpace0['combo_id']))),conCombo)
    timeExt0 = pd.read_sql_query("SELECT t_extinction, run_id FROM vextinctions WHERE run_id in ({})"
                            .format(', '.join(map(str, runIDs0['run_id']))), conSim).merge(runIDs0,on=['run_id'])\
                            .merge(comboSpace0,on=['combo_id']).drop_duplicates()
    durations0 = pd.read_sql_query("SELECT duration, run_id FROM microbial_peakwall_durations \
                                WHERE run_id in ({0}) AND lower_percent = {1}"
                            .format(', '.join(map(str, runIDs0['run_id'])),thresholdValDurations), conWalls).merge(runIDs0,on=['run_id'])\
                            .merge(comboSpace0,on=['combo_id'])
    # duration_0length = list(set(runIDs0['run_id']) - set(np.unique(durations0['run_id'].values)))
    # duration_0length = runIDs0[runIDs0['run_id'].isin(duration_0length)]
    # duration_0length.insert(0, "duration", len(duration_0length)*[0])
    # durations0 = pd.concat([durations0,duration_0length.merge(comboSpace0,on=['combo_id'])]).drop_duplicates() 
    numwalls0 = pd.read_sql_query("SELECT num_walls, run_id FROM microbial_peakwall_count \
                                WHERE run_id in ({0}) AND lower_percent = {1}"
                            .format(', '.join(map(str, runIDs0['run_id'])),thresholdValNumwalls), conWalls).merge(runIDs0,on=['run_id'])\
                            .merge(comboSpace0,on=['combo_id']).drop_duplicates()
    ##
    comboSpace = pd.read_sql_query(
    "SELECT combo_id, evofunctionScale \
    FROM param_combos WHERE init_bcomm_function = {0} \
            AND microbe_mutation_prob_per_spacer = {1} \
            AND evofunctionScale in (0,1,3,6,10) \
            AND n_particles_per_vstrain = 100 \
            ORDER BY combo_id" \
                .format(vdiv,mut),conCombo)
    runIDs = pd.read_sql_query(
        "SELECT combo_id, run_id FROM runs \
            WHERE combo_id in ({})"\
            .format(', '.join(map(str, comboSpace['combo_id']))),conCombo)
    timeExt = pd.read_sql_query("SELECT t_extinction, run_id FROM vextinctions WHERE run_id in ({})"
                            .format(', '.join(map(str, runIDs['run_id']))), conSim).merge(runIDs,on=['run_id'])\
                            .merge(comboSpace,on=['combo_id']).drop_duplicates()
    timeExt = pd.concat([timeExt,timeExt0]).drop_duplicates()
    timeExt['run_id'] = list(range(1,len(timeExt['run_id'])+1))
    timeExt = timeExt.pivot(index='run_id',columns='evofunctionScale',values='t_extinction')
    boxprops = dict(linewidth=1.5,color='darkblue')
    meanpointprops = dict(marker='D', markeredgecolor='black',
                    markerfacecolor='black') 
    medianprops = dict(linewidth=2,color='blue')
    whiskerprops = dict(linewidth=1,color='black')
    boxplot = timeExt.boxplot(ax=ax[0,j],grid=False,showfliers=False, 
                            boxprops=boxprops, medianprops=medianprops, whiskerprops=whiskerprops, #whis=[20,80],
                            meanprops=meanpointprops, meanline=False, showmeans=True)
    ax[0,j].tick_params(axis='x', labelsize=15)
    ax[0,j].tick_params(axis='y', labelsize=15)
    ####
    ##
    numwalls = pd.read_sql_query("SELECT num_walls, run_id FROM microbial_peakwall_count \
                            WHERE run_id in ({0}) AND lower_percent = {1}"
                        .format(', '.join(map(str, runIDs['run_id'])),thresholdValNumwalls), conWalls).merge(runIDs,on=['run_id'])\
                        .merge(comboSpace,on=['combo_id']).drop_duplicates()
    numwalls = pd.concat([numwalls,numwalls0]).drop_duplicates()
    numwalls['run_id'] = list(range(1,len(numwalls['run_id'])+1))
    numwalls = numwalls.pivot(index='run_id',columns='evofunctionScale',values='num_walls')
    ##
    boxprops = dict(linewidth=1.5,color='darkgreen') 
    medianprops = dict(linewidth=2,color='limegreen')
    whiskerprops = dict(linewidth=1,color='black')
    boxplotWalls = numwalls.boxplot(ax=ax[1,j],showfliers=False,grid=False, 
                            boxprops=boxprops, medianprops=medianprops, whiskerprops=whiskerprops, #whis=[20,80],
                            meanprops=meanpointprops, meanline=False, showmeans=True)
    ax[1,j].tick_params(axis='x', labelsize=15)
    ax[1,j].tick_params(axis='y', labelsize=15)
    ###
    ###
    durations = pd.read_sql_query("SELECT duration, run_id FROM microbial_peakwall_durations \
                            WHERE run_id in ({0}) AND lower_percent = {1}"
                        .format(', '.join(map(str, runIDs['run_id'])),thresholdValDurations), conWalls).merge(runIDs,on=['run_id'])\
                        .merge(comboSpace,on=['combo_id']).drop_duplicates()
    durations = pd.concat([durations,durations0]).drop_duplicates()
    # duration_0length = list(set(runIDs['run_id'].values) - set(np.unique(durations['run_id'].values)))
    # duration_0length = runIDs[runIDs['run_id'].isin(duration_0length)]
    # duration_0length.insert(0, "duration", len(duration_0length)*[0])
    # durations = pd.concat([durations,duration_0length.merge(comboSpace,on=['combo_id'])]).drop_duplicates() 
    durations['run_id'] = list(range(1,len(durations['run_id'])+1))
    durations = durations.pivot(index='run_id',columns='evofunctionScale',values='duration')
    ##
    boxprops = dict(linewidth=1.5,color='darkred') 
    medianprops = dict(linewidth=2,color='red')
    whiskerprops = dict(linewidth=1,color='black')
    boxplotWalls = durations.boxplot(ax=ax[2,j],showfliers=False,grid=False, 
                            boxprops=boxprops, medianprops=medianprops, whiskerprops=whiskerprops, #whis=[20,80],
                            meanprops=meanpointprops, meanline=False, showmeans=True)
    ax[2,j].tick_params(axis='x', labelsize=15)
    ax[2,j].tick_params(axis='y', labelsize=15)
    #    
    
ax[0,0].set_title('Monomorphic Treatment', fontsize=15)
ax[0,1].set_title('Polymorphic Treatment', fontsize=15)
ax[2,0].set_xlabel(xlabel=r'Selection Intensity $\sigma$', labelpad=15, fontsize=15)
ax[0,0].set_ylabel(ylabel=r'Ultimate Viral Extinction Time $\tau$', labelpad=15, fontsize=15)
ax[1,0].set_ylabel(ylabel=r'Number of\nViral Outbreaks', labelpad=15, fontsize=15)
ax[2,0].set_ylabel(ylabel=r'Duration of\nHost Control', labelpad=15, fontsize=15)
fig.savefig(os.path.join('/Volumes/Yadgah/sylvain-martin-collab',simDir,'mono-poly-extinction_.png'.format(treatment)),dpi=resolve)
plt.show()  



analyze(
    [df['One'], df['Two'], df['Three'], df['Four']],
    groups=['One', 'Two', 'Three', 'Four'],
    categories='Columns',
    title='Unstacked Oneway'
)
analyze(
    timeExt['t_extinction'],
    groups=timeExt['evofunctionScale'],
    categories='evofunctionScale',
    name='t_extinction',
    title='Oneway check'
)

test = timeExt.groupby('evofunctionScale')['t_extinction'].apply(list)

krusk = stats.kruskal(test[0], test[1], test[3], test[6], test[10])

anova = stats.f_oneway(test[0], test[1], test[3], test[6], test[10])

[stats.skew(test[i]) for i in [0,1,3,6,10]]









# ## NEW PAPER RESULT
# plt.close('all')
# thresholdValNumwalls = 0.55
# thresholdValDurations = 0.55
# maxRuns = 400
# resolve = 500
# treatment = 'Monomorphic'
# vdiv = 1
# comboSpace0 = pd.read_sql_query(
#     "SELECT combo_id, evofunctionScale \
#     FROM param_combos WHERE init_bcomm_function = {0} \
#     AND evofunctionScale = 0 \
#     AND n_particles_per_vstrain = 100 \
#     ORDER BY combo_id" \
#         .format(vdiv),conCombo)
# runIDs0 = pd.read_sql_query(
#     "SELECT combo_id, run_id FROM runs \
#         WHERE combo_id in ({})"\
#         .format(', '.join(map(str, comboSpace0['combo_id']))),conCombo)
# timeExt0 = pd.read_sql_query("SELECT t_extinction, run_id FROM vextinctions WHERE run_id in ({})"
#                         .format(', '.join(map(str, runIDs0['run_id']))), conSim).merge(runIDs0,on=['run_id'])\
#                         .merge(comboSpace0,on=['combo_id']).drop_duplicates()
# vstrains0 = pd.read_sql_query("SELECT vstrain_id, run_id FROM vstrains WHERE run_id in ({})"
#                         .format(', '.join(map(str, runIDs0['run_id']))), conRSim).merge(runIDs0,on=['run_id'])\
#                         .merge(comboSpace0,on=['combo_id']).drop_duplicates()
# vstrains0 = vstrains0.groupby(['combo_id', 'run_id', 'evofunctionScale']).agg(vstrains=('vstrain_id', 'size')).reset_index() 
# bstrains0 = pd.read_sql_query("SELECT bstrain_id, run_id FROM bstrains WHERE run_id in ({})"
#                         .format(', '.join(map(str, runIDs0['run_id']))), conRSim).merge(runIDs0,on=['run_id'])\
#                         .merge(comboSpace0,on=['combo_id']).drop_duplicates()
# bstrains0 = bstrains0.groupby(['combo_id', 'run_id', 'evofunctionScale']).agg(bstrains=('bstrain_id', 'size')).reset_index()
# # 
# vstrainsTree0 = pd.read_sql_query("SELECT vstrain_id, parent_vstrain_id, run_id FROM vstrains WHERE run_id in ({})"
#                         .format(', '.join(map(str, runIDs0['run_id']))), conRSim).merge(runIDs0,on=['run_id'])\
#                         .merge(comboSpace0,on=['combo_id']).drop_duplicates()
# vstrainsTree0 = vstrainsTree0.groupby(['combo_id', 'run_id', 'evofunctionScale', 'parent_vstrain_id']).agg(vstrains=('vstrain_id', 'size')).reset_index() 
# vstrainsTree0 = vstrainsTree0.groupby(['combo_id', 'run_id', 'evofunctionScale']).agg(vstrains=('vstrains', 'mean')).reset_index() 
# #
# bstrainsTree0 = pd.read_sql_query("SELECT bstrain_id, parent_bstrain_id, run_id FROM bstrains WHERE run_id in ({})"
#                         .format(', '.join(map(str, runIDs0['run_id']))), conRSim).merge(runIDs0,on=['run_id'])\
#                         .merge(comboSpace0,on=['combo_id']).drop_duplicates()
# bstrainsTree0 = bstrainsTree0.groupby(['combo_id', 'run_id', 'evofunctionScale', 'parent_bstrain_id']).agg(bstrains=('bstrain_id', 'size')).reset_index() 
# bstrainsTree0 = bstrainsTree0.groupby(['combo_id', 'run_id', 'evofunctionScale']).agg(bstrains=('bstrains', 'mean')).reset_index() 
# #
# durations0 = pd.read_sql_query("SELECT duration, run_id FROM microbial_peakwall_durations \
#                             WHERE run_id in ({0}) AND lower_percent = {1}"
#                         .format(', '.join(map(str, runIDs0['run_id'])),thresholdValDurations), conWalls).merge(runIDs0,on=['run_id'])\
#                         .merge(comboSpace0,on=['combo_id'])
# #
# numwalls0 = pd.read_sql_query("SELECT num_peaks, run_id FROM microbial_peakwall_count \
#                             WHERE run_id in ({0}) AND upper_percent = 0.98 AND lower_percent = {1}"
#                         .format(', '.join(map(str, runIDs0['run_id'])),thresholdValNumwalls), conWalls).merge(runIDs0,on=['run_id'])\
#                         .merge(comboSpace0,on=['combo_id']).drop_duplicates()
# fig, ax = plt.subplots(3, 2,layout="constrained",
#                     figsize=(10, 11),sharex=True)
# mut = 0
# comboSpace = pd.read_sql_query(
# "SELECT combo_id, evofunctionScale \
# FROM param_combos WHERE init_bcomm_function = {0} \
#         AND microbe_mutation_prob_per_spacer = {1} \
#         AND evofunctionScale in (0,1,3,6,10) \
#         AND n_particles_per_vstrain = 100 \
#         ORDER BY combo_id" \
#             .format(vdiv,mut),conCombo)
# runIDs = pd.read_sql_query(
#     "SELECT combo_id, run_id FROM runs \
#         WHERE combo_id in ({})"\
#         .format(', '.join(map(str, comboSpace['combo_id']))),conCombo)
# timeExt = pd.read_sql_query("SELECT t_extinction, run_id FROM vextinctions WHERE run_id in ({})"
#                         .format(', '.join(map(str, runIDs['run_id']))), conSim).merge(runIDs,on=['run_id'])\
#                         .merge(comboSpace,on=['combo_id']).drop_duplicates()
# timeExt = pd.concat([timeExt,timeExt0]).drop_duplicates()
# timeExt['run_id'] = list(range(1,len(timeExt['run_id'])+1))
# timeExt = timeExt.pivot(index='run_id',columns='evofunctionScale',values='t_extinction')
# boxprops = dict(linewidth=1.5,color='darkblue')
# meanpointprops = dict(marker='D', markeredgecolor='black',
#                 markerfacecolor='black') 
# medianprops = dict(linewidth=1.5,color='blue')
# whiskerprops = dict(linewidth=1,color='black')
# boxplot = timeExt.boxplot(ax=ax[0,0],showfliers=False,grid=False, 
#                         boxprops=boxprops, medianprops=medianprops, whiskerprops=whiskerprops, whis=[20,80],
#                         meanprops=meanpointprops, meanline=False, showmeans=True)
# ax[0,0].tick_params(axis='x', labelsize=15)
# ax[0,0].tick_params(axis='y', labelsize=15)
# ####
# ##
# numwalls = pd.read_sql_query("SELECT num_peaks, run_id FROM microbial_peakwall_count \
#                         WHERE run_id in ({0}) AND upper_percent = 0.98 AND lower_percent = {1}"
#                     .format(', '.join(map(str, runIDs['run_id'])),thresholdValNumwalls), conWalls).merge(runIDs,on=['run_id'])\
#                     .merge(comboSpace,on=['combo_id']).drop_duplicates()
# numwalls = pd.concat([numwalls,numwalls0]).drop_duplicates()
# numwalls['run_id'] = list(range(1,len(numwalls['run_id'])+1))
# numwalls = numwalls.pivot(index='run_id',columns='evofunctionScale',values='num_peaks')
# ##
# boxprops = dict(linewidth=1.5,color='darkgreen') 
# medianprops = dict(linewidth=1.5,color='limegreen')
# whiskerprops = dict(linewidth=1,color='black')
# boxplotWalls = numwalls.boxplot(ax=ax[0,1],showfliers=False,grid=False, 
#                         boxprops=boxprops, medianprops=medianprops, whiskerprops=whiskerprops, whis=[20,80],
#                         meanprops=meanpointprops, meanline=False, showmeans=True)
# ax[0,1].tick_params(axis='x', labelsize=15)
# ax[0,1].tick_params(axis='y', labelsize=15)
# ###
# ###
# durations = pd.read_sql_query("SELECT duration, run_id FROM microbial_peakwall_durations \
#                         WHERE run_id in ({0}) AND upper_percent = 0.98 AND lower_percent = {1}"
#                     .format(', '.join(map(str, runIDs['run_id'])),thresholdValDurations), conWalls).merge(runIDs,on=['run_id'])\
#                     .merge(comboSpace,on=['combo_id']).drop_duplicates()
# durations = pd.concat([durations,durations0]).drop_duplicates()
# # duration_0length = list(set(runIDs['run_id'].values) - set(np.unique(durations['run_id'].values)))
# # duration_0length = runIDs[runIDs['run_id'].isin(duration_0length)]
# # duration_0length.insert(0, "duration", len(duration_0length)*[0])
# # durations = pd.concat([durations,duration_0length.merge(comboSpace,on=['combo_id'])]).drop_duplicates()
# # .groupby(['combo_id', 't'])\
# durations = durations.groupby(['combo_id', 'run_id', 'evofunctionScale']).agg(duration=('duration', 'mean')).reset_index() 
# durations['run_id'] = list(range(1,len(durations['run_id'])+1))
# durations = durations.pivot(index='run_id',columns='evofunctionScale',values='duration')
# ##
# boxprops = dict(linewidth=1.5,color='darkorange') 
# medianprops = dict(linewidth=1.5,color='orange')
# whiskerprops = dict(linewidth=1,color='black')
# boxplotWalls = durations.boxplot(ax=ax[1,1],showfliers=False,grid=False, 
#                         boxprops=boxprops, medianprops=medianprops, whiskerprops=whiskerprops, whis=[20,80],
#                         meanprops=meanpointprops, meanline=False, showmeans=True)
# ax[1,1].tick_params(axis='x', labelsize=15)
# ax[1,1].tick_params(axis='y', labelsize=15)
# ###
# ###
# ###
# vstrainsTree = pd.read_sql_query("SELECT vstrain_id, parent_vstrain_id, run_id FROM vstrains WHERE run_id in ({})"
#                     .format(', '.join(map(str, runIDs['run_id']))), conRSim).merge(runIDs,on=['run_id'])\
#                     .merge(comboSpace,on=['combo_id']).drop_duplicates()
# vstrainsTree = vstrainsTree.groupby(['combo_id', 'run_id', 'evofunctionScale', 'parent_vstrain_id']).agg(vstrains=('vstrain_id', 'size')).reset_index() 
# vstrainsTree = vstrainsTree.groupby(['combo_id', 'run_id', 'evofunctionScale']).agg(vstrains=('vstrains', 'mean')).reset_index()    
# vstrainsTree = pd.concat([vstrainsTree,vstrainsTree0]).drop_duplicates()
# ratio = vstrainsTree.copy()
# vstrainsTree = vstrainsTree.pivot(index='run_id',columns='evofunctionScale',values='vstrains')
# ##
# boxprops = dict(linewidth=1.5,color='darkcyan') 
# medianprops = dict(linewidth=1.5,color='cyan')
# whiskerprops = dict(linewidth=1,color='black')
# boxplotStrains = vstrainsTree.boxplot(ax=ax[1,0],showfliers=False,grid=False, 
#                         boxprops=boxprops, medianprops=medianprops, whiskerprops=whiskerprops, whis=[20,80],
#                         meanprops=meanpointprops, meanline=False, showmeans=True)
# ax[1,0].tick_params(axis='x', labelsize=15)
# ax[1,0].tick_params(axis='y', labelsize=15)
# ###
# ###
# bstrainsTree = pd.read_sql_query("SELECT bstrain_id, parent_bstrain_id, run_id FROM bstrains WHERE run_id in ({})"
#                     .format(', '.join(map(str, runIDs['run_id']))), conRSim).merge(runIDs,on=['run_id'])\
#                     .merge(comboSpace,on=['combo_id']).drop_duplicates()
# bstrainsTree = bstrainsTree.groupby(['combo_id', 'run_id', 'evofunctionScale', 'parent_bstrain_id']).agg(bstrains=('bstrain_id', 'size')).reset_index() 
# bstrainsTree = bstrainsTree.groupby(['combo_id', 'run_id', 'evofunctionScale']).agg(bstrains=('bstrains', 'mean')).reset_index()    
# bstrainsTree = pd.concat([bstrainsTree,bstrainsTree0]).drop_duplicates()
# ratio = ratio.merge(bstrainsTree,on=['combo_id','run_id','evofunctionScale'])
# ratio['ratio'] = ratio['vstrains']/ratio['bstrains']
# ratio = ratio.pivot(index='run_id',columns='evofunctionScale',values='ratio')
# bstrainsTree = bstrainsTree.pivot(index='run_id',columns='evofunctionScale',values='bstrains')
# ##
# boxprops = dict(linewidth=1.5,color='mediumvioletred') 
# medianprops = dict(linewidth=1.5,color='deeppink')
# whiskerprops = dict(linewidth=1,color='black')
# boxplotStrains = bstrainsTree.boxplot(ax=ax[2,0],showfliers=False,grid=False, 
#                         boxprops=boxprops, medianprops=medianprops, whiskerprops=whiskerprops, whis=[20,80],
#                         meanprops=meanpointprops, meanline=False, showmeans=True)
# ax[2,0].tick_params(axis='x', labelsize=15)
# ax[2,0].tick_params(axis='y', labelsize=15)
# ###
# ###
# ###
# boxprops = dict(linewidth=1.5,color='darkred') 
# medianprops = dict(linewidth=1.5,color='red')
# whiskerprops = dict(linewidth=1,color='black')
# boxplotStrains = ratio.boxplot(ax=ax[2,1],showfliers=False,grid=False, 
#                         boxprops=boxprops, medianprops=medianprops, whiskerprops=whiskerprops, whis=[20,80],
#                         meanprops=meanpointprops, meanline=False, showmeans=True)
# ax[2,1].tick_params(axis='x', labelsize=15)
# ax[2,1].tick_params(axis='y', labelsize=15)

# ax[2,0].set_xlabel(xlabel=r'Selection Intensity $\sigma$', labelpad=15, fontsize=15)
# ax[0,0].set_ylabel(ylabel='Ultimate Viral Extinction Time', labelpad=15, fontsize=15)
# ax[0,1].set_ylabel(ylabel='Number of\nMajor Viral Epidemics', labelpad=15, fontsize=15)
# ax[1,1].set_ylabel(ylabel='Duration of\nHost Control', labelpad=15, fontsize=15)
# ax[1,0].set_ylabel(ylabel='Mean Number of\nViral Mutations\nper Outbreak', labelpad=15, fontsize=15)
# ax[2,0].set_ylabel(ylabel='Mean Number of\nSpacers Acquired\nper Outbreak', labelpad=15, fontsize=15)
# ax[2,1].set_ylabel(ylabel='Mean Ratio of\nDiversity Generated\nper Outbreak', labelpad=15, fontsize=15)
# fig.savefig(os.path.join('/Volumes/Yadgah/sylvain-martin-collab',simDir,'{0}-extinction_55_quantile_mean_ALL_result2.png'.format(treatment)),dpi=resolve)
