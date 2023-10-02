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


DBSIM_PATH = os.path.join('/Volumes', 'Yadgah','sylvain-martin-collab/12_MOI3/sweep_db_gathered.sqlite')
DBSIM_PATH2 = os.path.join('/Volumes', 'Yadgah','sylvain-martin-collab/14_MOI3/sweep_db_gathered.sqlite')
DBTREE_PATH = os.path.join('/Volumes', 'Yadgah','sylvain-martin-collab/12_MOI3/mrcaC27.sqlite')
DBTREE_PATH2 = os.path.join('/Volumes', 'Yadgah','sylvain-martin-collab/14_MOI3/treesC3.sqlite')
DBSIM_PATH = os.path.join('/project2/pascualmm/armun/crispr/vary-fitness/12_MOI3/simulation/sweep_db_gathered.sqlite')
DBSIM_PATH2 = os.path.join('/project2/pascualmm/armun/crispr/vary-fitness/14_MOI3/simulation/sweep_db_gathered.sqlite')
DBTREE_PATH = os.path.join('/project2/pascualmm/armun/crispr/vary-fitness/12_MOI3/isolated-runs/isolates/comboID27/mrcaC27.sqlite')
DBTREE_PATH2 = os.path.join('/project2/pascualmm/armun/crispr/vary-fitness/14_MOI3/isolated-runs/isolates/comboID3/treesC3.sqlite')

conSim = sqlite3.connect(DBSIM_PATH)
curSim = conSim.cursor()
conTree = sqlite3.connect(DBTREE_PATH)
curTree = conTree.cursor()
conSim2 = sqlite3.connect(DBSIM_PATH2)
curSim2 = conSim2.cursor()
conTree2 = sqlite3.connect(DBTREE_PATH2)
curTree2 = conTree2.cursor()

cID = 27
# runIDs = pd.read_sql_query(
#     "SELECT combo_id, replicate, run_id FROM runs \
#         WHERE combo_id in ({}) AND replicate = 181"
#     .format(cID), conSim)
runIDs = pd.read_sql_query(
    "SELECT combo_id, replicate, run_id FROM runs \
        WHERE combo_id in ({})"
    .format(cID), conSim)
times = list(range(0, 2001, 5))
# expectedVirusPairwiseMRCA = pd.DataFrame()
# # curTree.execute("CREATE INDEX time_index ON pairwise_MRCA_vstrains (t)")
# for runID in runIDs['run_id']:
#     print(runID)
#     pairwiseVirus = pd.read_sql_query("SELECT run_id, t, t_creation \
#                             FROM pairwise_MRCA_vstrains \
#                             WHERE run_id in ({0}) AND t in ({1})"\
#                             .format(runID, ', '.join(map(str,times))), conTree)
#     pairwiseVirus['t_to_creation'] = pairwiseVirus['t'] - pairwiseVirus['t_creation']
#     expectedVirusPairwiseMRCA = pd.concat([expectedVirusPairwiseMRCA, pairwiseVirus.groupby(['t', 'run_id'])
#                                 .agg(exp_tCreation=('t_creation','mean'),exp_tToCreation=('t_to_creation','mean'))\
#                                 .reset_index()])
        
# conTree.execute('CREATE TABLE expected_pairwise_MRCA_vstrains (run_id, t, exp_tCreation, exp_tToCreation)')
# expectedVirusPairwiseMRCA.to_sql('expected_pairwise_MRCA_vstrains',
#           conTree, if_exists='replace', index=False)
# conTree.commit()
# curTree.execute("CREATE INDEX expected_pairwise_MRCA_vstrains_index \
#                 ON expected_pairwise_MRCA_vstrains (run_id)")
# curTree.execute("CREATE INDEX exp_time_index \
#                 ON expected_pairwise_MRCA_vstrains (t)")

# expectedVirusPairwiseMRCA = pd.read_sql_query("SELECT run_id, t, exp_tCreation, exp_tToCreation \
#                             FROM expected_pairwise_MRCA_vstrains \
#                             WHERE run_id in ({0})".format(', '.join(map(str, runIDs['run_id']))), conTree)

expectedVirusPairwiseMRCA = pd.read_sql_query("SELECT run_id, t, exp_t_creation \
                            FROM expected_pairwise_MRCA_vstrains \
                            WHERE run_id in ({0})".format(', '.join(map(str, runIDs['run_id']))), conTree)
expectedVirusPairwiseMRCA = expectedVirusPairwiseMRCA[expectedVirusPairwiseMRCA.exp_t_creation>=0]
expectedVirusPairwiseMRCA = expectedVirusPairwiseMRCA[expectedVirusPairwiseMRCA['t'].isin(times)]
expectedVirusPairwiseMRCA = expectedVirusPairwiseMRCA.rename(columns={'exp_t_creation':'exp_tToCreation'})
expectedVirusPairwiseMRCA['exp_timeTo'] = expectedVirusPairwiseMRCA['t'] - expectedVirusPairwiseMRCA['exp_tToCreation']
    

microbe_total = pd.read_sql_query("SELECT run_id, t,microbial_abundance FROM summary WHERE run_id in ({})"
                                  .format(', '.join(map(str, runIDs['run_id']))), conSim)\
    .rename(columns={"microbial_abundance": "btotal"})
total = 'btotal'
MRCAbundance = expectedVirusPairwiseMRCA.copy()
# MRCAbundance['exp_timeTo'] = MRCAbundance['t']  - MRCAbundance['exp_tCreation']
MRCAbundance = MRCAbundance.merge(microbe_total,on=['t','run_id'])
MRCAbundance[total] = 1 - MRCAbundance[total]/(10**(5.5))
MRCAbundance[MRCAbundance.t<=200]
base = 2
bins = [base**i for i in range(-12,32+1,1)]
labels = [i for i in range(-12,32,1)]
expMRCAbundance = MRCAbundance.copy()
expMRCAbundance['btotalBinned'] = pd.cut(expMRCAbundance[total], bins=bins, labels=labels)
expMRCAbundance['btotalBinned'] = np.array([float(base)**i for i in np.array(expMRCAbundance['btotalBinned'])])
expMRCAbundance = expMRCAbundance.groupby(['btotalBinned']).agg(exp_timeTo=('exp_timeTo', 'mean'), 
                                                                std_timeTo=('exp_timeTo', 'std'))\
    .reset_index()

######
fig, ax = plt.subplots(1,sharex=True)
axes = [ax, ax.twinx()]
axes[0].plot(virus_total['t'],virus_total['vtotal'],linewidth=0,color='grey')
axes[0].fill_between(virus_total['t'],virus_total['vtotal'], color='grey',alpha=0.6)
axes[0].margins(x=0)
lim = axes[0].get_ylim()
axes[0].set_ylim(0,lim[1])
axes[0].set_yticklabels([])
axes[0].set_yticks([])
axes[1].plot(expectedVirusPairwiseMRCA['t'],\
             expectedVirusPairwiseMRCA['t']/expectedVirusPairwiseMRCA['t'] -
             expectedVirusPairwiseMRCA['exp_tCreation'] /
             expectedVirusPairwiseMRCA['t'],
                color = 'darkred',label='Virus',linewidth=1.5)
# axes[3].plot(expectedVirusPairwiseMRCA['t'],expectedVirusPairwiseMRCA['exp_tCreation'],\
#                 color = 'darkred', label='Virus',linewidth=1)
# axes[1].yaxis.tick_left()
# axes[3].yaxis.set_label_position("left")
axes[1].set_ylabel(ylabel ='Expected Time of MRCA',labelpad=15,fontsize=7)
ax.set_xlabel(xlabel = 'Time t',fontsize=7)
fig.tight_layout()
plt.show()


fig, ax = plt.subplots(1,sharex=True)
ax.scatter(expMRCAbundance['btotalBinned'],expMRCAbundance['exp_timeTo'],c='darkred')
ax.margins(x=0)
ax.set_xscale('log', base=2)
ax.set_yscale('log', base=2)
plt.show()
# ax.set_xscale('log')


######
######
cID2 = 3
# runIDs = pd.read_sql_query(
#     "SELECT combo_id, replicate, run_id FROM runs \
#         WHERE combo_id in ({}) AND replicate = 181"
#     .format(cID), conSim)
runIDs2 = pd.read_sql_query(
    "SELECT combo_id, replicate, run_id FROM runs \
        WHERE combo_id in ({})"
    .format(cID2), conSim2)
expectedVirusPairwiseMRCA2 = pd.DataFrame()
times = list(range(0, 2001, 5))
# curTree.execute("CREATE INDEX time_index ON pairwise_MRCA_vstrains (t)")
for runID in runIDs2['run_id']:
    print(runID)
    pairwiseVirus = pd.read_sql_query("SELECT run_id, t, t_creation \
                            FROM pairwise_MRCA_vstrains \
                            WHERE run_id in ({0}) AND t in ({1})"\
                            .format(runID, ', '.join(map(str,times))), conTree2)
    pairwiseVirus['t_to_creation'] = pairwiseVirus['t'] - pairwiseVirus['t_creation']
    expectedVirusPairwiseMRCA2 = pd.concat([expectedVirusPairwiseMRCA2, pairwiseVirus.groupby(['t', 'run_id'])
                                .agg(exp_tCreation=('t_creation','mean'),exp_tToCreation=('t_to_creation','mean'))\
                                .reset_index()])
        
conTree2.execute('CREATE TABLE expected_pairwise_MRCA_vstrains (run_id, t, exp_tCreation, exp_tToCreation)')
expectedVirusPairwiseMRCA2.to_sql('expected_pairwise_MRCA_vstrains',
          conTree2, if_exists='replace', index=False)
conTree2.commit()
curTree2.execute("CREATE INDEX expected_pairwise_MRCA_vstrains_index \
                ON expected_pairwise_MRCA_vstrains (run_id)")
curTree2.execute("CREATE INDEX exp_time_index \
                ON expected_pairwise_MRCA_vstrains (t)")
# expectedVirusPairwiseMRCA2 = pd.read_sql_query("SELECT run_id, t, exp_tCreation, exp_tToCreation \
#                             FROM expected_pairwise_MRCA_vstrains \
#                             WHERE run_id in ({0})".format(', '.join(map(str, runIDs['run_id']))), conTree2)

microbe_total2 = pd.read_sql_query("SELECT run_id, t,microbial_abundance FROM summary WHERE run_id in ({})"
                                  .format(', '.join(map(str, runIDs2['run_id']))), conSim2)\
    .rename(columns={"microbial_abundance": "btotal"})
total = 'btotal'
MRCAbundance2 = expectedVirusPairwiseMRCA2.copy()
MRCAbundance2['exp_timeTo'] = MRCAbundance2['t']  - MRCAbundance2['exp_tCreation']
MRCAbundance2 = MRCAbundance2.merge(microbe_total2,on=['t','run_id'])
MRCAbundance2[total] = 1 - MRCAbundance2[total]/(10**(5.5))
MRCAbundance2[MRCAbundance2.t<=200]
base = 2
bins = [base**i for i in range(-12,32+1,1)]
labels = [i for i in range(-12,32,1)]
expMRCAbundance2 = MRCAbundance2.copy()
expMRCAbundance2['btotalBinned'] = pd.cut(expMRCAbundance2[total], bins=bins, labels=labels)
expMRCAbundance2['btotalBinned'] = np.array([float(base)**i for i in np.array(expMRCAbundance2['btotalBinned'])])
expMRCAbundance2 = expMRCAbundance2.groupby(['btotalBinned']).agg(exp_timeTo=('exp_timeTo', 'mean'), 
                                                                std_timeTo=('exp_timeTo', 'std'))\
    .reset_index()

######
######


# fig, ax = plt.subplots(1, figsize=(10,10))
# # gs = fig.add_gridspec(1, 2, hspace=0, wspace=.35)
# ax.errorbar(expMRCAbundance2['btotalBinned'], expMRCAbundance2['exp_timeTo'],
#             yerr=expMRCAbundance2['std_timeTo'], fmt='o',
#             color='darkorange', capsize=2, label=r'$\sigma=0$', markersize=8, capthick=2, elinewidth=2)
# ax.errorbar(expMRCAbundance['btotalBinned'], expMRCAbundance['exp_timeTo'],
#             yerr=expMRCAbundance['std_timeTo'], fmt='o',
#             color='mediumblue', capsize=2, label=r'$\sigma=3$', markersize=8, capthick=2, elinewidth=2, alpha=0.6)
# ax.set_xscale('log', base=2)
# ax.set_yscale('log', base=2)
# ax.legend(loc='lower left', fontsize=15)
# # ax[0].set_ylim(.4, 1)
# # ax[0].set_yticks([np.round(i,4) for i in np.arange(.4,1.01,.1)])
# s = ['Expected Time to\n',r'Viral MRCA $t_{MRCA}$']
# ax.set_ylabel(ylabel=r''.join(s),
#                  labelpad=15, fontsize=15)
# ax.set_xlabel(xlabel=r'Host Displacement 1 - ${N/K}$', labelpad=15, fontsize=15)
# ax.tick_params(axis='x', labelsize=15)
# ax.tick_params(axis='y', labelsize=15)
# # fig.tight_layout()
# plt.show()
####
fig = plt.figure(figsize=(20, 8))
gs = fig.add_gridspec(2, 1, hspace=0, wspace=.35)
(ax1, ax2) = gs.subplots(sharex='col')
ax = [ax1, ax2]
ax[1].errorbar(expMRCAbundance['btotalBinned'], expMRCAbundance['exp_timeTo'],
               yerr=expMRCAbundance['std_timeTo'], fmt='o',
               color='mediumblue', capsize=2, label=r'$\sigma=3$', markersize=8, capthick=2, elinewidth=2)
ax[0].errorbar(expMRCAbundance2['btotalBinned'], expMRCAbundance2['exp_timeTo'],
               yerr=expMRCAbundance2['std_timeTo'], fmt='o',
               color='darkorange', capsize=2, label=r'$\sigma=0$', markersize=8, capthick=2, elinewidth=2)
ax[0].set_xscale('log', base=2)
ax[0].set_yscale('log', base=2)
ax[0].legend(loc='lower left', fontsize=15)
ax[1].set_xscale('log', base=2)
ax[1].set_yscale('log', base=2)
ax[1].legend(loc='lower left', fontsize=15)

# ax[0].set_ylim(.4, 1)
# ax[0].set_yticks([np.round(i,4) for i in np.arange(.4,1.01,.1)])
s = ['Expected Time to\n', r'Viral MRCA $t_{MRCA}$']
ax[0].set_ylabel(ylabel=r''.join(s),
                 labelpad=15, fontsize=15)
ax[1].set_xlabel(
    xlabel=r'Host Displacement 1 - ${N/K}$', labelpad=15, fontsize=15)
ax[1].tick_params(axis='x', labelsize=15)
ax[0].tick_params(axis='y', labelsize=15)
ax[1].tick_params(axis='y', labelsize=15)
# fig.tight_layout()
plt.show()
