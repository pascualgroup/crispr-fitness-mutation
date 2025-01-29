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
import math
import subprocess
import random as rand
import figure_functions as ff

#############
#############

simDir = '26_MOI3'
simDirNoV = '26_MOI3'
DBSIM_PATH = os.path.join('/Volumes/Yadgah/sylvain-martin-collab/{}/sweep_db_gathered.sqlite'.format(simDir))
DBSIM_PATHnoV = os.path.join('/Volumes/Yadgah/sylvain-martin-collab/{}/sweep_db_gathered.sqlite'.format(simDirNoV))
conSim = sqlite3.connect(DBSIM_PATH)
curSim = conSim.cursor()
conSimNoV = sqlite3.connect(DBSIM_PATHnoV)
curSimNoV = conSimNoV.cursor()


numTrees = 1
scale = 3
hosts_per_strain = 100
viruses_per_strain = 100
carryCap = 400000
micMutRep = 0
# combos = [(1,0),(1,1),(2,0),(2,1)]
combos = [(1,0)]
# combos = [(1,0),(2,0)]
vthreshold = 5
bthreshold = 5
##
resolve = 500
figxy = (10, 12)  # setting for tree abundance figure
hratio = [1, 3]  # setting for tree abundance figure
maxticksize = 100  # setting for abundances on individual branches of tree
treepalette = 'turbo'  # Color palette for tree: use contiguous color palette
abundthresholdM = bthreshold/100
abundthresholdV = vthreshold/100
filterStrains = False
spacingM, spacingV = 0, 0
for i in range(0,len(combos)):
    (bcomm,micMutSpacer) = combos[i]
    comboSpace = pd.read_sql_query(
        "SELECT combo_id, evofunctionScale \
            FROM param_combos WHERE \
            microbe_mutation_prob_per_spacer = {0} \
            AND microbe_mutation_prob_per_replication = {1} \
            AND init_bcomm_function = {2} \
            AND n_hosts_per_bstrain >= {3} \
            AND n_particles_per_vstrain = {4} \
            ORDER BY combo_id"
        .format(micMutSpacer, micMutRep, bcomm, hosts_per_strain, viruses_per_strain),
        conSim)
    comboSpaceNoV = pd.read_sql_query(
        "SELECT combo_id, evofunctionScale \
            FROM param_combos WHERE \
            microbe_mutation_prob_per_spacer = 0 \
            AND microbe_mutation_prob_per_replication = {0} \
            AND init_bcomm_function = {1} \
            AND n_hosts_per_bstrain >= {2} \
            AND n_particles_per_vstrain = 0 \
            ORDER BY combo_id"
        .format(micMutRep, bcomm, hosts_per_strain),
        conSimNoV)
    cID = comboSpace[comboSpace.evofunctionScale==scale]['combo_id'].values[0]
    cIDNoV = comboSpaceNoV[comboSpaceNoV.evofunctionScale==scale]['combo_id'].values[0] 
    ##
    runIDs = pd.read_sql_query(
        "SELECT combo_id, run_id FROM runs \
            WHERE combo_id in ({})"
        .format(cID), conSim)
    runIDsNoV = pd.read_sql_query(
        "SELECT combo_id, run_id FROM runs \
            WHERE combo_id in ({})"
        .format(cIDNoV), conSimNoV)
    ##
    for i in range(0,numTrees):
        subprocess.run(["/Volumes/Yadgah/sylvain-martin-collab/trees.jl","{}".format(cID),simDir,"{}".format(bthreshold),"{}".format(vthreshold)]) 
        with open('/Volumes/Yadgah/sylvain-martin-collab/{0}/runID-cID{1}.txt'.format(simDir,cID)) as f:
            run_id = f.readlines()
        ##  
        run_id = int(run_id[0])
        ##
        DBTREE_PATH = os.path.join(
            '/Volumes', 'Yadgah', 'sylvain-martin-collab', simDir,'trees_output_cID{0}-runID{1}-bthresh{2}-vthresh{3}.sqlite'.format(cID,run_id,bthreshold,vthreshold))
        fig, axesTree, figInset, _, speciesColorDict, strainTimes, _, _ = \
            ff.speciesMullerTreeHostVirus(False,run_id, DBSIM_PATH, DBTREE_PATH, maxticksize, 
                                abundthresholdM, abundthresholdV,  filterStrains, 
                                spacingM, spacingV, figxy)
        fig.savefig(os.path.join('/Volumes/Yadgah/sylvain-martin-collab',simDir,'trees-cID{0}-runID{1}.pdf'.format(cID,run_id)),dpi=resolve)
        figInset.savefig(os.path.join('/Volumes/Yadgah/sylvain-martin-collab',simDir,'inset_trees-cID{0}-runID{1}.pdf'.format(cID,run_id)),dpi=resolve)
        plt.close('all')
    # plt.show() 

##

def fgm(s,x):
    return np.exp(-s*x**2)


domain = list(np.linspace(-2, 2, 100))
selection = [3]
alleles = [0, 0.12, 0.24, 0.36, 0.48, 0.6, 0.72, 0.84]
ymax = 1.01
xmax = 1.5
for s in selection:
    fig, ax = plt.subplots(1,figsize=(8,8/1.5))
    ax.plot(domain,list(map(lambda x: fgm(s,x), domain)),\
            color='black',linewidth=6)
    ax.set_ylim(0,ymax)
    ax.set_xlim(-xmax, xmax)
    # ax.legend()
    ax.set_xlabel(xlabel='Trait Allele',labelpad=10,fontsize=30)
    ax.set_ylabel(ylabel='Demographic Rate\n(Competitive Ability)', labelpad=10, fontsize=30)
    ax.tick_params(axis='x', labelsize=30)
    ax.tick_params(axis='y', labelsize=30)
    ax.arrow(0, np.exp(-1/s), -1/s + .075, 0, width = 0.01, ec =None, fc='black')
    ax.text(-1/s*1/2, np.exp(-1/s) - .15, r"$\frac{1}{\sigma}$",fontsize=30, horizontalalignment='center')
    ax.text(1, 0.75, r"$\sigma = 3$",fontsize=30, horizontalalignment='center')  
    # sign = 1
    for i in range(0,len(alleles)):
        
        ax.axvline(x=alleles[i], ymax=fgm(s,alleles[i])/ymax, linestyle='dotted', \
                   color=speciesColorDict[strainTimes[strainTimes.tree_parent_bstrain_id==0]['tree_bstrain_id'].values[i]], linewidth=6)
        # sign = -1 * sign
    fig.tight_layout()
    fig.savefig(os.path.join('/Volumes/Yadgah/sylvain-martin-collab',simDir,'fgmSigma{0}.pdf'.format(s)),dpi=500)





##
numTrees = 1
scale = 3
hosts_per_strain = 100
viruses_per_strain = 100
carryCap = 400000
micMutRep = 0
combos = [(1,0)]
# combos = [(1,0)]
vthreshold = 5
bthreshold = 5
ti = 10
tf = 1000
##
resolve = 500
figxy = (10, 12)  # setting for tree abundance figure
hratio = [1, 3]  # setting for tree abundance figure
maxticksize = 100  # setting for abundances on individual branches of tree
treepalette = 'turbo'  # Color palette for tree: use contiguous color palette
abundthresholdM = bthreshold/100
abundthresholdV = vthreshold/100
filterStrains = False
spacingM, spacingV = 0, 0
rateName = 'Competitive Ability'
diversity = "Immune Richness"
# diversity = "Shannon Immune Diversity"
propTitle = 'Initial Demographic Rates'
strains = 'bstrains'
strain = 'bstrain'
strain_id = 'bstrain_id'
abundanceTitle = 'Host\nAbundance'
strainTitle = 'Host Phylogeny'
strain_abundance = "microbial_abundance"
straintotal = "btotal"
abundance = 'babundance'
tree_strain_id = "tree_bstrain_id"
tree_parent_strain_id = "tree_parent_bstrain_id"
parent_strain_id = "parent_bstrain_id"

for scale in [3]:
    for i in range(0,len(combos)):
        (bcomm,micMutSpacer) = combos[i]
        comboSpace = pd.read_sql_query(
            "SELECT combo_id, evofunctionScale \
                FROM param_combos WHERE \
                microbe_mutation_prob_per_spacer = {0} \
                AND microbe_mutation_prob_per_replication = {1} \
                AND init_bcomm_function = {2} \
                AND n_hosts_per_bstrain >= {3} \
                AND n_particles_per_vstrain = {4} \
                ORDER BY combo_id"
            .format(micMutSpacer, micMutRep, bcomm, hosts_per_strain, viruses_per_strain),
            conSim)
        comboSpaceNoV = pd.read_sql_query(
            "SELECT combo_id, evofunctionScale \
                FROM param_combos WHERE \
                microbe_mutation_prob_per_replication = {0} \
                AND init_bcomm_function = {1} \
                AND n_hosts_per_bstrain >= {2} \
                AND n_particles_per_vstrain = 0 \
                ORDER BY combo_id"
            .format(micMutRep, bcomm, hosts_per_strain),
            conSimNoV)
        ##
        cID = comboSpace[comboSpace.evofunctionScale==scale]['combo_id'].values[0]
        cIDNoV = comboSpaceNoV[comboSpaceNoV.evofunctionScale==scale]['combo_id'].values[0] 
        ##
        runIDs = pd.read_sql_query(
            "SELECT combo_id, run_id FROM runs \
                WHERE combo_id in ({})"
            .format(cID), conSim)
        runIDsNoV = pd.read_sql_query(
            "SELECT combo_id, run_id FROM runs \
                WHERE combo_id in ({})"
            .format(cIDNoV), conSimNoV)
        ##
        with open('/Volumes/Yadgah/sylvain-martin-collab/{0}/all-figures/runID-cID{1}.txt'.format(simDir,cID)) as f:
            run_id = f.readlines()
        ##  
        sampleRunID = int(run_id[0])
        print(cID)
        DBTREE_PATH = os.path.join(
        '/Volumes', 'Yadgah', 'sylvain-martin-collab', simDir,'all-figures','trees_output_cID{0}-runID{1}-bthresh{2}-vthresh{3}.sqlite'\
            .format(cID,sampleRunID,vthreshold,bthreshold))
        _, _, _, _, speciesColorDict, strainTimes, _, microbeStrainsDF = \
            ff.speciesMullerTreeHostVirus(True,sampleRunID, DBSIM_PATH, DBTREE_PATH, maxticksize, 
                                abundthresholdM, abundthresholdV,  filterStrains, 
                                spacingM, spacingV, figxy)
        immuneIDs = ff.findImmuneIDs(runIDs,conSim)
        immuneIDsNoV = ff.findImmuneIDs(runIDsNoV,conSimNoV)
        cladeIDs = ff.findCladeIDs(runIDs,conSim)
        ##########
        noVdist, meanNoVDist = ff.computeTraitDistribution(runIDsNoV,cIDNoV,conSimNoV)
        dist, meanDist = ff.computeTraitDistribution(runIDs,cID,conSim)
        ## DIVERSITY
        diversityStats = ff.computeSpacerDistribution(runIDs,immuneIDs,conSim)
        diversityStatsNoV = ff.computeSpacerDistribution(runIDsNoV,immuneIDsNoV,conSimNoV)
        ## Clade Immune Proportions
        cladeDivAllN = ff.computeCladeDiversity(runIDs,immuneIDs,cladeIDs,conSim)
        cladeDiversity = cladeDivAllN[cladeDivAllN['n']>9]
        #
        vstats = ff.computeVstats(runIDs,conSim)
        ####
        ####
        # FIGURE
        figsizeQuad = (15,7)
        fig = plt.figure(figsize=figsizeQuad)
        gs = plt.GridSpec(3, 2, figure=fig, hspace=.065, wspace=.6, height_ratios=[1,1,4])
        ax1 = fig.add_subplot(gs[0:4, 0])
        ax2 = fig.add_subplot(gs[0:1, 1])
        ax3 = fig.add_subplot(gs[1:2, 1], sharex=ax2)
        ax4 = fig.add_subplot(gs[2:3, 1], sharex=ax2)
        axes = [ax1, ax2, ax3, ax4,
                ax1.twinx(),ax2.twinx(), ax3.twinx()]
        for i in range(0,3):
            # i = 0
            axes[i].fill_between(vstats[vstats.t!=0]['t'],
                                vstats[vstats.t!=0]['exp_vtotal'] -
                                vstats[vstats.t!=0]['std_vtotal'],
                                vstats[vstats.t!=0]['exp_vtotal'] +
                                vstats[vstats.t!=0]['std_vtotal'], color='grey', alpha=0.1)
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
                # print(i)
                axes[i].set_yticks([])
        #
        # axes[0].tick_params(axis='x', labelcolor='w', top=False,
        #                     bottom=False, left=False, right=False)
        axes[1].tick_params(axis='x', labelcolor='w', top=False,
                            bottom=False, left=False, right=False)
        axes[2].tick_params(axis='x', labelcolor='w', top=False,
                        bottom=False, left=False, right=False)
        axes[4].fill_between(meanDist[meanDist['combo_id'] == cID]['t'],
                            meanDist[meanDist['combo_id'] == cID]['meanExp'] -
                            meanDist[meanDist['combo_id'] == cID]['stdExp'],
                            meanDist[meanDist['combo_id'] == cID]['meanExp'] +
                            meanDist[meanDist['combo_id'] == cID]['stdExp'], color='mediumblue', alpha=0.3)
        axes[4].plot(meanDist[meanDist['combo_id'] == cID]['t'],
                    meanDist[meanDist['combo_id'] == cID]['meanExp'],
                    linewidth=2, color='mediumblue', label=' '.join([r'$\sigma =$','{}'.format(scale)]))
        axes[4].fill_between(meanNoVDist[meanNoVDist['combo_id'] == cIDNoV]['t'],
                            meanNoVDist[meanNoVDist['combo_id'] == cIDNoV]['meanExp'] -
                            meanNoVDist[meanNoVDist['combo_id']
                                        == cIDNoV]['stdExp'],
                            meanNoVDist[meanNoVDist['combo_id'] == cIDNoV]['meanExp'] +
                            meanNoVDist[meanNoVDist['combo_id'] == cIDNoV]['stdExp'], color='lime', alpha=0.3)
        axes[4].plot(meanNoVDist[meanNoVDist['combo_id'] == cIDNoV]['t'],
                    meanNoVDist[meanNoVDist['combo_id'] == cIDNoV]['meanExp'],
                    linewidth=2, color='lime', label=' '.join([r'$\sigma =$','{},'.format(scale),'no virus']), linestyle='dashed')
        ##
        cladeIDX = strainTimes[strainTimes[tree_parent_strain_id] == 0]\
                    .merge(microbeStrainsDF,on=['tree_bstrain_id']).sort_values(by=['tree_bstrain_id'])
        cladeIDX = cladeIDX[['tree_bstrain_id','bstrain_id','growth_rate']]
        for cladeID in cladeIDX['bstrain_id']:
            speciesID = microbeStrainsDF[microbeStrainsDF.bstrain_id == cladeID][tree_strain_id].values[0]
            axes[3].fill_between(cladeDiversity[cladeDiversity.clade_id == cladeID]['t'],
                            cladeDiversity[cladeDiversity.clade_id == cladeID]['prop_richness_mean'] -
                            cladeDiversity[cladeDiversity.clade_id ==
                                            cladeID]['prop_richness_std'],
                            cladeDiversity[cladeDiversity.clade_id == cladeID]['prop_richness_mean'] +
                            cladeDiversity[cladeDiversity.clade_id == cladeID]['prop_richness_std'], 
                            color=speciesColorDict[speciesID], alpha=0.3)
        #
        for cladeID in cladeIDX['bstrain_id']:
            gRate = cladeIDX[cladeIDX.bstrain_id == cladeID]['growth_rate'].values[0]
            speciesID = microbeStrainsDF[microbeStrainsDF.bstrain_id == cladeID][tree_strain_id].values[0]
            axes[3].plot(cladeDiversity[cladeDiversity.clade_id == cladeID]['t'],
                    cladeDiversity[cladeDiversity.clade_id ==
                                    cladeID]['prop_richness_mean'],
                    color=speciesColorDict[speciesID], \
                        label=' '.join([r" $\sim$","{:.3e}".format(np.round(gRate,12))]),\
                        linewidth=2, linestyle='solid',alpha=0.75)
        #
        axes[5].fill_between(diversityStats['t'],
                            diversityStats['simpson_mean'] -
                            diversityStats['simpson_std'],
                            diversityStats['simpson_mean'] +
                            diversityStats['simpson_std'], color='mediumblue', alpha=0.3)
        axes[5].plot(diversityStats['t'],
                    diversityStats['simpson_mean'],
                    linewidth=2, color='mediumblue')
        axes[5].fill_between(diversityStatsNoV['t'],
                            diversityStatsNoV['simpson_mean'] -
                            diversityStatsNoV['simpson_std'],
                            diversityStatsNoV['simpson_mean'] +
                            diversityStatsNoV['simpson_std'], color='lime', alpha=0.3)
        axes[5].plot(diversityStatsNoV['t'],
                    diversityStatsNoV['simpson_mean'],
                    linewidth=2, color='lime', linestyle='dashed')
        #
        axes[6].fill_between(diversityStats['t'],
                            diversityStats['richness_mean'] -
                            diversityStats['richness_std'],
                            diversityStats['richness_mean'] +
                            diversityStats['richness_std'], color='mediumblue', alpha=0.3)
        axes[6].plot(diversityStats['t'],
                    diversityStats['richness_mean'],
                    linewidth=2, color='mediumblue')
        axes[6].fill_between(diversityStatsNoV['t'],
                            diversityStatsNoV['richness_mean'] -
                            diversityStatsNoV['richness_std'],
                            diversityStatsNoV['richness_mean'] +
                            diversityStatsNoV['richness_std'], color='lime', alpha=0.3)
        axes[6].plot(diversityStatsNoV['t'],
                    diversityStatsNoV['richness_mean'],
                    linewidth=2, color='lime', linestyle='dashed')
            #
            #
        axes[0].set_xscale('log', base=10)
        axes[1].set_xscale('log', base=10)
        axes[3].tick_params(axis='x', labelsize=15)
        axes[3].tick_params(axis='y', labelsize=15)
        axes[4].tick_params(axis='x', labelsize=15)
        axes[4].tick_params(axis='y', labelsize=15)
        axes[5].tick_params(axis='x', labelsize=15)
        axes[5].tick_params(axis='y', labelsize=10)
        axes[6].tick_params(axis='x', labelsize=15)
        axes[6].tick_params(axis='y', labelsize=10)
        axes[4].yaxis.set_label_position("left")
        axes[4].yaxis.tick_left()
        axes[5].yaxis.set_label_position("left")
        axes[5].yaxis.tick_left()
        axes[6].yaxis.set_label_position("left")
        axes[6].yaxis.tick_left()
        s = ['Expected\n{0} Mean '.format(rateName), r'$\mathbb{E}[\bar{r}]$']
        axes[4].set_ylabel(
            ylabel=r''.join(s), labelpad=10, fontsize=15)
        axes[3].set_ylabel(ylabel='Fraction of\n{}'.format(diversity), labelpad=10, fontsize=15)
        axes[6].set_ylabel(ylabel = 'Immune\nRichness',labelpad=10,fontsize=10)
        axes[5].set_ylabel(ylabel = 'Immune\nSimpson\nIndex',labelpad=10,fontsize=10)
        handles = []
        labels = []
        handle, label = axes[0].get_legend_handles_labels()
        handles.extend(handle)
        labels.extend(label)
        handle, label = axes[4].get_legend_handles_labels()
        handles.extend(handle)
        labels.extend(label)
        axes[4].legend(handles, labels, loc='lower left', fontsize=12)
        handles = []
        labels = []
        handle, label = axes[3].get_legend_handles_labels()
        handles.extend(handle[::-1])
        labels.extend(label[::-1])
        axes[3].legend(handles,labels,loc='upper left', fontsize=10, title=propTitle, title_fontsize=10)
        axes[3].set_xlabel(xlabel='Time', fontsize=15, labelpad=15)
        axes[4].set_xlabel(xlabel='Time', fontsize=15, labelpad=15)
        axes[0].set_xlabel(xlabel='Time', fontsize=15, labelpad=15)
        axes[0].set_xlim(ti, tf)
        axes[1].set_xlim(ti, tf)
        axes[4].set_xlim(ti, tf)
        axes[5].set_xlim(ti, tf)
        axes[2].set_xlim(ti, tf)
        axes[3].set_xlim(ti, tf)
        axes[5].set_xlim(ti, tf)
        if cID == 108:
            lim = axes[3].get_ylim()
            axes[3].set_ylim(0, lim[1])
        fig.savefig(os.path.join('/Volumes/Yadgah/sylvain-martin-collab',simDir,'all-figures','diversity-cID{0}_scale{1}.pdf'.format(cID,scale)),dpi=resolve)
        plt.close('all')
    ##
##
plt.show()


## SUPPLLEMENTARY FIGURE
numTrees = 1
scale = 3
hosts_per_strain = 100
viruses_per_strain = 100
carryCap = 400000
micMutRep = 0
combos = [(1,1),(2,0),(2,1)]
# combos = [(1,0)]
vthreshold = 5
bthreshold = 5
ti = 10
tf = 1000
##
resolve = 500
figxy = (10, 12)  # setting for tree abundance figure
hratio = [1, 3]  # setting for tree abundance figure
maxticksize = 100  # setting for abundances on individual branches of tree
treepalette = 'turbo'  # Color palette for tree: use contiguous color palette
abundthresholdM = bthreshold/100
abundthresholdV = vthreshold/100
filterStrains = False
spacingM, spacingV = 0, 0
rateName = 'Competitive Ability'
diversity = "Immune Richness"
# diversity = "Shannon Immune Diversity"
propTitle = 'Host-intrinsic Growth Rates'
strains = 'bstrains'
strain = 'bstrain'
strain_id = 'bstrain_id'
abundanceTitle = 'Host\nAbundance'
strainTitle = 'Host Phylogeny'
strain_abundance = "microbial_abundance"
straintotal = "btotal"
abundance = 'babundance'
tree_strain_id = "tree_bstrain_id"
tree_parent_strain_id = "tree_parent_bstrain_id"
parent_strain_id = "parent_bstrain_id"

for scale in [3]:
    figsizeQuad = (17,22)
    fig = plt.figure(figsize=figsizeQuad)
    outer_grid = fig.add_gridspec(3, 1, wspace=1.2, hspace=.3)
    for i in range(0,len(combos)):
        (bcomm,micMutSpacer) = combos[i]
        comboSpace = pd.read_sql_query(
            "SELECT combo_id, evofunctionScale \
                FROM param_combos WHERE \
                microbe_mutation_prob_per_spacer = {0} \
                AND microbe_mutation_prob_per_replication = {1} \
                AND init_bcomm_function = {2} \
                AND n_hosts_per_bstrain >= {3} \
                AND n_particles_per_vstrain = {4} \
                ORDER BY combo_id"
            .format(micMutSpacer, micMutRep, bcomm, hosts_per_strain, viruses_per_strain),
            conSim)
        comboSpaceNoV = pd.read_sql_query(
            "SELECT combo_id, evofunctionScale \
                FROM param_combos WHERE \
                microbe_mutation_prob_per_replication = {0} \
                AND init_bcomm_function = {1} \
                AND n_hosts_per_bstrain >= {2} \
                AND n_particles_per_vstrain = 0 \
                ORDER BY combo_id"
            .format(micMutRep, bcomm, hosts_per_strain),
            conSimNoV)
        ##
        cID = comboSpace[comboSpace.evofunctionScale==scale]['combo_id'].values[0]
        cIDNoV = comboSpaceNoV[comboSpaceNoV.evofunctionScale==scale]['combo_id'].values[0] 
        ##
        runIDs = pd.read_sql_query(
            "SELECT combo_id, run_id FROM runs \
                WHERE combo_id in ({})"
            .format(cID), conSim)
        runIDsNoV = pd.read_sql_query(
            "SELECT combo_id, run_id FROM runs \
                WHERE combo_id in ({})"
            .format(cIDNoV), conSimNoV)
        ##
        with open('/Volumes/Yadgah/sylvain-martin-collab/{0}/all-figures/runID-cID{1}.txt'.format(simDir,cID)) as f:
            run_id = f.readlines()
        ##  
        sampleRunID = int(run_id[0])
        print(cID)
        DBTREE_PATH = os.path.join(
        '/Volumes', 'Yadgah', 'sylvain-martin-collab', simDir,'all-figures','trees_output_cID{0}-runID{1}-bthresh{2}-vthresh{3}.sqlite'\
            .format(cID,sampleRunID,vthreshold,bthreshold))
        _, _, _, _, speciesColorDict, strainTimes, _, microbeStrainsDF = \
            ff.speciesMullerTreeHostVirus(True,sampleRunID, DBSIM_PATH, DBTREE_PATH, maxticksize, 
                                abundthresholdM, abundthresholdV,  filterStrains, 
                                spacingM, spacingV, figxy)
        immuneIDs = ff.findImmuneIDs(runIDs,conSim)
        immuneIDsNoV = ff.findImmuneIDs(runIDsNoV,conSimNoV)
        cladeIDs = ff.findCladeIDs(runIDs,conSim)
        ##########
        noVdist, meanNoVDist = ff.computeTraitDistribution(runIDsNoV,cIDNoV,conSimNoV)
        dist, meanDist = ff.computeTraitDistribution(runIDs,cID,conSim)
        ## DIVERSITY
        diversityStats = ff.computeSpacerDistribution(runIDs,immuneIDs,conSim)
        diversityStatsNoV = ff.computeSpacerDistribution(runIDsNoV,immuneIDsNoV,conSimNoV)
        ## Clade Immune Proportions
        cladeDivAllN = ff.computeCladeDiversity(runIDs,immuneIDs,cladeIDs,conSim)
        cladeDiversity = cladeDivAllN[cladeDivAllN['n']>9]
        #
        vstats = ff.computeVstats(runIDs,conSim)
        ####
        ####
        # FIGURE
        inner_grid = outer_grid[i].subgridspec(3, 2, hspace=.095, wspace=.6, height_ratios=[1,1,4])
        ax1 = fig.add_subplot(inner_grid[0:4, 0])
        ax2 = fig.add_subplot(inner_grid[0:1, 1])
        ax3 = fig.add_subplot(inner_grid[1:2, 1], sharex=ax2)
        ax4 = fig.add_subplot(inner_grid[2:3, 1], sharex=ax2)
        axes = [ax1, ax2, ax3, ax4,
                ax1.twinx(),ax2.twinx(), ax3.twinx()]
        for j in range(0,3):
            # i = 0
            axes[j].fill_between(vstats[vstats.t!=0]['t'],
                                vstats[vstats.t!=0]['exp_vtotal'] -
                                vstats[vstats.t!=0]['std_vtotal'],
                                vstats[vstats.t!=0]['exp_vtotal'] +
                                vstats[vstats.t!=0]['std_vtotal'], color='grey', alpha=0.1)
            axes[j].plot(vstats['t'],
                        vstats['exp_vtotal'],
                        linewidth=2, color='grey', label='Viral Abund.', linestyle='solid', alpha=0.75)
            lim = axes[j].get_ylim()
            axes[j].set_ylim(0, lim[1])
            axes[j].yaxis.tick_right()
            # axes[i].set_xlim(0, 2000)
            axes[j].tick_params(axis='x', labelsize=15)
            # axes[i].set_xscale('log', base=10)
            if j == 0:
                axes[j].set_ylabel(ylabel='Viral Abundance',
                                rotation=270, labelpad=25, fontsize=15)
                axes[j].tick_params(axis='y', labelsize=15)
                axes[j].yaxis.get_offset_text().set_fontsize(15)
                axes[j].yaxis.set_label_position("right")
            if j != 0:
                # print(i)
                axes[j].set_yticks([])
        #
        # axes[0].tick_params(axis='x', labelcolor='w', top=False,
        #                     bottom=False, left=False, right=False)
        axes[1].tick_params(axis='x', labelcolor='w', top=False,
                            bottom=False, left=False, right=False)
        axes[2].tick_params(axis='x', labelcolor='w', top=False,
                        bottom=False, left=False, right=False)
        axes[4].fill_between(meanDist[meanDist['combo_id'] == cID]['t'],
                            meanDist[meanDist['combo_id'] == cID]['meanExp'] -
                            meanDist[meanDist['combo_id'] == cID]['stdExp'],
                            meanDist[meanDist['combo_id'] == cID]['meanExp'] +
                            meanDist[meanDist['combo_id'] == cID]['stdExp'], color='mediumblue', alpha=0.3)
        axes[4].plot(meanDist[meanDist['combo_id'] == cID]['t'],
                    meanDist[meanDist['combo_id'] == cID]['meanExp'],
                    linewidth=2, color='mediumblue', label=' '.join([r'$\sigma =$','{}'.format(scale)]))
        axes[4].fill_between(meanNoVDist[meanNoVDist['combo_id'] == cIDNoV]['t'],
                            meanNoVDist[meanNoVDist['combo_id'] == cIDNoV]['meanExp'] -
                            meanNoVDist[meanNoVDist['combo_id']
                                        == cIDNoV]['stdExp'],
                            meanNoVDist[meanNoVDist['combo_id'] == cIDNoV]['meanExp'] +
                            meanNoVDist[meanNoVDist['combo_id'] == cIDNoV]['stdExp'], color='lime', alpha=0.3)
        axes[4].plot(meanNoVDist[meanNoVDist['combo_id'] == cIDNoV]['t'],
                    meanNoVDist[meanNoVDist['combo_id'] == cIDNoV]['meanExp'],
                    linewidth=2, color='lime', label=' '.join([r'$\sigma =$','{},'.format(scale),'no virus']), linestyle='dashed')
        ##
        cladeIDX = strainTimes[strainTimes[tree_parent_strain_id] == 0]\
                    .merge(microbeStrainsDF,on=['tree_bstrain_id']).sort_values(by=['tree_bstrain_id'])
        cladeIDX = cladeIDX[['tree_bstrain_id','bstrain_id','growth_rate']]
        for cladeID in cladeIDX['bstrain_id']:
            speciesID = microbeStrainsDF[microbeStrainsDF.bstrain_id == cladeID][tree_strain_id].values[0]
            axes[3].fill_between(cladeDiversity[cladeDiversity.clade_id == cladeID]['t'],
                            cladeDiversity[cladeDiversity.clade_id == cladeID]['prop_richness_mean'] -
                            cladeDiversity[cladeDiversity.clade_id ==
                                            cladeID]['prop_richness_std'],
                            cladeDiversity[cladeDiversity.clade_id == cladeID]['prop_richness_mean'] +
                            cladeDiversity[cladeDiversity.clade_id == cladeID]['prop_richness_std'], 
                            color=speciesColorDict[speciesID], alpha=0.3)
        #
        for cladeID in cladeIDX['bstrain_id']:
            gRate = cladeIDX[cladeIDX.bstrain_id == cladeID]['growth_rate'].values[0]
            speciesID = microbeStrainsDF[microbeStrainsDF.bstrain_id == cladeID][tree_strain_id].values[0]
            axes[3].plot(cladeDiversity[cladeDiversity.clade_id == cladeID]['t'],
                    cladeDiversity[cladeDiversity.clade_id ==
                                    cladeID]['prop_richness_mean'],
                    color=speciesColorDict[speciesID], \
                        label=' '.join([r" $\sim$","{:.3e}".format(np.round(gRate,12))]),\
                        linewidth=2, linestyle='solid',alpha=0.75)
        #
        axes[5].fill_between(diversityStats['t'],
                            diversityStats['simpson_mean'] -
                            diversityStats['simpson_std'],
                            diversityStats['simpson_mean'] +
                            diversityStats['simpson_std'], color='mediumblue', alpha=0.3)
        axes[5].plot(diversityStats['t'],
                    diversityStats['simpson_mean'],
                    linewidth=2, color='mediumblue')
        axes[5].fill_between(diversityStatsNoV['t'],
                            diversityStatsNoV['simpson_mean'] -
                            diversityStatsNoV['simpson_std'],
                            diversityStatsNoV['simpson_mean'] +
                            diversityStatsNoV['simpson_std'], color='lime', alpha=0.3)
        axes[5].plot(diversityStatsNoV['t'],
                    diversityStatsNoV['simpson_mean'],
                    linewidth=2, color='lime', linestyle='dashed')
        #
        axes[6].fill_between(diversityStats['t'],
                            diversityStats['richness_mean'] -
                            diversityStats['richness_std'],
                            diversityStats['richness_mean'] +
                            diversityStats['richness_std'], color='mediumblue', alpha=0.3)
        axes[6].plot(diversityStats['t'],
                    diversityStats['richness_mean'],
                    linewidth=2, color='mediumblue')
        axes[6].fill_between(diversityStatsNoV['t'],
                            diversityStatsNoV['richness_mean'] -
                            diversityStatsNoV['richness_std'],
                            diversityStatsNoV['richness_mean'] +
                            diversityStatsNoV['richness_std'], color='lime', alpha=0.3)
        axes[6].plot(diversityStatsNoV['t'],
                    diversityStatsNoV['richness_mean'],
                    linewidth=2, color='lime', linestyle='dashed')
            #
            #
        axes[0].set_xscale('log', base=10)
        axes[1].set_xscale('log', base=10)
        axes[3].tick_params(axis='x', labelsize=15)
        axes[3].tick_params(axis='y', labelsize=15)
        axes[4].tick_params(axis='x', labelsize=15)
        axes[4].tick_params(axis='y', labelsize=15)
        axes[5].tick_params(axis='x', labelsize=15)
        axes[5].tick_params(axis='y', labelsize=10)
        axes[6].tick_params(axis='x', labelsize=15)
        axes[6].tick_params(axis='y', labelsize=10)
        axes[4].yaxis.set_label_position("left")
        axes[4].yaxis.tick_left()
        axes[5].yaxis.set_label_position("left")
        axes[5].yaxis.tick_left()
        axes[6].yaxis.set_label_position("left")
        axes[6].yaxis.tick_left()
        s = ['Expected\n{0} Mean '.format(rateName), r'$\mathbb{E}[\bar{r}]$']
        axes[4].set_ylabel(
            ylabel=r''.join(s), labelpad=10, fontsize=15)
        axes[3].set_ylabel(ylabel='Fraction of\n{}'.format(diversity), labelpad=10, fontsize=15)
        axes[6].set_ylabel(ylabel = 'Immune\nRichness',labelpad=10,fontsize=10)
        axes[5].set_ylabel(ylabel = 'Immune\nSimpson\nIndex',labelpad=10,fontsize=10)

        handles = []
        labels = []
        handle, label = axes[0].get_legend_handles_labels()
        handles.extend(handle)
        labels.extend(label)
        handle, label = axes[4].get_legend_handles_labels()
        handles.extend(handle)
        labels.extend(label)
        axes[4].legend(handles, labels, loc='lower left', fontsize=12)
        handles = []
        labels = []
        handle, label = axes[3].get_legend_handles_labels()
        handles.extend(handle[::-1])
        labels.extend(label[::-1])
        axes[3].legend(handles,labels,loc='upper left', fontsize=10, title=propTitle, title_fontsize=10)        
        axes[3].set_xlabel(xlabel='Time', fontsize=15, labelpad=8)
        axes[4].set_xlabel(xlabel='Time', fontsize=15, labelpad=8)
        axes[0].set_xlabel(xlabel='Time', fontsize=15, labelpad=8)
        axes[0].set_xlim(ti, tf)
        axes[1].set_xlim(ti, tf)
        axes[4].set_xlim(ti, tf)
        axes[5].set_xlim(ti, tf)
        axes[2].set_xlim(ti, tf)
        axes[3].set_xlim(ti, tf)
        axes[5].set_xlim(ti, tf)
        if cID == 108:
            lim = axes[3].get_ylim()
            axes[3].set_ylim(0, lim[1])
        axes[4].text(10.5, 1, r"$\mu_s =$ {}".format(micMutSpacer),
                fontsize=20, horizontalalignment='left', verticalalignment='center')
    fig.savefig(os.path.join('/Volumes/Yadgah/sylvain-martin-collab',simDir,'all-figures','supp_diversities-scale{0}.pdf'.format(cID)),dpi=resolve)
    plt.close('all')
    ##
##







##
## SHANNON!
numTrees = 3
scale = 12
hosts_per_strain = 100
viruses_per_strain = 100
carryCap = 400000
micMutRep = 0
combos = [(1,0),(1,1),(2,0),(2,1)]
rateName = 'Competitive Ability'
diversity = "Immune Richness"
# diversity = "Shannon Immune Diversity"
propTitle = 'Initial Demographic Rates'
strains = 'bstrains'
strain = 'bstrain'
strain_id = 'bstrain_id'
abundanceTitle = 'Host\nAbundance'
strainTitle = 'Host Phylogeny'
strain_abundance = "microbial_abundance"
straintotal = "btotal"
abundance = 'babundance'
tree_strain_id = "tree_bstrain_id"
tree_parent_strain_id = "tree_parent_bstrain_id"
parent_strain_id = "parent_bstrain_id"
for i in range(0,len(combos)):
    (bcomm,micMutSpacer) = combos[i]
    comboSpace = pd.read_sql_query(
        "SELECT combo_id, evofunctionScale \
            FROM param_combos WHERE \
            microbe_mutation_prob_per_spacer = {0} \
            AND microbe_mutation_prob_per_replication = {1} \
            AND init_bcomm_function = {2} \
            AND n_hosts_per_bstrain >= {3} \
            AND n_particles_per_vstrain = {4} \
            ORDER BY combo_id"
        .format(micMutSpacer, micMutRep, bcomm, hosts_per_strain, viruses_per_strain),
        conSim)
    comboSpaceNoV = pd.read_sql_query(
        "SELECT combo_id, evofunctionScale \
            FROM param_combos WHERE \
            microbe_mutation_prob_per_replication = {0} \
            AND init_bcomm_function = {1} \
            AND n_hosts_per_bstrain >= {2} \
            AND n_particles_per_vstrain = 0 \
            ORDER BY combo_id"
        .format(micMutRep, bcomm, hosts_per_strain),
        conSimNoV)
    ##
    cID = comboSpace[comboSpace.evofunctionScale==scale]['combo_id'].values[0]
    cIDNoV = comboSpaceNoV[comboSpaceNoV.evofunctionScale==scale]['combo_id'].values[0] 
    ##
    runIDs = pd.read_sql_query(
        "SELECT combo_id, run_id FROM runs \
            WHERE combo_id in ({})"
        .format(cID), conSim)
    runIDsNoV = pd.read_sql_query(
        "SELECT combo_id, run_id FROM runs \
            WHERE combo_id in ({})"
        .format(cIDNoV), conSimNoV)
    ##
    with open('/Volumes/Yadgah/sylvain-martin-collab/{0}/runID-cID{1}.txt'.format(simDir,cID)) as f:
        run_id = f.readlines()
    ##  
    sampleRunID = int(run_id[0])
    print(cID)
    DBTREE_PATH = os.path.join(
    '/Volumes', 'Yadgah', 'sylvain-martin-collab', simDir,'trees_output_cID{0}-runID{1}-bthresh{2}-vthresh{3}.sqlite'\
        .format(cID,sampleRunID,vthreshold,bthreshold))
    _, _, _, _, speciesColorDict, strainTimes, _, microbeStrainsDF = \
        ff.speciesMullerTreeHostVirus(True,sampleRunID, DBSIM_PATH, DBTREE_PATH, maxticksize, 
                            abundthresholdM, abundthresholdV,  filterStrains, 
                            spacingM, spacingV, figxy)
    immuneIDs = ff.findImmuneIDs(runIDs,conSim)
    immuneIDsNoV = ff.findImmuneIDs(runIDsNoV,conSimNoV)
    cladeIDs = ff.findCladeIDs(runIDs,conSim)
    ##########
    noVdist, meanNoVDist = ff.computeTraitDistribution(runIDsNoV,cIDNoV,conSimNoV)
    dist, meanDist = ff.computeTraitDistribution(runIDs,cID,conSim)
    ## DIVERSITY
    diversityStats = ff.computeSpacerDistribution(runIDs,immuneIDs,conSim)
    diversityStatsNoV = ff.computeSpacerDistribution(runIDsNoV,immuneIDsNoV,conSimNoV)
    ## Clade Immune Proportions
    cladeDiversity = ff.computeCladeDiversity(runIDs,immuneIDs,cladeIDs,conSim)
    #
    vstats = ff.computeVstats(runIDs,conSim)
    ####
    ####
    # FIGURE
    figsizeQuad = (15,10)
    fig = plt.figure(figsize=figsizeQuad)
    gs = plt.GridSpec(2, 2, figure=fig, hspace=.065, wspace=.6)
    ax1 = fig.add_subplot(gs[0:1, 0])
    ax2 = fig.add_subplot(gs[0:1, 1])
    ax3 = fig.add_subplot(gs[1:2, 0], sharex=ax1)
    ax4 = fig.add_subplot(gs[1:2, 1], sharex=ax2)
    axes = [ax1, ax2, ax3, ax4,
            ax1.twinx(),ax2.twinx(),ax3.twinx()]
    for i in range(0,3):
        # i = 0
        axes[i].fill_between(vstats[vstats.t!=0]['t'],
                            vstats[vstats.t!=0]['exp_vtotal'] -
                            vstats[vstats.t!=0]['std_vtotal'],
                            vstats[vstats.t!=0]['exp_vtotal'] +
                            vstats[vstats.t!=0]['std_vtotal'], color='grey', alpha=0.1)
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
            # print(i)
            axes[i].set_yticks([])
    #
    axes[0].tick_params(axis='x', labelcolor='w', top=False,
                        bottom=False, left=False, right=False)
    axes[1].tick_params(axis='x', labelcolor='w', top=False,
                        bottom=False, left=False, right=False)
    axes[4].fill_between(meanDist[meanDist['combo_id'] == cID]['t'],
                        meanDist[meanDist['combo_id'] == cID]['meanExp'] -
                        meanDist[meanDist['combo_id'] == cID]['stdExp'],
                        meanDist[meanDist['combo_id'] == cID]['meanExp'] +
                        meanDist[meanDist['combo_id'] == cID]['stdExp'], color='mediumblue', alpha=0.3)
    axes[4].plot(meanDist[meanDist['combo_id'] == cID]['t'],
                meanDist[meanDist['combo_id'] == cID]['meanExp'],
                linewidth=2, color='mediumblue', label=' '.join([r'$\sigma =$','{}'.format(scale)]))
    axes[4].fill_between(meanNoVDist[meanNoVDist['combo_id'] == cIDNoV]['t'],
                        meanNoVDist[meanNoVDist['combo_id'] == cIDNoV]['meanExp'] -
                        meanNoVDist[meanNoVDist['combo_id']
                                    == cIDNoV]['stdExp'],
                        meanNoVDist[meanNoVDist['combo_id'] == cIDNoV]['meanExp'] +
                        meanNoVDist[meanNoVDist['combo_id'] == cIDNoV]['stdExp'], color='lime', alpha=0.3)
    axes[4].plot(meanNoVDist[meanNoVDist['combo_id'] == cIDNoV]['t'],
                meanNoVDist[meanNoVDist['combo_id'] == cIDNoV]['meanExp'],
                linewidth=2, color='lime', label=' '.join([r'$\sigma =$','{},'.format(scale),'no virus']), linestyle='dashed')
    axes[6].fill_between(meanDist[meanDist['combo_id'] == cID]['t'],
                        meanDist[meanDist['combo_id'] == cID]['meanStd'] -
                        meanDist[meanDist['combo_id'] == cID]['stdStd'],
                        meanDist[meanDist['combo_id'] == cID]['meanStd'] +
                        meanDist[meanDist['combo_id'] == cID]['stdStd'], color='mediumblue', alpha=0.3)
    axes[6].plot(meanDist[meanDist['combo_id'] == cID]['t'],
                meanDist[meanDist['combo_id'] == cID]['meanStd'],
                linewidth=2, color='mediumblue', label=' '.join([r'$\sigma =$','{}'.format(scale)]))
    axes[6].fill_between(meanNoVDist[meanNoVDist['combo_id'] == cIDNoV]['t'],
                        meanNoVDist[meanNoVDist['combo_id'] == cIDNoV]['meanStd'] -
                        meanNoVDist[meanNoVDist['combo_id']
                                    == cIDNoV]['stdStd'],
                        meanNoVDist[meanNoVDist['combo_id'] == cIDNoV]['meanStd'] +
                        meanNoVDist[meanNoVDist['combo_id'] == cIDNoV]['stdStd'], color='lime', alpha=0.3)
    axes[6].plot(meanNoVDist[meanNoVDist['combo_id'] == cIDNoV]['t'],
                meanNoVDist[meanNoVDist['combo_id'] == cIDNoV]['meanStd'],
                linewidth=2, color='lime', label= ' '.join([r'$\sigma =$','{},'.format(scale),'no virus']), linestyle='dashed')
    ##
    if diversity == "Immune Richness":
        cladeIDX = strainTimes[strainTimes[tree_parent_strain_id] == 0]\
                    .merge(microbeStrainsDF,on=['tree_bstrain_id']).sort_values(by=['tree_bstrain_id'])
        cladeIDX = cladeIDX[['tree_bstrain_id','bstrain_id','growth_rate']]
        for cladeID in cladeIDX['bstrain_id']:
            speciesID = microbeStrainsDF[microbeStrainsDF.bstrain_id == cladeID][tree_strain_id].values[0]
            axes[3].fill_between(cladeDiversity[cladeDiversity.clade_id == cladeID]['t'],
                            cladeDiversity[cladeDiversity.clade_id == cladeID]['prop_richness_mean'] -
                            cladeDiversity[cladeDiversity.clade_id ==
                                            cladeID]['prop_richness_std'],
                            cladeDiversity[cladeDiversity.clade_id == cladeID]['prop_richness_mean'] +
                            cladeDiversity[cladeDiversity.clade_id == cladeID]['prop_richness_std'], 
                            color=speciesColorDict[speciesID], alpha=0.3)
        #
        for cladeID in cladeIDX['bstrain_id']:
            gRate = cladeIDX[cladeIDX.bstrain_id == cladeID]['growth_rate'].values[0]
            speciesID = microbeStrainsDF[microbeStrainsDF.bstrain_id == cladeID][tree_strain_id].values[0]
            axes[3].plot(cladeDiversity[cladeDiversity.clade_id == cladeID]['t'],
                    cladeDiversity[cladeDiversity.clade_id ==
                                    cladeID]['prop_richness_mean'],
                    color=speciesColorDict[speciesID], \
                        label=' '.join([r" $\sim$","{:.3e}".format(np.round(gRate,12))]),\
                        linewidth=2, linestyle='solid',alpha=0.75)
        #
        axes[5].fill_between(diversityStats['t'],
                            diversityStats['richness_mean'] -
                            diversityStats['richness_std'],
                            diversityStats['richness_mean'] +
                            diversityStats['richness_std'], color='mediumblue', alpha=0.3)
        axes[5].plot(diversityStats['t'],
                    diversityStats['richness_mean'],
                    linewidth=2, color='mediumblue')
        axes[5].fill_between(diversityStatsNoV['t'],
                            diversityStatsNoV['richness_mean'] -
                            diversityStatsNoV['richness_std'],
                            diversityStatsNoV['richness_mean'] +
                            diversityStatsNoV['richness_std'], color='lime', alpha=0.3)
        axes[5].plot(diversityStatsNoV['t'],
                    diversityStatsNoV['richness_mean'],
                    linewidth=2, color='lime', linestyle='dashed')
        #
    if diversity == "Shannon Immune Diversity":
        cladeIDX = strainTimes[strainTimes[tree_parent_strain_id] == 0]\
                    .merge(microbeStrainsDF,on=['tree_bstrain_id']).sort_values(by=['tree_bstrain_id'])
        cladeIDX = cladeIDX[['tree_bstrain_id','bstrain_id','growth_rate']]
        for cladeID in cladeIDX['bstrain_id']:
            speciesID = microbeStrainsDF[microbeStrainsDF.bstrain_id == cladeID][tree_strain_id].values[0]
            axes[3].fill_between(cladeDiversity[cladeDiversity.clade_id == cladeID]['t'],
                            cladeDiversity[cladeDiversity.clade_id == cladeID]['prop_shannon_mean'] -
                            cladeDiversity[cladeDiversity.clade_id ==
                                            cladeID]['prop_shannon_std'],
                            cladeDiversity[cladeDiversity.clade_id == cladeID]['prop_shannon_mean'] +
                            cladeDiversity[cladeDiversity.clade_id == cladeID]['prop_shannon_std'], 
                            color=speciesColorDict[speciesID], alpha=0.3)
        #
        for cladeID in cladeIDX['bstrain_id']:
            gRate = cladeIDX[cladeIDX.bstrain_id == cladeID]['growth_rate'].values[0]
            speciesID = microbeStrainsDF[microbeStrainsDF.bstrain_id == cladeID][tree_strain_id].values[0]
            axes[3].plot(cladeDiversity[cladeDiversity.clade_id == cladeID]['t'],
                    cladeDiversity[cladeDiversity.clade_id ==
                                    cladeID]['prop_shannon_mean'],
                    color=speciesColorDict[speciesID], \
                        label=' '.join([r" $\sim$","{:.3e}".format(np.round(gRate,12))]),\
                        linewidth=2, linestyle='solid',alpha=0.75)
        #
        axes[5].fill_between(diversityStats['t'],
                            diversityStats['shannon_mean'] -
                            diversityStats['shannon_std'],
                            diversityStats['shannon_mean'] +
                            diversityStats['shannon_std'], color='mediumblue', alpha=0.3)
        axes[5].plot(diversityStats['t'],
                    diversityStats['shannon_mean'],
                    linewidth=2, color='mediumblue')
        axes[5].fill_between(diversityStatsNoV['t'],
                            diversityStatsNoV['shannon_mean'] -
                            diversityStatsNoV['shannon_std'],
                            diversityStatsNoV['shannon_mean'] +
                            diversityStatsNoV['shannon_std'], color='lime', alpha=0.3)
        axes[5].plot(diversityStatsNoV['t'],
                    diversityStatsNoV['shannon_mean'],
                    linewidth=2, color='lime', linestyle='dashed')
        #
    axes[0].set_xscale('log', base=10)
    axes[1].set_xscale('log', base=10)
    axes[2].set_xlabel(xlabel='Time t', fontsize=15, labelpad=15)
    axes[3].set_xlabel(xlabel='Time t', fontsize=15, labelpad=15)
    axes[3].tick_params(axis='x', labelsize=15)
    axes[3].tick_params(axis='y', labelsize=15)
    axes[4].tick_params(axis='x', labelsize=15)
    axes[4].tick_params(axis='y', labelsize=15)
    axes[5].tick_params(axis='x', labelsize=15)
    axes[5].tick_params(axis='y', labelsize=15)
    axes[6].tick_params(axis='x', labelsize=15)
    axes[6].tick_params(axis='y', labelsize=15)
    axes[0].set_xlim(1, 2000)
    axes[1].set_xlim(1, 2000)
    axes[2].set_xlim(1, 2000)
    lim = axes[4].get_xlim()
    axes[4].set_xlim(1, 2000)
    lim = axes[5].get_xlim()
    axes[5].set_xlim(1, 2000)
    lim = axes[6].get_xlim()
    axes[6].set_xlim(1, 2000)
    axes[3].set_xlim(1, 2000)
    axes[4].yaxis.set_label_position("left")
    axes[4].yaxis.tick_left()
    axes[5].yaxis.set_label_position("left")
    axes[5].yaxis.tick_left()
    axes[6].yaxis.set_label_position("left")
    axes[6].yaxis.tick_left()
    s = ['Expected\n{0} Mean '.format(rateName), r'$\mathbb{E}[\bar{r}]$']
    axes[4].set_ylabel(
        ylabel=r''.join(s), labelpad=10, fontsize=15)
    s = ['Expected\n{} SD '.format(rateName), r'$\mathbb{E}[\sigma_f]$']
    axes[3].set_ylabel(ylabel='Fraction of\n{}'.format(diversity), labelpad=10, fontsize=15)
    axes[5].set_ylabel(ylabel = diversity,labelpad=10,fontsize=15)
    axes[6].set_ylabel(ylabel=r''.join(s), labelpad=10, fontsize=15)
    axes[6].set_xlabel(xlabel = 'Time t',fontsize=15,labelpad=15)
    handles = []
    labels = []
    handle, label = axes[0].get_legend_handles_labels()
    handles.extend(handle)
    labels.extend(label)
    handle, label = axes[4].get_legend_handles_labels()
    handles.extend(handle)
    labels.extend(label)
    axes[4].legend(handles, labels, loc='lower left', fontsize=12)
    handles = []
    labels = []
    handle, label = axes[3].get_legend_handles_labels()
    handles.extend(handle[::-1])
    labels.extend(label[::-1])
    axes[3].legend(handles,labels,loc='upper left', fontsize=10, title=propTitle, title_fontsize=10, ncol=2)
    fig.savefig(os.path.join('/Volumes/Yadgah/sylvain-martin-collab',simDir,'diversity-cID{0}.pdf'.format(cID)),dpi=resolve)
    plt.close('all')
##
##
plt.show()










####
### Mean Competitive Ability of Immune Types (NOT Averaged over all simulations; just one)
       # Run code block above first 
runID = 1 #CHOOSE RUN ID       
conSimTree = sqlite3.connect(DBSIM_PATH)
curSimTree = conSimTree.cursor()
bfreqs = pd.read_sql_query("SELECT t, bstrain_id, abundance, run_id FROM babundance \
                            WHERE run_id in ({0})"
                                .format(run_id), conSimTree)\
                .merge(pd.read_sql_query("SELECT t, microbial_abundance, run_id FROM summary \
                                        WHERE run_id in ({0})"
                                            .format(run_id), conSimTree), on=['t','run_id'])\
                .merge(pd.read_sql_query("SELECT bstrain_id, growth_rate, run_id FROM bgrowthrates \
                                        WHERE run_id in ({0})"
                                            .format(run_id), conSimTree), on=['bstrain_id','run_id'])

bfreqs['bfreq'] = bfreqs['abundance']/bfreqs['microbial_abundance']
bfreqs = bfreqs.drop(columns=['microbial_abundance','abundance'])
# bfreqs['growth_rate'] = bfreqs['growth_rate']*bfreqs['bfreq']
bfreqs = bfreqs.merge(immuneIDs,on=['bstrain_id','run_id'])
norm = bfreqs.groupby(['t','immune_id','run_id'])\
                .agg(norm=('bfreq','sum'))\
                    .reset_index()
bfreqs = bfreqs.merge(norm,on=['run_id','t','immune_id'])
bfreqs['expWeight'] = bfreqs['growth_rate']*bfreqs['bfreq']/bfreqs['norm']
norm = bfreqs.groupby(['t','immune_id','run_id'])\
                .agg(exp=('expWeight','sum'))\
                    .reset_index()
bfreqs = bfreqs.merge(norm,on=['run_id','t','immune_id'])
bfreqs['varWeight'] = ((bfreqs['growth_rate'] - bfreqs['exp'])**2)*(bfreqs['bfreq']/bfreqs['norm'])
bfreqs = bfreqs.groupby(['t','immune_id','run_id'])\
                .agg(exp_growth_rate = ('expWeight','sum'), var_growth_rate = ('varWeight','sum'), bfreq=('bfreq','sum'))\
                    .reset_index()
fig, ax = plt.subplots(1)
axes = [ax]
iPalette = 'viridis'
# iCmap = cm.get_cmap(iPalette)
iCmap = plt.colormaps[iPalette]
cNorm = Normalize(vmin=float(1), 
                vmax=float(len(np.unique(immuneIDs[immuneIDs.run_id==runID]['immune_id']))))
colorID = 1
for iID in sorted(bfreqs.immune_id.drop_duplicates().values):
    color = iCmap(cNorm(float(colorID)))
    print(color)
    colorID += 1
    axes[0].fill_between(bfreqs[bfreqs.immune_id == iID]['t'],
                    np.array(bfreqs[bfreqs.immune_id == iID]['exp_growth_rate']) -
                    np.sqrt(bfreqs[bfreqs.immune_id == iID]['var_growth_rate']),
                    np.array(bfreqs[bfreqs.immune_id == iID]['exp_growth_rate']) +
                    np.sqrt(bfreqs[bfreqs.immune_id == iID]['var_growth_rate']), 
                    color=color, alpha=0.25)
    #
colorID = 1
for iID in sorted(bfreqs.immune_id.drop_duplicates().values):
    color = iCmap(cNorm(float(colorID)))
    axes[0].plot(bfreqs[bfreqs.immune_id == iID]['t'],
            bfreqs[bfreqs.immune_id == iID]['exp_growth_rate'],
            color=color, label='{}'.format(colorID), linewidth=2, linestyle='solid',alpha=0.75)
    colorID += 1
#
axes[0].set_ylabel(ylabel ='Competitive Ability',labelpad=10,fontsize=15)
axes[0].set_xlabel(xlabel = 'Time t',fontsize=15,labelpad=15)
handles = []
labels = []
handle, label = axes[0].get_legend_handles_labels()
handles.extend(handle)
labels.extend(label)
axes[0].legend(handles,labels,loc='lower right', fontsize=10, title='Immune ID', title_fontsize=10,ncol=20)
#############
#############
#############






print('SQLite Query: virus abundance data')
virus_total = pd.read_sql_query("SELECT t, viral_abundance \
    FROM summary WHERE run_id = {}".format(run_id), conSim)\
    .rename(columns={'viral_abundance': 'vtotal'})
virus_total = virus_total[virus_total['vtotal'] > 0]
t = max(virus_total['t'].values)
print('SQLite Query: microbe abundance data')
microbe_stacked = pd.read_sql_query("SELECT t,bstrain_id,abundance \
    FROM babundance WHERE run_id = {}".format(run_id), conSim)
microbe_stacked = microbe_stacked[microbe_stacked.t <= t]
microbe_total = pd.read_sql_query("SELECT t, microbial_abundance \
    FROM summary WHERE run_id = {}".format(run_id), conSim)\
    .rename(columns={'microbial_abundance': 'btotal'})
microbe_total = microbe_total[microbe_total['t'] <= t]
maxIDs = microbe_stacked.set_index('t').groupby(['bstrain_id']).agg(t = ('abundance','idxmax'),\
                            maxAbund = ('abundance','max')).reset_index()
maxIDs = microbe_total.merge(maxIDs,on=['t'])
maxIDs['btotal'] = bthreshold*np.array(maxIDs['btotal'])
keepStrains = list(maxIDs[maxIDs['maxAbund']>maxIDs['btotal']]['bstrain_id'].values)
microbe_stacked = microbe_stacked[[(i in keepStrains) for i in microbe_stacked['bstrain_id']]].reset_index(drop=True)
#
virus_stacked = pd.read_sql_query("SELECT t,vstrain_id,abundance \
    FROM vabundance WHERE run_id = {}".format(run_id), conSim)
maxIDs = virus_stacked.set_index('t').groupby(['vstrain_id']).agg(t = ('abundance','idxmax'),\
                            maxAbund = ('abundance','max')).reset_index()
maxIDs = virus_total.merge(maxIDs,on=['t'])
maxIDs['vtotal'] = vthreshold*np.array(maxIDs['vtotal'])
keepStrains = list(maxIDs[maxIDs['maxAbund']>maxIDs['vtotal']]['vstrain_id'].values)
virus_stacked = virus_stacked[[(i in keepStrains) for i in virus_stacked['vstrain_id']]].reset_index(drop=True)
#
microbe_stacked = microbe_stacked[microbe_stacked.t % 1 == 0]
virus_stacked =virus_stacked[virus_stacked.t % 1 == 0]
#
microbe_stacked = microbe_stacked.pivot(index='t',columns='bstrain_id',values='abundance')
virus_stacked = virus_stacked.pivot(index='t',columns='vstrain_id',values='abundance')


print('Compiling viral strain time series plot')
fig, ax = plt.subplots(1)
pal = sns.color_palette("tab20b")
virus_stacked.plot.area(ax=ax, stacked=True, legend=False, linewidth=0,color=pal)
ax.set_ylabel(ylabel ='Viral Strain Abundances',labelpad=15,fontsize=7)
ax.set_xlabel(xlabel = 'Time t',fontsize=7)
ax.ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(25))
lim = ax.get_ylim()
ax.set_ylim(0,lim[1])
plt.tight_layout()
plt.show()

fig, ax = plt.subplots(2,sharex=True)
#fig = plt.figure()
# fig.suptitle('Strain Abundances (run{0}-c{1}-r{2})'.format(run_id,combo_id,replicate))
axes = [ax[0], ax[1]]
microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=pal)
axes[0].set_ylabel(ylabel ='Microbial Immune Abundances N_i',labelpad=15,fontsize=7)
axes[0].set_xlabel(xlabel = 'Time t',fontsize=7)
axes[0].ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(25))
lim = axes[0].get_ylim()
axes[0].set_ylim(0,lim[1])
virus_stacked.plot.area(ax = axes[1],stacked=True,legend=False, linewidth=0,color=pal)
axes[1].set_ylabel(ylabel ='Viral Strain Abundances V_i',labelpad=15,fontsize=7)
axes[1].set_xlabel(xlabel = 'Time t',fontsize=7)
axes[1].ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
axes[1].xaxis.set_minor_locator(ticker.MultipleLocator(25))
lim = axes[1].get_ylim()
axes[1].set_ylim(0,lim[1])
fig.tight_layout()
plt.show()