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
import importlib
from mpl_toolkits.mplot3d import Axes3D
importlib.reload(ff)

#############
if len(sys.argv) == 1:
    simDir = os.getcwd()
    simDirNoV = os.getcwd()
else:
    simDir = sys.argv[1]
    simDirNoV = sys.argv[1]
## 
# simDir = '/Volumes/Yadgah/sylvain-martin-collab/26_MOI3'
# simDirNoV = '/Volumes/Yadgah/sylvain-martin-collab/26_MOI3'
#
DBSIM_PATH = os.path.join(simDir,'sweep_db_gathered.sqlite')
DBSIM_PATHnoV = os.path.join(simDirNoV,'sweep_db_gathered.sqlite')
if not os.path.exists(DBSIM_PATH):
    print("File {} does not exist. Please request from authors or simulate on a cluster.".format(DBSIM_PATH))
    sys.exit(1)
if not os.path.exists(DBSIM_PATHnoV):
    print("File {} does not exist. Please request from authors or simulate on a cluster.".format(DBSIM_PATHnoV))
    sys.exit(1)
conSim = sqlite3.connect(DBSIM_PATH)
curSim = conSim.cursor()
conSimNoV = sqlite3.connect(DBSIM_PATHnoV)
curSimNoV = conSimNoV.cursor()

####
# FIGURE 2 AND SUPP FIGURE 1 (FIGURE 1 CODE BLOCK FOLLOWS THAT OF FIGURE 3 BELOW)
####
numTrees = 1
scale = 3
hosts_per_strain = 100
viruses_per_strain = 100
carryCap = 400000
micMutRep = 0
combos = [(1,0),(2,0)]
vthreshold = 5
bthreshold = 0
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
        subprocess.run(["/trees.jl","{}".format(cID),simDir,"{}".format(bthreshold),"{}".format(vthreshold)]) 
        with open('{0}/runID-cID{1}.txt'.format(simDir,cID)) as f:
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


####
# FIGURE 3
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
propTitle = 'Intrinsic Host Growth Rates (Competitive Abilities)'
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
        vstats,meanDist,meanNoVDist,cladeDiversity,diversityStats,diversityStatsNoV,cID,cIDNoV,scale,\
                speciesColorDict,strainTimes,microbeStrainsDF = \
            ff.make_DFs(i,combos,scale,simDir,DBSIM_PATH,conSim,conSimNoV,figxy,vthreshold,bthreshold,micMutRep,hosts_per_strain, viruses_per_strain,\
                abundthresholdM,abundthresholdV,filterStrains,spacingM,spacingV,maxticksize)
        #
        ff.make_figure(simDir,vstats,meanDist,meanNoVDist,cladeDiversity,diversityStats,diversityStatsNoV,ti,tf,cID,cIDNoV,scale,\
                resolve,speciesColorDict,strainTimes,microbeStrainsDF,rateName,diversity,propTitle,tree_parent_strain_id,tree_strain_id) 


####
# FIGURE 1
####
def fgm(s,x):
    return np.exp(-s*x**2)
####
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
    ax.set_xlabel(xlabel=r'Host Competition Trait $x_i$',labelpad=10,fontsize=20)
    ax.set_ylabel(ylabel=r'Intrinsic Growth Rate $r_i$'+'\n(Competitive Ability)', labelpad=10, fontsize=20)
    ax.tick_params(axis='x', labelsize=20)
    ax.tick_params(axis='y', labelsize=20)
    ax.arrow(0, np.exp(-1/s), -1/s + .075, 0, width = 0.01, ec =None, fc='black')
    ax.text(-1/s*1/2, np.exp(-1/s) - .15, r"$\frac{1}{\sigma}$",fontsize=20, horizontalalignment='center')
    # ax.text(1, 0.75, r"$\sigma = 3$",fontsize=20, horizontalalignment='center')  
    # sign = 1
    for i in range(0,len(alleles)):
        ax.axvline(x=alleles[i], ymax=fgm(s,alleles[i])/ymax, linestyle='dotted', \
                   color=speciesColorDict[strainTimes[strainTimes.tree_parent_bstrain_id==0]['tree_bstrain_id'].values[i]], linewidth=6)
        # sign = -1 * sign
    ax.set_xticks([-1,0,1])
    ax.set_yticks([0, 1])
    ax.set_yticklabels(['m','1+m'])
    fig.tight_layout()
    fig.savefig(os.path.join(simDir,'fgmSigma{0}.pdf'.format(s)),dpi=500,transparent=True)
######
######
s = 3
# Generate data
x = np.linspace(-1, 1, 1000)
y = np.linspace(0, 1, 1000)
X, Y = np.meshgrid(x, y)
Z = np.exp(-s*X**2) + Y**2 - 1# Example function
# Set z-axis range
Z_truncated = np.where((Z >= 0), Z, np.nan)
# Create a figure
fig = plt.figure(figsize=(8, 6),facecolor='none')
ax = fig.add_subplot(111, projection='3d')
# Plot surface
ax.plot_surface(X, Y, Z_truncated, cmap='bone', edgecolor='none',alpha=0.8)
# Labels
ax.set_xlabel(r'Host Competition Trait $x_i$',labelpad=14,fontsize=10)
ax.set_ylabel(r'Number of Viral Matches $\mathcal{M}_i$',labelpad=14,fontsize=10)
ax.set_zlabel(r'Total Fitness $f_i$',labelpad=14,fontsize=10,rotation=180)
# ax.set_title('3D Surface Plot')
ax.set_zlim(0, 1.1)
ax.set_ylim(0, 1)
ax.set_xlim(-1, 1)
ax.set_xticks([-1,  0,  1])  # Set specific x-axis tick locations
ax.set_yticks([0, 0.5, 1])
ax.set_zticks([])
ax.set_yticklabels([0, r'$\frac{N}{2}$', r'$N$'])
ax.set_zticklabels([])
# Add a 2D curve
curve_x = np.linspace(-1, 1, 1000)
curve_y = np.zeros_like(curve_x) + 1  # y = 0 for a 2D curve in the x-z plane
curve_z = np.exp(-s*curve_x**2) + 0.01 # Example curve function
ax.plot(curve_x, curve_y, curve_z, color='black', linewidth=3)
# Remove background
ax.xaxis.pane.fill = False
ax.yaxis.pane.fill = False
ax.zaxis.pane.fill = False
# Remove grid lines
ax.grid(False)
ax.view_init(elev=12, azim=-120)
fig.savefig(os.path.join(simDir,'3d-two-traits.pdf'),dpi=500,transparent=True)



####
# SUPP FIGURE 2
####
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
propTitle = 'Initial Host-intrinsic Growth Rates (Competitive Abilities)'
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
        vstats,meanDist,meanNoVDist,cladeDiversity,diversityStats,diversityStatsNoV,cID,cIDNoV,\
                        micMutSpacer, strainTimes,microbeStrainsDF,speciesColorDict = \
                        ff.make_suppDFs2(i,combos,scale,simDir,DBSIM_PATH,conSim,conSimNoV,figxy,\
                                            vthreshold,bthreshold,micMutRep,hosts_per_strain, viruses_per_strain,\
                                            abundthresholdM,abundthresholdV,filterStrains,spacingM,spacingV,maxticksize)
        # FIGURE
        fig = ff.make_suppfigure2(i,outer_grid,fig,resolve,simDir,scale,vstats,meanDist,meanNoVDist,cladeDiversity,\
                            diversityStats,diversityStatsNoV,ti,tf,cID,cIDNoV,\
                        micMutSpacer, strainTimes,microbeStrainsDF,tree_parent_strain_id,tree_strain_id,speciesColorDict,\
                            rateName,diversity,propTitle)
    fig.savefig(os.path.join(simDir,'supp_diversities-scale{0}.pdf'.format(scale)),dpi=resolve)
    plt.close('all')   
    ##
##
