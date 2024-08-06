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

resolve = 500
imgTypes = ["pdf"]
treeImgTypes = ["pdf"]
# graphImgTypes = ["pdf"]
figxy = (15, 15)  # setting for tree abundance figure
hratio = [1, 3]  # setting for tree abundance figure
maxticksize = 100  # setting for abundances on individual branches of tree
treepalette = 'turbo'  # Color palette for tree: use contiguous color palette
abundthreshold = 0.01

spacingM, spacingV = 0, 0

DBSIMTree_PATH = os.path.join(
    '/Volumes', 'Yadgah', 'sylvain-martin-collab/12_MOI3/isolates/runID4209-c15-r9/runID4209-c15-r9.sqlite')
DBTREE_PATH = os.path.join(
    '/Volumes', 'Yadgah', 'sylvain-martin-collab/12_MOI3/isolates/runID4209-c15-r9/trees_output.sqlite')
conSimTree = sqlite3.connect(DBSIMTree_PATH)
curSimTree = conSimTree.cursor()
run_id = curSimTree.execute('SELECT DISTINCT run_id FROM summary').fetchall()
run_id = run_id[0][0]

def get_inner_order(children_by_parent,identity):
    children_identities = children_by_parent.get(identity, [])
    if len(children_identities) == 0:
        return [identity, identity]
    inner = [get_inner_order(children_by_parent,c) for c in sorted(children_identities)]
    return [identity] + sum(inner, []) + [identity]

def speciesMullerPlotBoth2(run_id, DBSIMTree_PATH, DBTREE_PATH, maxticksize, abundthreshold, spacingM, spacingV):
    figTree = plt.figure()
    gs = figTree.add_gridspec(7, 1, hspace=0, wspace=.15,
                        height_ratios=[1, 0.1, 1, 0.3, 3, .1, 3])
    ax1, ax6, ax2, ax3, ax4, ax7, ax5 = gs.subplots(sharex='col')
    ax3.remove()
    ax6.remove()
    ax7.remove()
    axesTree = [ax1, ax2, ax4, ax5]
    ##
    conSimTree = sqlite3.connect(DBSIMTree_PATH)
    curSimTree = conSimTree.cursor()
    ID = curSimTree.execute(
        'SELECT combo_id,replicate FROM runs WHERE run_id = {}'.format(run_id)).fetchall()
    combo_id = ID[0][0]
    replicate = ID[0][1]
    RUN_DIR = os.path.join(
        'runID{0}-c{1}-r{2}'.format(run_id, combo_id, replicate))
    for species in ['microbe','virus']:
        if species == 'virus':
            strains = 'vstrains'
            strain = 'vstrain'
            strain_id = 'vstrain_id'
            abundanceTitle = 'Viral\nAbundance'
            strainTitle = "Viral Phylogeny"
            strain_abundance = "viral_abundance"
            straintotal = "vtotal"
            abundance = 'vabundance'
            tree_strain_id = "tree_vstrain_id"
            tree_parent_strain_id = "tree_parent_vstrain_id"
            parent_strain_id = "parent_vstrain_id"
            i = 1
        if species == 'microbe':
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
            treepalette = 'turbo'
            Vcmap = cm.get_cmap(treepalette)
            i = 0
        species_stacked = pd.read_sql_query(
            "SELECT t,{0},abundance FROM {1} WHERE run_id = {2}".format(strain_id, abundance, run_id), conSimTree)
        t = max(pd.read_sql_query(
            "SELECT t FROM vabundance WHERE run_id = {}".format(run_id), conSimTree).t.values)
        species_stacked = species_stacked[species_stacked.t <= t]
        species_total = pd.read_sql_query("SELECT t,{0} FROM summary WHERE run_id = {1}".format(strain_abundance, run_id), conSimTree)\
            .rename(columns={strain_abundance: straintotal})
        species_total = species_total[species_total.t <= max(species_stacked.t)]
        maxIDs = species_stacked.set_index('t').groupby([strain_id]).agg(t=('abundance', 'idxmax'),
                                                                        maxAbund=('abundance', 'max')).reset_index()
        maxIDs = species_total.merge(maxIDs, on=['t'])
        maxIDs[straintotal] = abundthreshold*np.array(maxIDs[straintotal])
        keepStrains = list(
            maxIDs[maxIDs['maxAbund'] > maxIDs[straintotal]][strain_id].values)
        conTree = sqlite3.connect(DBTREE_PATH)
        curTree = conTree.cursor()
        print('SQLite Query: tree data')
        parentTrack = keepStrains.copy()
        keepStrainsNew = keepStrains.copy()
        while len(parentTrack) != 0:
            parentTrack = list(pd.read_sql_query(
                "SELECT parent_{1} \
            FROM {0} WHERE {1} in ({2})".format(strains, strain_id, ', '.join(map(str, parentTrack))), conSimTree)
                [parent_strain_id].values)
            # print(parentTrack)
            keepStrainsNew.extend(parentTrack)
            keepStrainsNew = list(np.unique(keepStrainsNew))
            # print(keepStrainsNew)
            parentTrack = list(np.array(parentTrack)[
                np.array(parentTrack) != 0])
        keepStrainsNew = list(*np.array(keepStrainsNew)[keepStrainsNew != 0])
        keepStrainsNew.sort()
        species_stacked = species_stacked[[
            (i in keepStrainsNew) for i in species_stacked[strain_id]]].reset_index(drop=True)
        keepTreeStrainsDF = pd.read_sql_query(
            "SELECT tree_{0}, {0} \
        FROM tree_{1}_order WHERE {0} in ({2})".format(strain_id, strain, ', '.join(map(str, keepStrainsNew))), conTree)
        keepTreeStrains = list(
            keepTreeStrainsDF['tree_{}'.format(strain_id)].values)
        keepTreeStrains = list(np.array(keepTreeStrains)[
            np.array(keepTreeStrains) != 0])
        keepTreeStrains.sort()

        if species == 'microbe':
            strainTimes = pd.read_sql_query(
                "SELECT tree_{0}, t_creation, t_extinction, tree_parent_{0} \
            FROM tree_{1}_creation_extinction WHERE tree_{0} in ({2})".format(
                    strain_id, strain, ', '.join(map(str, keepTreeStrains))), conTree)
            strainTimesNoExt = pd.read_sql_query(
                "SELECT {0}, t_creation, parent_{0} \
            FROM {1} WHERE {0} in ({2})".format(
                    strain_id, strains, ', '.join(map(str, keepTreeStrainsDF[strain_id]))), conSimTree)\
                .merge(keepTreeStrainsDF, on=strain_id).drop(columns=[strain_id])\
                .rename(columns={'parent_{0}'.format(strain_id): strain_id})\
                .merge(keepTreeStrainsDF.rename(columns={tree_strain_id: 'tree_parent_{0}'.format(strain_id)}), on=strain_id)\
                .drop(columns=[strain_id])
            strainTimesNoExt = strainTimesNoExt[[(i not in list(strainTimes[tree_strain_id]))
                                                for i in list(strainTimesNoExt[tree_strain_id])]]
            strainTimesNoExt['t_extinction'] = len(strainTimesNoExt)*[t]
            strainTimes = pd.concat(
                [strainTimes, strainTimesNoExt]).reset_index(drop=True)
        if species == 'virus':
            strainTimes = pd.read_sql_query(
                "SELECT tree_{0}, t_creation, t_extinction, tree_parent_{0} \
            FROM tree_{1}_creation_extinction WHERE tree_{0} in ({2})".format(
                    strain_id, strain, ', '.join(map(str, keepTreeStrains))), conTree)
            infectionTimes = pd.read_sql_query(
                "SELECT {0}, t_creation, parent_{0}, infected_bstrain_id \
            FROM {1} WHERE {0} in ({2})".format(
                    strain_id, strains, ', '.join(map(str, keepTreeStrainsDF[strain_id]))), conSimTree)\
                .merge(keepTreeStrainsDF, on=strain_id).drop(columns=[strain_id])\
                .rename(columns={'parent_{0}'.format(strain_id): strain_id})\
                .merge(keepTreeStrainsDF.rename(columns={tree_strain_id: 'tree_parent_{0}'.format(strain_id)}), on=strain_id)\
                .drop(columns=[strain_id]).rename(columns={'infected_bstrain_id':'bstrain_id'}).merge(microbeStrainsDF,on='bstrain_id')\
                .drop(columns=['bstrain_id'])
            strainTimesNoExt = infectionTimes.drop(columns=['tree_bstrain_id'])
            strainTimesNoExt = strainTimesNoExt[[(i not in list(strainTimes[tree_strain_id]))
                                                for i in list(strainTimesNoExt[tree_strain_id])]]
            strainTimesNoExt['t_extinction'] = len(strainTimesNoExt)*[t]
            strainTimes = pd.concat(
                [strainTimes, strainTimesNoExt]).reset_index(drop=True)
        #
        if species == 'microbe':
            bstrainRates = pd.read_sql_query("SELECT bstrain_id, locus_allele, growth_rate FROM bgrowthrates", conSimTree)\
                                .merge(keepTreeStrainsDF, on=strain_id).drop(columns=[strain_id])
            bstrainRates['growth_rate'] = np.array(bstrainRates['growth_rate']) - .025
            # bstrainRates['growth_rate'] = np.round(bstrainRates['growth_rate'],3)
            strainTimes = strainTimes.merge(bstrainRates, on=tree_strain_id)
        #
        treeAbundances = pd.read_sql_query(
            "SELECT t, tree_{0}, abundance \
        FROM tree_{1} WHERE tree_{0} in ({2})"
            .format(strain_id, abundance, ', '.join(map(str, keepTreeStrains))), conTree)
        newTreeOrder = {}
        newTreeID = 1
        for treeStrain in keepTreeStrains[::-1]:
            newTreeOrder[treeStrain] = newTreeID
            newTreeID += 1
        if species == 'microbe':
            microbeTreeOrder = newTreeOrder
        treeAbundances.replace({tree_strain_id: newTreeOrder}, inplace=True)
        strainTimes.replace({tree_strain_id: newTreeOrder,
                            tree_parent_strain_id: newTreeOrder}, inplace=True)
        keepTreeStrainsDF['old_{}'.format(tree_strain_id)] = keepTreeStrainsDF[tree_strain_id]
        keepTreeStrainsDF.replace({tree_strain_id: newTreeOrder}, inplace=True)
        if species == 'microbe':
            microbeStrainsDF = keepTreeStrainsDF
            bstrainTimes = strainTimes.copy()
        numStrains = len(keepTreeStrains)
        if species == 'microbe':
            clades = strainTimes[strainTimes[tree_parent_strain_id]==0].sort_values(by=tree_strain_id,ascending=False)
            maxDesc = max(np.append(np.array(list(clades[tree_strain_id])[0:-1]) - np.array(list(clades[tree_strain_id])[1:]),
                        list(clades[tree_strain_id])[-1])) 
            spacing = spacingM
            norm = Normalize(vmin=float(1), vmax=float(len(clades)*(maxDesc+1)+spacing*(len(clades)-1)))
        if species == 'virus':
            clades = strainTimes[strainTimes[tree_parent_strain_id]==0].sort_values(by=tree_strain_id,ascending=False)
            maxDesc = max(np.append(np.array(list(clades[tree_strain_id])[0:-1]) - np.array(list(clades[tree_strain_id])[1:]),
                        list(clades[tree_strain_id])[-1]))
            spacing = spacingV
            norm = Normalize(vmin=float(1), vmax=float(len(clades)*(maxDesc+1)+spacing*(len(clades)-1)))
        speciesColorDict = {}
        maxAbundanceSpecies = np.max(treeAbundances.abundance.values)
        markerIncSpecies = maxticksize/maxAbundanceSpecies
        markerColorsSpecies = []
        hlinecSpecies = []
        vlinecSpecies = []
        hcolorsSpecies = []
        vcolorsSpecies = []
        print('Compiling stacked species abundances and tree plots')
        if species == 'microbe':
            for (l,j,k) in zip(list(clades[tree_strain_id]),np.append(list(clades[tree_strain_id])[1:],0),np.arange(len(clades[tree_strain_id]),0,-1)):
                # print(i,j,keepTreeStrainsDF)
                desc = strainTimes[(strainTimes[tree_strain_id] < l) & (strainTimes[tree_strain_id] > j)]\
                        .sort_values(by=tree_strain_id,ascending=False)[tree_strain_id].values
                color = Vcmap(norm(float(k*(maxDesc+1) + spacing*(k-1))))
                color =  tuple((float(color[0] + 0.2*(1-color[0])), float(color[1] + 0.2*(1-color[1])), float(color[2] + 0.2*(1-color[2])), color[3]))
                # color =  tuple((float(color[0]), float(color[1]), float(color[2] + 0.5*(1-color[2])), color[3]))
                speciesColorDict[l] = color
                if len(desc) == 0:
                    continue
                for a in range(0,len(desc)):
                    color = Vcmap(norm(float(k*(maxDesc+1) + spacing*(k-1) - (a + 1))))
                    color =  tuple((float(color[0] - 0.2*color[0]), float(color[1] - 0.2*color[1]), float(color[2] - 0.2*color[2]), color[3]))
                    print(color)
                    speciesColorDict[desc[a]] = color
            microbeColorDict = speciesColorDict
        if species == 'virus':
            infectionTimes.replace({tree_strain_id: newTreeOrder}, inplace=True)
            infectionTimes.replace({tree_parent_strain_id: newTreeOrder}, inplace=True)
            # infectionTimes.replace({'tree_bstrain_id': microbeTreeOrder}, inplace=True)

            for vstrainID in infectionTimes['tree_vstrain_id']:
                # print(vstrainID)
                if vstrainID not in list(infectionTimes['tree_parent_vstrain_id']):
                    print(vstrainID)
                    origvstrainID = keepTreeStrainsDF[keepTreeStrainsDF[tree_strain_id]==vstrainID]['vstrain_id'].values[0]
                    spacers = [spacer for (spacer,) in curSimTree.execute("SELECT spacer_id FROM vpspacers WHERE vstrain_id = {}"\
                    .format(origvstrainID)).fetchall()]
                    maxAbund = max([abund for (abund,) in curSimTree.execute("SELECT abundance FROM vabundance WHERE vstrain_id = {}"\
                    .format(origvstrainID)).fetchall()])
                    tAbund = max([t for (t,) in curSimTree.execute("SELECT t FROM vabundance WHERE vstrain_id = {} AND abundance = {}"\
                    .format(origvstrainID, maxAbund)).fetchall()])
                    bstrainsT = pd.read_sql_query("SELECT bstrain_id, abundance \
                        FROM babundance WHERE t = {}".format(tAbund-3), conSimTree)
                    iBstrains = [bstrainID for (bstrainID,) in \
                                curSimTree.execute("SELECT bstrain_id FROM bspacers \
                                                    WHERE bstrain_id in ({0}) AND spacer_id in ({1})"\
                    .format(', '.join(map(str, bstrainsT['bstrain_id'])),', '.join(map(str, spacers)))).fetchall()]
                    bsus = list(set(bstrainsT['bstrain_id']) - set(iBstrains))
                    bstrainsT = bstrainsT[bstrainsT['bstrain_id'].isin(bsus)]
                    bstrainsT = bstrainsT[bstrainsT['abundance'] == bstrainsT.max().abundance]
                    bstrainsT = bstrainsT['bstrain_id'].values[0]
                    # maxSusBstrains.append(bstrainsT)
                    speciesColorDict[vstrainID] = \
                        microbeColorDict[microbeStrainsDF[microbeStrainsDF['bstrain_id']==bstrainsT]['tree_bstrain_id'].values[0]]
                    # speciesColorDict[vstrainID] = \
                    #     microbeColorDict[1]
                    continue
                infectID = max(set(infectionTimes[infectionTimes.tree_parent_vstrain_id==vstrainID]['tree_bstrain_id']), \
                               key=list(infectionTimes[infectionTimes.tree_parent_vstrain_id==vstrainID]['tree_bstrain_id']).count) 
                speciesColorDict[vstrainID] = \
                    microbeColorDict[infectID]    
                         
        for strainID in sorted(treeAbundances[tree_strain_id].unique()):
            # print('{0}: {1}'.format(species,strainID))
            tCreate = strainTimes[strainTimes[tree_strain_id]
                                == strainID].t_creation.values[0]
            tExtinct = strainTimes[strainTimes[tree_strain_id]
                                == strainID].t_extinction.values[0]
            parent = strainTimes[strainTimes[tree_strain_id]
                                == strainID][tree_parent_strain_id].values[0]
            hlinecSpecies.append([[tCreate, strainID], [tExtinct, strainID]])
            vlinecSpecies.append([[tCreate, parent], [tCreate, strainID]])  
            hcolorsSpecies.append(speciesColorDict[strainID])
            vcolorsSpecies.append(speciesColorDict[strainID])
            markerColorsSpecies.append(speciesColorDict[strainID])
            axesTree[i+2].scatter(treeAbundances[treeAbundances[tree_strain_id] == strainID]['t'],
                            treeAbundances[treeAbundances[tree_strain_id]
                                        == strainID][tree_strain_id],
                            lw=2,
                            s=(treeAbundances[treeAbundances[tree_strain_id] ==
                                            strainID]['abundance'].values)*markerIncSpecies,
                            color=speciesColorDict[strainID], marker='|')
        strainLineages = LineCollection(
            hlinecSpecies, linestyles='solid', colors=hcolorsSpecies, linewidths=(1))
        creationLines = LineCollection(
            vlinecSpecies, linestyles='solid', colors=vcolorsSpecies, linewidths=(1))
        axesTree[i+2].add_collection(strainLineages)
        axesTree[i+2].add_collection(creationLines)
        axesTree[i+2].set_ylim(0, np.max(strainTimes[tree_strain_id])+1)
        axesTree[i+2].set_xlim(0, t)
        axesTree[i+2].xaxis.set_minor_locator(ticker.MultipleLocator(25))
        axesTree[i+2].set_yticks([])
        muller_total = treeAbundances.groupby(['t']).agg(
            total=('abundance', 'sum')).reset_index()
        mullerAbunds = treeAbundances.merge(muller_total, on=['t'])
        mullerAbunds = mullerAbunds[mullerAbunds.t <= t]
        mullerAbunds['freq'] = mullerAbunds['abundance']/mullerAbunds['total']
        mullerAbunds = mullerAbunds.drop(columns=['abundance', 'total'])
        muller_total['total'] = np.log10(muller_total['total'])
        mullerAbunds = mullerAbunds.merge(muller_total, on=['t'])
        mullerAbunds['freq'] = mullerAbunds['freq']*mullerAbunds['total']
        mullerAbunds = mullerAbunds.pivot(
            index='t', columns=tree_strain_id, values='freq')
        mullerAbunds = mullerAbunds/2
        check = strainTimes.copy().drop(columns=['t_creation', 't_extinction'])\
            .rename(columns={tree_strain_id: 'Identity', tree_parent_strain_id: 'Parent'})
        check = check[check.Parent != 0]
        children_by_parent = check.groupby(
            'Parent')['Identity'].apply(lambda x: list(sorted(x)))
        order = []
        for mrca in sorted(treeAbundances[treeAbundances.t == 0][tree_strain_id].values):
            if mrca not in order:
                order += get_inner_order(children_by_parent, mrca)
        mullerAbunds[order].plot.area(
            ax=axesTree[i], stacked=True, legend=False, linewidth=0, color=speciesColorDict)
        axesTree[i].set_ylabel(ylabel=abundanceTitle, labelpad=15, fontsize=10)
        axesTree[i+2].set_ylabel(ylabel=strainTitle, labelpad=15, fontsize=10)
        axesTree[i+2].set_xlabel(xlabel='Time t', labelpad=10, fontsize=10)
        maxO = int(np.ceil(np.log10(max(species_total[straintotal]))))
        axesTree[i].set_yticks([i for i in range(0, maxO+1, 1)])
        yticklabels = [r'$10^{}$'.format(i)
                        for i in range(0, maxO, 1)]
        yticklabels.append('')
        axesTree[i].set_yticklabels(yticklabels)
        axesTree[i].xaxis.set_minor_locator(ticker.MultipleLocator(25))
        lim = axesTree[i].get_ylim()
        axesTree[i].set_ylim(0, lim[1])
        lim = axesTree[i].get_xlim()
        axesTree[i].set_xlim(lim[0], t)
        axesTree[i+2].tick_params(axis='x', labelsize=12)
        axesTree[i].tick_params(axis='y', labelsize=7)
    # plt.show()    
    return figTree, axesTree, microbeColorDict, bstrainTimes

spacingM, spacingV = 100, 0
_, _, speciesColorDict, strainTimes = speciesMullerPlotBoth2(run_id, DBSIMTree_PATH, DBTREE_PATH, maxticksize, abundthreshold, spacingM, spacingV)
DBSIM_PATH = os.path.join('/Volumes/Yadgah/sylvain-martin-collab/20_MOI3/sweep_db_gathered.sqlite')
DBTRAIT_PATH = os.path.join('/Volumes/Yadgah/sylvain-martin-collab/20_MOI3/sweep_db_gathered.sqlite')
DBSIM_PATHnoV = os.path.join('/Volumes/Yadgah/sylvain-martin-collab/19_MOI3/sweep_db_gathered.sqlite')
DBTRAIT_PATHnoV = os.path.join('/Volumes/Yadgah/sylvain-martin-collab/19_MOI3/sweep_db_gathered.sqlite')
conSim = sqlite3.connect(DBSIM_PATH)
curSim = conSim.cursor()
conTrait = sqlite3.connect(DBTRAIT_PATH)
curTrait = conTrait.cursor()
conSimNoV = sqlite3.connect(DBSIM_PATHnoV)
curSimNoV = conSimNoV.cursor()
conTraitNoV = sqlite3.connect(DBTRAIT_PATHnoV)
curTraitNoV = conTraitNoV.cursor()
micMutProb = 0.001
viralMutProb = 1.04e-6
viralDecay = .1
failureProb = 0
carryCap = 400000
micDeath = 0.025
comboSpace = pd.read_sql_query(
    "SELECT combo_id, evofunctionScale \
FROM param_combos WHERE microbe_mutation_prob_per_replication = {0} AND viral_mutation_rate = {1} \
        AND viral_decay_rate = {2} AND crispr_failure_prob = {3} AND microbe_carrying_capacity = {4} \
        AND microbe_death_rate = {5} \
        ORDER BY combo_id"
    .format(micMutProb, viralMutProb, viralDecay, failureProb, carryCap, micDeath),
    conSim)
comboSpaceNoV = pd.read_sql_query(
    "SELECT combo_id, evofunctionScale \
FROM param_combos WHERE microbe_mutation_prob_per_replication = {0} AND viral_mutation_rate = {1} \
        AND viral_decay_rate = {2} AND crispr_failure_prob = {3} AND microbe_carrying_capacity = {4} \
        AND microbe_death_rate = {5} \
        ORDER BY combo_id"
    .format(micMutProb, viralMutProb, viralDecay, failureProb, carryCap, micDeath),
    conTraitNoV)
if micMutProb == .001:
    cIDNoV = 28 #28 for mutation of .001
    cID = 28 #28 for mutation of .001
if micMutProb == 0:
    cIDNoV = 25 #28 for mutation of .001
    cID = 25 #28 for mutation of .001
# runIDs = pd.read_sql_query(
#     "SELECT combo_id, run_id FROM runs \
#         WHERE combo_id in ({})"
#     .format(', '.join(map(str, comboSpace['combo_id']))), conSim)
runIDs = pd.read_sql_query(
    "SELECT combo_id, run_id FROM runs \
        WHERE combo_id in ({})"
    .format(cID), conSim)
runIDsNoV = pd.read_sql_query(
    "SELECT combo_id, run_id FROM runs \
        WHERE combo_id in ({})"
    .format(cIDNoV), conTraitNoV)
##########
noVDist = pd.DataFrame({'t' : [], 'run_id' : [], 'expectation' : [], 'variance' : []})
for runID in runIDsNoV['run_id']:
    dist = pd.read_sql_query("SELECT t, bstrain_id, abundance, run_id FROM babundance \
                            WHERE run_id in ({0})"
                                .format(runID), conTraitNoV)\
                .merge(pd.read_sql_query("SELECT t, microbial_abundance, run_id FROM summary \
                                        WHERE run_id in ({0})"
                                            .format(runID), conTraitNoV), on=['t'])\
                .merge(pd.read_sql_query("SELECT bstrain_id, growth_rate, run_id FROM bgrowthrates \
                                        WHERE run_id in ({0})"
                                            .format(runID), conTraitNoV), on=['bstrain_id'])
    dist['bfreq'] = dist['abundance']/dist['microbial_abundance']
    dist = dist.drop(columns=['abundance','microbial_abundance'])
    dist['expectation'] = dist['growth_rate']*dist['bfreq']
    dist['variance'] = ((dist['growth_rate'] - dist['expectation'])**2)*dist['bfreq']
    dist = dist.groupby(['t', 'run_id']).agg(expectation=('expectation','sum'), variance=('variance','sum'))\
                .reset_index()
    noVDist = pd.concat([noVDist,dist],ignore_index=True)

noVDist['variance'] = np.array(noVDist['variance']**(1/2))
noVDist = noVDist.rename(columns={'variance': 'std'})
# times = list(range(0,2001,1))
# noVDist = noVDist.sort_values(by=['run_id','t'])
noVDist = noVDist.merge(runIDsNoV, on=['run_id']).merge(comboSpaceNoV, on=['combo_id'])\
    .sort_values(by=['evofunctionScale', 'run_id'])
meanNoVDist = noVDist.groupby(['combo_id', 'evofunctionScale','t'])\
    .agg(meanExp=('expectation', 'mean'), stdExp=('expectation', 'std'), \
         meanStd=('std', 'mean'), stdStd=('std', 'std')).reset_index()
numRuns = noVDist.groupby(['t','combo_id', 'evofunctionScale'])\
    .agg(n=('expectation', 'size')).reset_index()
meanNoVDist = meanNoVDist.merge(numRuns, on=['combo_id', 'evofunctionScale','t'])
########
dist = pd.DataFrame({'t' : [], 'run_id' : [], 'expectation' : [], 'variance' : []})
for runID in runIDs['run_id']:
    subDist = pd.read_sql_query("SELECT t, bstrain_id, abundance, run_id FROM babundance \
                            WHERE run_id in ({0})"
                                .format(runID), conTrait)\
                .merge(pd.read_sql_query("SELECT t, microbial_abundance, run_id FROM summary \
                                        WHERE run_id in ({0})"
                                            .format(runID), conTrait), on=['t'])\
                .merge(pd.read_sql_query("SELECT bstrain_id, growth_rate, run_id FROM bgrowthrates \
                                        WHERE run_id in ({0})"
                                            .format(runID), conTrait), on=['bstrain_id'])
    subDist['bfreq'] = subDist['abundance']/subDist['microbial_abundance']
    subDist = subDist.drop(columns=['abundance','microbial_abundance'])
    subDist['expectation'] = subDist['growth_rate']*subDist['bfreq']
    subDist['variance'] = ((subDist['growth_rate'] - subDist['expectation'])**2)*subDist['bfreq']
    subDist = subDist.groupby(['t', 'run_id']).agg(expectation=('expectation','sum'), variance=('variance','sum'))\
                .reset_index()
    dist = pd.concat([dist,subDist],ignore_index=True)

dist['variance'] = np.array(dist['variance']**(1/2))
dist = dist.rename(columns={'variance': 'std'})
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

immuneIDs = pd.DataFrame({'bstrain_id' : [], 'immune_id' : [], 'run_id' : []})
for runID in runIDs['run_id']:
    immuneID = 1
    arrays = {}
    bspacers = pd.read_sql_query("SELECT bstrain_id, spacer_id FROM bspacers \
                            WHERE run_id in ({0})"
                                .format(runID), conSim)
    
    ## PIVOT AAND THEN DROP DUPLICATES
     #bspacers
     #virus_stacked.pivot(index='t',columns='vstrain_id',values='abundance')
    ##
    for bstrainID in np.unique(bspacers['bstrain_id']):
        spacers = [tuple(bspacers[bspacers.bstrain_id==bstrainID]['spacer_id'].values)]

        if set(spacers).issubset(set(arrays.keys())):
            immuneIDs = pd.concat([immuneIDs, \
                                   pd.DataFrame({'bstrain_id' : [bstrainID], \
                                                 'immune_id' : [arrays[spacers[0]]], \
                                                    'run_id' : [runID]})],\
                                                ignore_index=True)
        else:
            arrays[spacers[0]] = immuneID
            immuneID += 1
            immuneIDs = pd.concat([immuneIDs, \
                        pd.DataFrame({'bstrain_id' : [bstrainID], \
                                        'immune_id' : [arrays[spacers[0]]], \
                                        'run_id' : [runID]})],\
                                    ignore_index=True)
            

fitnessIDs = pd.DataFrame({'bstrain_id' : [], 'fitness_id' : [], 'growth_rate' : [], 'locus_allele' : [], 'run_id' : []})
for runID in runIDs['run_id']:
    immuneID = 1
    arrays = {}
    bspacers = pd.read_sql_query("SELECT bstrain_id, spacer_id, run_id FROM bspacers \
                            WHERE run_id in ({0})"
                                .format(runID), conSim)
    for bstrainID in np.unique(bspacers['bstrain_id']):
        spacers = [tuple(bspacers[bspacers.bstrain_id==bstrainID]['spacer_id'].values)]

        if set(spacers).issubset(set(arrays.keys())):
            immuneIDs = pd.concat([immuneIDs, \
                                   pd.DataFrame({'bstrain_id' : [bstrainID], \
                                                 'immune_id' : [arrays[spacers[0]]], \
                                                    'run_id' : [runID]})],\
                                                ignore_index=True)
        else:
            arrays[spacers[0]] = immuneID
            immuneID += 1
            immuneIDs = pd.concat([immuneIDs, \
                        pd.DataFrame({'bstrain_id' : [bstrainID], \
                                        'immune_id' : [arrays[spacers[0]]], \
                                        'run_id' : [runID]})],\
                                    ignore_index=True)
            





cladeImmune = pd.DataFrame({'t' : [], 'bstrain_id' : [], 'locus_allele' : [], 'growth_rate' : [], 
                            'bfreq' : [], 'shannon' : [], 'simpson' : [], 'run_id' : []})

for runID in runIDs['run_id']:
    subDF = pd.read_sql_query("SELECT bstrain_id, locus_allele, growth_rate, run_id FROM bgrowthrates \
                            WHERE run_id in ({0})"
                                .format(runID), conSim)\
                                .merge(\
                                pd.read_sql_query("SELECT t, bstrain_id, abundance, run_id FROM babundance \
                                                        WHERE run_id in ({0})"
                                                    .format(runID), conSim), on = ['run_id','bstrain_id'] 
                                )\
                                .merge(\
                                pd.read_sql_query("SELECT t, microbial_abundance, run_id FROM summary \
                                                        WHERE run_id in ({0})"
                                                    .format(runID), conSim), on = ['run_id','t'] 
                                )                          
    subDF['bfreq'] = subDF['abundance']/subDF['microbial_abundance']
    subDF = subDF.drop(columns=['abundance','microbial_abundance'])
    subDF['shannon'] = np.exp(-1*np.array(subDF['bfreq'])*np.log(subDF['bfreq']))
    subDF['simpson'] = np.array(subDF['bfreq'])**2
    cladeImmune = pd.concat([cladeImmune,subDF],ignore_index=True)



shannon = cladeImmune.groupby(['t','run_id']).agg(stotal=('shannon','prod')).reset_index()
simpson = cladeImmune.groupby(['t','run_id']).agg(ltotal=('simpson','sum')).reset_index()
simpson['ltotal'] = 1/simpson['ltotal']
richness = cladeImmune[['run_id','t','bstrain_id']].groupby(['t', 'run_id']).agg(
    rtotal=('bstrain_id', 'size')).reset_index()
shannonStats = shannon.groupby(['t']).agg(
    shanMean=('stotal', 'mean'), shanStd=('stotal', 'std')).reset_index()
simpsonStats = simpson.groupby(['t']).agg(
    simpMean=('ltotal', 'mean'), simpStd=('ltotal', 'std')).reset_index()
richnessStats = richness.groupby(['t']).agg(
    richMean=('rtotal', 'mean'), richStd=('rtotal', 'std')).reset_index()




divNoV = pd.DataFrame({'t' : [], 'bstrain_id' : [], 'bfreq' : [], 'shannon' : [], 'simpson' : [], 'run_id' : []})
for runID in runIDsNoV['run_id']:
    subDF = pd.read_sql_query("SELECT t, bstrain_id, abundance, run_id FROM babundance \
                            WHERE run_id in ({0})".format(runID), conSimNoV)\
                            .merge(
                            pd.read_sql_query("SELECT t, microbial_abundance, run_id FROM summary \
                                                    WHERE run_id in ({0})"
                                                .format(runID), conSimNoV), on = ['run_id','t'] 
                            )                          
    subDF['bfreq'] = subDF['abundance']/subDF['microbial_abundance']
    subDF = subDF.drop(columns=['abundance','microbial_abundance'])
    subDF['shannon'] = np.exp(-1*np.array(subDF['bfreq'])*np.log(subDF['bfreq']))
    subDF['simpson'] = np.array(subDF['bfreq'])**2
    divNoV = pd.concat([divNoV ,subDF],ignore_index=True)


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
cladeImmune = cladeImmune.groupby(['run_id','t','locus_allele','growth_rate'])\
    .agg(shannon=('shannon', 'prod'), richness=('bstrain_id', 'size')).reset_index()\
                        .merge(shannon,on=['run_id','t'])\
                        .merge(richness, on=['run_id', 't'])
cladeImmune['shanP'] = cladeImmune['shannon']/cladeImmune['stotal']
cladeImmune['richP'] = cladeImmune['richness']/cladeImmune['rtotal']

cladeIstats = cladeImmune.groupby(['t','locus_allele','growth_rate'])\
                .agg(shanPmean=('shanP','mean'),n=('shanP','size'),
                     shanPstd=('shanP','std'),
                     shanPmedian=('shanP','median'),
                     richPmean=('richP', 'mean'),
                     richPstd=('richP', 'std'),
                     richPmedian=('richP', 'median'))\
                .reset_index()
cladeIstats.loc[:, 'growth_rate'] = cladeIstats['growth_rate'] - .025
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
vstats = vstats[vstats.n >= 50]

####
####
####
####
####

rateName = 'Competitive Ability'
diversity = "richness"
propTitle = 'Host Clades'
if diversity == "richness":
    propDivTitle = 'Fraction of\nImmune Richness'
if diversity == "shannon":
    propDivTitle = 'Fraction of\nShannon Immune Diversity'
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
fig = plt.figure()
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
        print(i)
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
            linewidth=2, color='mediumblue', label=r'$\sigma = 3$')
axes[4].fill_between(meanNoVDist[meanNoVDist['combo_id'] == cIDNoV]['t'],
                    meanNoVDist[meanNoVDist['combo_id'] == cIDNoV]['meanExp'] -
                    meanNoVDist[meanNoVDist['combo_id']
                                == cIDNoV]['stdExp'],
                    meanNoVDist[meanNoVDist['combo_id'] == cIDNoV]['meanExp'] +
                    meanNoVDist[meanNoVDist['combo_id'] == cIDNoV]['stdExp'], color='lime', alpha=0.3)
axes[4].plot(meanNoVDist[meanNoVDist['combo_id'] == cIDNoV]['t'],
            meanNoVDist[meanNoVDist['combo_id'] == cIDNoV]['meanExp'],
            linewidth=2, color='lime', label=r'$\sigma = 3$, no virus', linestyle='dashed')
axes[6].fill_between(meanDist[meanDist['combo_id'] == cID]['t'],
                    meanDist[meanDist['combo_id'] == cID]['meanStd'] -
                    meanDist[meanDist['combo_id'] == cID]['stdStd'],
                    meanDist[meanDist['combo_id'] == cID]['meanStd'] +
                    meanDist[meanDist['combo_id'] == cID]['stdStd'], color='mediumblue', alpha=0.3)
axes[6].plot(meanDist[meanDist['combo_id'] == cID]['t'],
            meanDist[meanDist['combo_id'] == cID]['meanStd'],
            linewidth=2, color='mediumblue', label=r'$\sigma=3$')
axes[6].fill_between(meanNoVDist[meanNoVDist['combo_id'] == cIDNoV]['t'],
                    meanNoVDist[meanNoVDist['combo_id'] == cIDNoV]['meanStd'] -
                    meanNoVDist[meanNoVDist['combo_id']
                                == cIDNoV]['stdStd'],
                    meanNoVDist[meanNoVDist['combo_id'] == cIDNoV]['meanStd'] +
                    meanNoVDist[meanNoVDist['combo_id'] == cIDNoV]['stdStd'], color='lime', alpha=0.3)
axes[6].plot(meanNoVDist[meanNoVDist['combo_id'] == cIDNoV]['t'],
            meanNoVDist[meanNoVDist['combo_id'] == cIDNoV]['meanStd'],
            linewidth=2, color='lime', label=r'$\sigma = 3$, no virus', linestyle='dashed')
##
if diversity == "richness":
    speciesIDX = strainTimes[strainTimes[tree_parent_strain_id] == 0][[tree_strain_id,'locus_allele']].sort_values(by=[tree_strain_id])
    for _, row in cladeIstatsTrunc[['locus_allele','growth_rate']].drop_duplicates().sort_values(by=['growth_rate']).iterrows():
        gRate = row['growth_rate']
        speciesID = speciesIDX[speciesIDX.locus_allele==row['locus_allele']][tree_strain_id].values[0]
        axes[3].fill_between(cladeIstatsTrunc[cladeIstatsTrunc.growth_rate == gRate]['t'],
                        cladeIstatsTrunc[cladeIstatsTrunc.growth_rate == gRate]['richPmean'] -
                        cladeIstatsTrunc[cladeIstatsTrunc.growth_rate ==
                                        gRate]['richPstd'],
                        cladeIstatsTrunc[cladeIstatsTrunc.growth_rate == gRate]['richPmean'] +
                        cladeIstatsTrunc[cladeIstatsTrunc.growth_rate == gRate]['richPstd'], 
                        color=speciesColorDict[speciesID], alpha=0.3)
        # break
    for _, row in cladeIstatsTrunc[['locus_allele','growth_rate']].drop_duplicates().sort_values(by=['growth_rate']).iterrows():
        gRate = row['growth_rate']
        speciesID = speciesIDX[speciesIDX.locus_allele==row['locus_allele']][tree_strain_id].values[0]
        axes[3].plot(cladeIstatsTrunc[cladeIstatsTrunc.growth_rate == gRate]['t'],
                cladeIstatsTrunc[cladeIstatsTrunc.growth_rate ==
                                gRate]['richPmean'],
                color=speciesColorDict[speciesID], label='{}'.format(gRate), linewidth=2, linestyle='solid',alpha=0.75)
        # break



    axes[5].fill_between(richnessStats['t'],
                        richnessStats['richMean'] -
                        richnessStats['richStd'],
                        richnessStats['richMean'] +
                        richnessStats['richStd'], color='mediumblue', alpha=0.3)
    axes[5].plot(richnessStats['t'],
                richnessStats['richMean'],
                linewidth=2, color='mediumblue')
    axes[5].fill_between(divNoV['t'],
                        divNoV['richMean'] -
                        divNoV['richStd'],
                        divNoV['richMean'] +
                        divNoV['richStd'], color='lime', alpha=0.3)
    axes[5].plot(divNoV['t'],
                divNoV['richMean'],
                linewidth=2, color='lime', linestyle='dashed')
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
    axes[3].set_ylabel(ylabel=propDivTitle, labelpad=10, fontsize=15)
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
    s = ['Expected\n{0} Mean '.format(rateName), r'$\mathbb{E}[\bar{f}]$']
    axes[4].set_ylabel(
        ylabel=r''.join(s), labelpad=10, fontsize=15)
    s = ['Expected\n{} SD '.format(rateName), r'$\mathbb{E}[\sigma_f]$']
    axes[5].set_ylabel(ylabel ='Immune Richness',labelpad=10,fontsize=15)
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
    axes[3].legend(handles,labels,loc='upper left', fontsize=10, title=propTitle, title_fontsize=10)
    plt.show()


if diversity == "shannon":
    speciesIDX = strainTimes[strainTimes[tree_parent_strain_id] == 0][[tree_strain_id,'locus_allele']]
    for _, row in cladeIstatsTrunc[['locus_allele','growth_rate']].drop_duplicates().sort_values(by=['growth_rate']).iterrows():
        gRate = row['growth_rate']
        speciesID = speciesIDX[speciesIDX.locus_allele==row['locus_allele']][tree_strain_id].values[0]
        axes[3].fill_between(cladeIstatsTrunc[cladeIstatsTrunc.growth_rate == gRate]['t'],
                        cladeIstatsTrunc[cladeIstatsTrunc.growth_rate == gRate]['shanPmean'] -
                        cladeIstatsTrunc[cladeIstatsTrunc.growth_rate ==
                                        gRate]['shanPstd'],
                        cladeIstatsTrunc[cladeIstatsTrunc.growth_rate == gRate]['shanPmean'] +
                        cladeIstatsTrunc[cladeIstatsTrunc.growth_rate == gRate]['shanPstd'], 
                        color=speciesColorDict[speciesID], alpha=0.3)
        # break
    for _, row in cladeIstatsTrunc[['locus_allele','growth_rate']].drop_duplicates().sort_values(by=['growth_rate']).iterrows():
        gRate = row['growth_rate']
        speciesID = speciesIDX[speciesIDX.locus_allele==row['locus_allele']][tree_strain_id].values[0]
        axes[3].plot(cladeIstatsTrunc[cladeIstatsTrunc.growth_rate == gRate]['t'],
                cladeIstatsTrunc[cladeIstatsTrunc.growth_rate ==
                                gRate]['shanPmean'],
                color=speciesColorDict[speciesID], label='{}'.format(gRate), linewidth=2, linestyle='solid',alpha=0.75)
        # break
    axes[5].fill_between(shannonStats['t'],
                        shannonStats['shanMean'] -
                        shannonStats['shanStd'],
                        shannonStats['shanMean'] +
                        shannonStats['shanStd'], color='mediumblue', alpha=0.3)
    axes[5].plot(shannonStats['t'],
                shannonStats['shanMean'],
                linewidth=2, color='mediumblue')
    axes[5].fill_between(divNoV['t'],
                        divNoV['shanMean'] -
                        divNoV['shanStd'],
                        divNoV['shanMean'] +
                        divNoV['shanStd'], color='lime', alpha=0.3)
    axes[5].plot(divNoV['t'],
                divNoV['shanMean'],
                linewidth=2, color='lime', linestyle='dashed')
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
    axes[3].set_ylabel(ylabel=propDivTitle, labelpad=10, fontsize=15)
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
    s = ['Expected\n{0} Mean '.format(rateName), r'$\mathbb{E}[\bar{f}]$']
    axes[4].set_ylabel(
        ylabel=r''.join(s), labelpad=10, fontsize=15)
    s = ['Expected\n{} SD '.format(rateName), r'$\mathbb{E}[\sigma_f]$']
    axes[5].set_ylabel(ylabel ='Shannon\nImmune Diversity',labelpad=10,fontsize=15)
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
    axes[3].legend(handles,labels,loc='upper left', fontsize=10, title=propTitle, title_fontsize=10)
    plt.show()