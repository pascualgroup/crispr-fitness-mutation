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
                            lw=1,
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
        axesTree[i].set_ylabel(ylabel=abundanceTitle, labelpad=15, fontsize=12)
        axesTree[i+2].set_ylabel(ylabel=strainTitle, labelpad=15, fontsize=12)
        axesTree[i+2].set_xlabel(xlabel='Time t', labelpad=10, fontsize=12)
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