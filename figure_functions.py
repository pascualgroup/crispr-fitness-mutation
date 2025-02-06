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

def get_inner_order(children_by_parent,identity):
    children_identities = children_by_parent.get(identity, [])
    if len(children_identities) == 0:
        return [identity, identity]
    inner = [get_inner_order(children_by_parent,c) for c in sorted(children_identities)]
    return [identity] + sum(inner, []) + [identity]

def keepWhichStrains(conSimTree, species_stacked, species_total, strain_id, straintotal, parent_strain_id, abundthreshold):
        maxIDs = species_stacked.set_index('t').groupby([strain_id]).agg(t=('abundance', 'idxmax'),
                                                                        maxFreq=('abundance', 'max')).reset_index()
        maxIDs = species_total.merge(maxIDs, on=['t'])
        maxIDs['maxFreq'] = maxIDs['maxFreq']/maxIDs[straintotal] 
        maxIDs = maxIDs.drop(columns=[straintotal])
        keepStrains = list(
            maxIDs[np.array(maxIDs['maxFreq']) >= abundthreshold][strain_id].values)
        #
        parentTrack = keepStrains.copy()
        keepStrainsNew = keepStrains.copy()
        while len(parentTrack) != 0:
            parentTrack = list(pd.read_sql_query(
                "SELECT parent_{1} \
            FROM {0} WHERE run_id = {3} AND {1} in ({2})".format(strains, strain_id, ', '.join(map(str, parentTrack)),run_id), conSimTree)
                [parent_strain_id].values)
            # print(parentTrack)
            keepStrainsNew.extend(parentTrack)
            keepStrainsNew = list(np.unique(keepStrainsNew))
            # print(keepStrainsNew)
            parentTrack = list(np.array(parentTrack)[
                np.array(parentTrack) != 0])
        return list(np.array(keepStrainsNew)[np.array(keepStrainsNew) != 0])

def speciesMullerTreeHostVirus(justStrains, run_id, DBSIM_PATH, DBTREE_PATH, maxticksize, abundthresholdM, abundthresholdV, filterStrains, spacingM, spacingV, figsize):
    if not justStrains:
        figTree = plt.figure(figsize=figsize)
        gs = figTree.add_gridspec(7, 1, hspace=0, wspace=.15,
                            height_ratios=[1, 0.1, 1, 0.3, 3, .3, 3])
        ax1, ax6, ax2, ax3, ax4, ax7, ax5 = gs.subplots(sharex='col')
        ax3.remove()
        ax6.remove()
        ax7.remove()
        axesTree = [ax1, ax2, ax4, ax5]
        #
        fig, ax = plt.subplots(2,sharex=True)
        axes = [ax[0], ax[1], ax[0].twinx().twiny(), ax[1].twinx().twiny()]
    ##
    conSimTree = sqlite3.connect(DBSIM_PATH)
    curSimTree = conSimTree.cursor()
    # metadata for saving plot file
    # ID = curSimTree.execute(
    #     'SELECT combo_id,replicate FROM runs WHERE run_id = {}'.format(run_id)).fetchall()
    # combo_id = ID[0][0]
    # replicate = ID[0][1]
    # RUN_DIR = os.path.join(
    #     'runID{0}-c{1}-r{2}'.format(run_id, combo_id, replicate))
    for species in ['microbe','virus']:
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
            cmap = cm.get_cmap(treepalette)
            i = 0
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
        species_stacked = pd.read_sql_query(
            "SELECT t,{0},abundance FROM {1} WHERE run_id = {2}".format(strain_id, abundance, run_id), conSimTree)
        t = max(pd.read_sql_query(
            "SELECT t FROM vabundance WHERE run_id = {}".format(run_id), conSimTree).t.values)
        species_stacked = species_stacked[species_stacked.t <= t]
        species_total = pd.read_sql_query("SELECT t,{0} FROM summary WHERE run_id = {1}".format(strain_abundance, run_id), conSimTree)\
            .rename(columns={strain_abundance: straintotal})
        species_total = species_total[species_total.t <= max(species_stacked.t)]
        ##
        # print('SQLite Query: tree data')
        conTree = sqlite3.connect(DBTREE_PATH)
        curTree = conTree.cursor()
        # select strains that will be plotted
        if filterStrains == True:
            if species == 'microbe':
                abundthreshold = abundthresholdM
            if species == 'virus':
                abundthreshold = abundthresholdV
            keepTreeStrains = keepWhichStrains(conSimTree, species_stacked, species_total, 
                                           strain_id, straintotal, parent_strain_id, 
                                           abundthreshold)
            species_stacked = species_stacked[[
            (i in keepTreeStrains) for i in species_stacked[strain_id]]].reset_index(drop=True)
            keepTreeStrainsDF = pd.read_sql_query(
                                    "SELECT tree_{0}, {0} \
                                    FROM tree_{1}_order WHERE {0} in ({2})"\
                                    .format(strain_id, strain, ', '.join(map(str, keepTreeStrains))), 
                                    conTree)
            keepTreeStrains = list(keepTreeStrainsDF['tree_{}'.format(strain_id)].values)
        else:
            keepTreeStrainsDF = pd.read_sql_query(
                                    "SELECT tree_{0}, {0} \
                                    FROM tree_{1}_order \
                                    WHERE {0} != 0"\
                                    .format(strain_id, strain), 
                                    conTree)
        keepTreeStrains = list(keepTreeStrainsDF['tree_{}'.format(strain_id)].values)
        numStrains = len(keepTreeStrains)
        keepTreeStrains.sort()      
        #
        if species == 'microbe':
            # compiles creation and extinction times of strains. 
            # also makes "extinction" times of strains that don't go extinct as final simulation time
            strainTimes = pd.read_sql_query(
                "SELECT tree_{0}, t_creation, t_extinction, tree_parent_{0} \
            FROM tree_{1}_creation_extinction WHERE tree_{0} in ({2})".format(
                    strain_id, strain, ', '.join(map(str, keepTreeStrains))), conTree)
            strainTimesNoExt = pd.read_sql_query(
                "SELECT {0}, t_creation, parent_{0} \
            FROM {1} WHERE run_id = {3} AND {0} in ({2})".format(
                    strain_id, strains, ', '.join(map(str, keepTreeStrainsDF[strain_id])),run_id), conSimTree)\
                .merge(keepTreeStrainsDF, on=strain_id).drop(columns=[strain_id])\
                .rename(columns={'parent_{0}'.format(strain_id): strain_id})\
                .merge(keepTreeStrainsDF.rename(columns={tree_strain_id: 'tree_parent_{0}'.format(strain_id)}), on=strain_id)\
                .drop(columns=[strain_id])
            strainTimesNoExt = strainTimesNoExt[[(i not in list(strainTimes[tree_strain_id]))
                                                for i in list(strainTimesNoExt[tree_strain_id])]]
            strainTimesNoExt['t_extinction'] = len(strainTimesNoExt)*[t]
            strainTimes = pd.concat(
                [strainTimes, strainTimesNoExt]).reset_index(drop=True)
            bstrainRates = pd.read_sql_query("SELECT bstrain_id, growth_rate FROM bgrowthrates WHERE run_id = {}".format(run_id), conSimTree)\
                    .merge(keepTreeStrainsDF, on=strain_id).drop(columns=[strain_id])
            bstrainRates['growth_rate'] = np.array(bstrainRates['growth_rate'])
            # bstrainRates['growth_rate'] = np.round(bstrainRates['growth_rate'],3)
            strainTimes = strainTimes.merge(bstrainRates, on=tree_strain_id)
        if species == 'virus':
            strainTimes = pd.read_sql_query(
                "SELECT tree_{0}, t_creation, t_extinction, tree_parent_{0} \
                FROM tree_{1}_creation_extinction WHERE tree_{0} in ({2})"\
                    .format(strain_id, strain, ', '.join(map(str, keepTreeStrains))), conTree)
            strainTimesNoExt = pd.read_sql_query(
                                "SELECT {0}, t_creation, parent_{0} \
                                FROM {1} WHERE run_id = {3} AND parent_{0} in ({2})"\
                                .format(strain_id, strains, ', '.join(map(str, keepTreeStrainsDF[strain_id])),run_id), 
                                        conSimTree)\
                                .merge(keepTreeStrainsDF, on=strain_id)\
                                .drop(columns=['vstrain_id'])\
                                .rename(columns={'parent_vstrain_id':'vstrain_id'})\
                                .merge(keepTreeStrainsDF.rename(columns={'tree_vstrain_id':'tree_parent_vstrain_id'}), on=strain_id)\
                                .drop(columns=['vstrain_id'])
            strainTimesNoExt = strainTimesNoExt[[(i not in list(strainTimes[tree_strain_id]))
                                                for i in list(strainTimesNoExt[tree_strain_id])]]
            strainTimesNoExt['t_extinction'] = len(strainTimesNoExt)*[t]
            strainTimes = pd.concat(
                [strainTimes, strainTimesNoExt]).reset_index(drop=True)
            infectionTimes = pd.read_sql_query(
                "SELECT {0}, t_creation, parent_{0}, infected_bstrain_id \
                FROM {1} WHERE run_id = {3} AND parent_{0} in ({2})"\
                .format(strain_id, strains, ', '.join(map(str, keepTreeStrainsDF[strain_id])),run_id), 
                        conSimTree)\
                .drop(columns=[strain_id])\
                .rename(columns={'parent_{0}'.format(strain_id): strain_id})\
                .merge(keepTreeStrainsDF, on=strain_id)\
                .drop(columns=[strain_id]).rename(columns={'infected_bstrain_id':'bstrain_id'})\
                .merge(microbeStrainsDF.drop(columns=['old_tree_bstrain_id']),on='bstrain_id')\
                .drop(columns=['bstrain_id'])\
                    .rename(columns={'tree_{0}'.format(strain_id): 'tree_parent_vstrain_id'})         
        # call abundances with tree strain IDs
        treeAbundances = pd.read_sql_query(
                            "SELECT t, tree_{0}, abundance \
                             FROM tree_{1} WHERE tree_{0} in ({2})"
                            .format(strain_id, abundance, ', '.join(map(str, keepTreeStrains))), 
                            conTree)
        # reorder strain IDs and relabel DFs
        newTreeOrder = {}
        newTreeID = 1
        for treeStrain in keepTreeStrains[::-1]:
            newTreeOrder[treeStrain] = newTreeID
            newTreeID += 1
        strainTimes.replace({tree_strain_id: newTreeOrder,
                            tree_parent_strain_id: newTreeOrder}, 
                            inplace=True)
        treeAbundances.replace({tree_strain_id: newTreeOrder}, 
                               inplace=True)
        keepTreeStrainsDF['old_{}'.format(tree_strain_id)] = keepTreeStrainsDF[tree_strain_id]
        keepTreeStrainsDF.replace({tree_strain_id: newTreeOrder}, 
                                  inplace=True)
        if species == 'virus':
            infectionTimes.replace({tree_parent_strain_id: newTreeOrder}, 
                                inplace=True)    
        # determine color spacing of strains
        if species == 'microbe':
            ## save for viral tree
            microbeStrainsDF = keepTreeStrainsDF 
            bstrainTimes = strainTimes.copy()
            bMRCAs = sorted(list(strainTimes[strainTimes['tree_parent_bstrain_id']==0]['tree_bstrain_id']))
            ##
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
                color = cmap(norm(float(k*(maxDesc+1) + spacing*(k-1))))
                color =  tuple((float(color[0] + 0.2*(1-color[0])), float(color[1] + 0.2*(1-color[1])), float(color[2] + 0.2*(1-color[2])), color[3]))
                # color =  tuple((float(color[0]), float(color[1]), float(color[2] + 0.5*(1-color[2])), color[3]))
                speciesColorDict[l] = color
                if len(desc) == 0:
                    continue
                for a in range(0,len(desc)):
                    color = cmap(norm(float(k*(maxDesc+1) + spacing*(k-1) - (a + 1))))
                    color =  tuple((float(color[0] - 0.2*color[0]), float(color[1] - 0.2*color[1]), float(color[2] - 0.2*color[2]), color[3]))
                    # print(color)
                    speciesColorDict[desc[a]] = color
            microbeColorDict = speciesColorDict
        if species == 'virus':
            for vstrainID in keepTreeStrainsDF['tree_vstrain_id']:
                # print(vstrainID)
                if vstrainID not in list(infectionTimes['tree_parent_vstrain_id']):
                    # print(vstrainID)
                    origvstrainID = keepTreeStrainsDF[keepTreeStrainsDF[tree_strain_id]==vstrainID]['vstrain_id'].values[0]
                    spacers = [spacer for (spacer,) in curSimTree.execute("SELECT spacer_id FROM vpspacers WHERE run_id = {0} AND vstrain_id = {1}"\
                    .format(run_id, origvstrainID)).fetchall()]
                    # maxAbund = max([abund for (abund,) in curSimTree.execute("SELECT abundance FROM vabundance WHERE run_id = {0} AND vstrain_id = {1}"\
                    # .format(run_id, origvstrainID)).fetchall()])
                    # tAbund = max([t for (t,) in curSimTree.execute("SELECT t FROM vabundance WHERE run_id = {0} AND vstrain_id = {1} AND abundance = {2}"\
                    # .format(run_id, origvstrainID, maxAbund)).fetchall()])
                    times = [t for (t,) in curSimTree.execute("SELECT t FROM vabundance WHERE run_id = {0} AND vstrain_id = {1}"\
                    .format(run_id, origvstrainID)).fetchall()]
                    times = [min(times)]
                    # bstrainsT = pd.read_sql_query("SELECT bstrain_id, abundance \
                    #     FROM babundance WHERE run_id = {0} AND t = {1}".format(run_id, tAbund-1), conSimTree)
                    # iBstrains = [bstrainID for (bstrainID,) in \
                    #             curSimTree.execute("SELECT bstrain_id FROM bspacers \
                    #                                 WHERE run_id = {0} AND bstrain_id in ({1}) AND spacer_id in ({2})"\
                    # .format(run_id, ', '.join(map(str, bstrainsT['bstrain_id'])),', '.join(map(str, spacers)))).fetchall()]
                    # bsus = list(set(bstrainsT['bstrain_id']) - set(iBstrains))
                    bstrainsTimes = pd.read_sql_query("SELECT bstrain_id, abundance \
                        FROM babundance WHERE run_id = {0} AND t in ({1})".format(run_id, ', '.join(map(str, times))), conSimTree)
                    iBstrains = [bstrainID for (bstrainID,) in \
                                curSimTree.execute("SELECT bstrain_id FROM bspacers \
                                                    WHERE run_id = {0} AND bstrain_id in ({1}) AND spacer_id in ({2})"\
                    .format(run_id, ', '.join(map(str, bstrainsTimes['bstrain_id'])),', '.join(map(str, spacers)))).fetchall()]
                    bsus = list(set(bstrainsTimes['bstrain_id']) - set(iBstrains))
                    bsus = list(microbeStrainsDF[microbeStrainsDF['bstrain_id'].isin(bsus)]['bstrain_id'])
                    # print('this is msDF {}'.format(microbeStrainsDF[microbeStrainsDF['bstrain_id'].isin(bsus)]['bstrain_id']))
                    # print('this is bsus {}'.format(bsus))
                    bstrainsT = bstrainsTimes[bstrainsTimes['bstrain_id'].isin(bsus)]
                    bstrainsT = bstrainsT[bstrainsT['abundance'] == bstrainsT.max().abundance]
                    bstrainsT = list(bstrainsT['bstrain_id'].values)
                    bstrainsT = bstrainsT[0]
                    # maxSusBstrains.append(bstrainsT)
                    # print('bstrainsT is {}'.format(bstrainsT))
                    # print('this is in question: {}'.format(microbeStrainsDF[microbeStrainsDF['bstrain_id']==bstrainsT]['tree_bstrain_id'].values))
                    speciesColorDict[vstrainID] = \
                        microbeColorDict[microbeStrainsDF[microbeStrainsDF['bstrain_id']==bstrainsT]['tree_bstrain_id'].values[0]]
                    # speciesColorDict[vstrainID] = \
                    #     microbeColorDict[1]
                    continue
                infectID = max(set(infectionTimes[infectionTimes.tree_parent_vstrain_id==vstrainID]['tree_bstrain_id']), \
                               key=list(infectionTimes[infectionTimes.tree_parent_vstrain_id==vstrainID]['tree_bstrain_id']).count) 
                speciesColorDict[vstrainID] = \
                    microbeColorDict[infectID]    
        if not justStrains:                 
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
                if species == 'virus':
                    axesTree[i+2].scatter(treeAbundances[treeAbundances[tree_strain_id] == strainID]['t'],
                                    treeAbundances[treeAbundances[tree_strain_id]
                                                == strainID][tree_strain_id],
                                    lw=2,
                                    s=(treeAbundances[treeAbundances[tree_strain_id] ==
                                                    strainID]['abundance'].values)*markerIncSpecies,
                                    color=speciesColorDict[strainID], marker='|')
                if species == 'microbe':
                    # if strainID in bMRCAs:
                    #     gRate = np.round(strainTimes[strainTimes['tree_bstrain_id']==strainID]['growth_rate'][0],4)
                    #     print(gRate)
                    #     axesTree[i+2].scatter(treeAbundances[treeAbundances[tree_strain_id] == strainID]['t'],
                    #                     treeAbundances[treeAbundances[tree_strain_id]
                    #                                 == strainID][tree_strain_id],
                    #                     lw=2,
                    #                     s=(treeAbundances[treeAbundances[tree_strain_id] ==
                    #                                     strainID]['abundance'].values)*markerIncSpecies,
                    #                     color=speciesColorDict[strainID], marker='|',label= "{:e}".format(gRate))
                    # else:
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
            muller_total['linear_total'] = muller_total['total'].copy()
            muller_total['total'] = np.log10(muller_total['total'])
            mullerAbunds = mullerAbunds.merge(muller_total, on=['t'])
            mullerAbunds['freq'] = mullerAbunds['freq']*mullerAbunds['total']
            mullerAbunds['linear_freq'] = mullerAbunds['freq']*mullerAbunds['linear_total']
            mullerLinearAbunds = mullerAbunds.pivot(
                index='t', columns=tree_strain_id, values='linear_freq')
            mullerAbunds = mullerAbunds.pivot(
                index='t', columns=tree_strain_id, values='freq')
            mullerLinearAbunds = mullerLinearAbunds/2
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
            if species == 'microbe':
                axes[i].plot(species_total['t'], species_total['btotal'],linewidth=0,color='grey')
                axes[i].fill_between(species_total['t'],species_total['btotal'], color='grey',alpha=0.3)
                mullerLinearAbunds[order].plot.area(
                ax=axes[i+2], stacked=True, legend=False, linewidth=0, color=speciesColorDict)
                lim = axes[i].get_ylim()
                axes[i].set_ylim(0,lim[1])
                lim = axes[i+2].get_ylim()
                axes[i+2].set_ylim(0,lim[1])
                axes[i+2].set_xlabel(xlabel='')
                axes[i+2].set_yticklabels([])
                axes[i+2].set_yticks([])
                axes[i+2].set_xticklabels([])
                axes[i+2].set_xticks([])
                axes[i].margins(x=0)
                axes[i+2].margins(x=0)
                axes[i].ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
                axes[i].yaxis.get_offset_text().set_fontsize(12)
                axes[i].tick_params(axis='x', labelsize=12)
                axes[i].tick_params(axis='y', labelsize=12)
                #
                mullerAbunds[order].plot.area(
                    ax=axesTree[i], stacked=True, legend=False, linewidth=0, color=speciesColorDict)
                handles, label = axesTree[i].get_legend_handles_labels()
                labels = [int(k) for k in label]
                legendLabelsDF = pd.DataFrame({'handle': handles,'tree_bstrain_id':labels})
                legendLabelsDF = legendLabelsDF[legendLabelsDF['tree_bstrain_id'].isin(bMRCAs)]
                handles = []
                legendLabels = []
                for k in bMRCAs[::-1]:
                    # print(k)
                    handles.extend(list([legendLabelsDF[legendLabelsDF['tree_bstrain_id']==k]['handle'].values[0]]))
                    legendLabels.extend([' '.join([r" $\sim$","{}".format(np.round(\
                                    strainTimes[strainTimes['tree_bstrain_id']==k]['growth_rate'].values[0],3))])])
                # handles = np.array(handle)[list(range(0,2*len(bMRCAs),2))][::-1]
                # labels = np.array(label)[list(range(0,2*len(bMRCAs),2))][::-1]
                # legendLabels = ["{:e}".format(np.round(\
                #                     strainTimes[strainTimes['tree_bstrain_id']==int(strain)]['growth_rate'].values[0],4)
                #                     ) for strain in labels]
                print(legendLabels)
                legend = axesTree[2].legend(handles, legendLabels, loc='lower right', title="Intrinsic Host\nGrowth Rates", fontsize=10, title_fontsize=10)
                legend.get_frame().set_alpha(0.8)
            if species == 'virus':
                axes[i].plot(species_total['t'], species_total['vtotal'],linewidth=0,color='grey')
                axes[i].fill_between(species_total['t'],species_total['vtotal'], color='grey',alpha=0.3)
                mullerLinearAbunds[order].plot.area(
                ax=axes[i+2], stacked=True, legend=False, linewidth=0, color=speciesColorDict)
                lim = axes[i].get_ylim()
                axes[i].set_ylim(0,lim[1])
                lim = axes[i+2].get_ylim()
                axes[i+2].set_ylim(0,lim[1])
                axes[i+2].set_xlabel(xlabel='')
                axes[i+2].set_yticklabels([])
                axes[i+2].set_yticks([])
                axes[i+2].set_xticklabels([])
                axes[i+2].set_xticks([])
                axes[i].margins(x=0)
                axes[i+2].margins(x=0)
                axes[i].ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
                axes[i].yaxis.get_offset_text().set_fontsize(12)
                axes[i].tick_params(axis='x', labelsize=12)
                axes[i].tick_params(axis='y', labelsize=12)
                #
                mullerAbunds[order].plot.area(
                    ax=axesTree[i], stacked=True, legend=False, linewidth=0, color=speciesColorDict)
            #
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
            axesTree[i].tick_params(axis='y', labelsize=10)
            #
            axes[i].set_ylabel(ylabel=abundanceTitle, labelpad=15, fontsize=12)
            axes[i].set_xlabel(xlabel='Time t', labelpad=10, fontsize=12)
    # keepTreeStrainsDF is viral tree IDs    
    if not justStrains:        
        return figTree, axesTree, fig, axes, microbeColorDict, bstrainTimes, keepTreeStrainsDF, microbeStrainsDF
    if justStrains:
        return [], [], [], [], microbeColorDict, bstrainTimes, keepTreeStrainsDF, microbeStrainsDF

def func(i):
    return not math.isnan(i)

def findImmuneIDs(runIDs,conSim):
    print('...finding immunes IDs...')
    immuneIDs = pd.DataFrame({'bstrain_id' : [], 'immune_id' : [], 'run_id' : []})
    for runID in runIDs['run_id']:
        # print(runID)
        bstrains = pd.read_sql_query("SELECT bstrain_id, parent_bstrain_id, \
                                    infecting_vstrain_id, run_id FROM bstrains \
                                    WHERE run_id in ({0})"
                                    .format(runID), conSim)
        immuneBstrains = bstrains[list(map(lambda i:func(i), bstrains.infecting_vstrain_id))]
        bstrains = bstrains[list(map(math.isnan, bstrains.infecting_vstrain_id))]
        # not looking for repetitions because 
        ## we assumee unique time of spacer acqusitions will make genotype unique despite overlap
        # spacers = [[sID for (sID,) in 
        #     curSim.execute(
        #         "SELECT spacer_id FROM bspacers \
        #             WHERE run_id = {0} AND bstrain_id = {1}"\
        #             .format(runID,bstrainID)).fetchall()] for bstrainID in immuneBstrains['bstrain_id'][immuneBstrains.bstrain_id!=1]]
        # spacers.count()
        for immuneID in immuneBstrains['bstrain_id']:
            numDescendents = 1
            lineage = [immuneID]
            ancestors = [immuneID]
            descendents = []
            while numDescendents > 0: 
                for id in ancestors:
                    desc_ids = list(bstrains[bstrains.parent_bstrain_id == id]['bstrain_id'].values)
                    # print(desc_ids)
                    # desc_ids = [id for (id,) in execute(dbTemp,"SELECT bstrain_id FROM bstrains \
                    #         WHERE parent_bstrain_id in ({}) ORDER BY bstrain_id".format(id)).fetchall()]
                    descendents.extend(desc_ids)
                    lineage.extend(desc_ids)
                    lineage = list(np.unique(lineage))
                numDescendents = len(descendents)
                ancestors = descendents
                descendents = []
            # print(pd.DataFrame({'bstrain_id' : lineage, 'immune_id' : [immuneID]*len(lineage), 'run_id' : [runID]*len(lineage)}))    
            immuneIDs = pd.concat([immuneIDs, 
                    pd.DataFrame({'bstrain_id' : list(map(int,lineage)), 
                                    'immune_id' : [int(immuneID)]*len(lineage), 
                                    'run_id' : [int(runID)]*len(lineage)})],
                        ignore_index=True)
    return immuneIDs 

def findCladeIDs(runIDs,conSim):
    print('...finding clade IDs...')
    cladeIDs = pd.DataFrame({'bstrain_id' : [], 'clade_id' : [], 'run_id' : []})
    for runID in runIDs['run_id']:
        # print(runID)
        bstrains = pd.read_sql_query("SELECT bstrain_id, parent_bstrain_id, \
                                    infecting_vstrain_id, run_id FROM bstrains \
                                    WHERE run_id in ({0})"
                                    .format(runID), conSim)
        cladeStrains = bstrains[bstrains.parent_bstrain_id == 0]
        for cladeID in cladeStrains['bstrain_id']:
            numDescendents = 1
            lineage = [cladeID]
            ancestors = [cladeID]
            descendents = []
            while numDescendents > 0: 
                for id in ancestors:
                    desc_ids = list(bstrains[bstrains.parent_bstrain_id == id]['bstrain_id'].values)
                    # print(desc_ids)
                    # desc_ids = [id for (id,) in execute(dbTemp,"SELECT bstrain_id FROM bstrains \
                    #         WHERE parent_bstrain_id in ({}) ORDER BY bstrain_id".format(id)).fetchall()]
                    descendents.extend(desc_ids)
                    lineage.extend(desc_ids)
                    lineage = list(np.unique(lineage))
                numDescendents = len(descendents)
                ancestors = descendents
                descendents = []
            cladeIDs = pd.concat([cladeIDs, 
                    pd.DataFrame({'bstrain_id' : list(map(int,lineage)), 
                                    'clade_id' : [int(cladeID)]*len(lineage), 
                                    'run_id' : [int(runID)]*len(lineage)})],
                        ignore_index=True) 
    return cladeIDs

def computeTraitDistribution(runIDs,cID,conSim):
    print('...computing trait distributions...')
    comboSpace = pd.read_sql_query(
    "SELECT combo_id, evofunctionScale \
        FROM param_combos WHERE combo_id = {}"
    .format(cID),conSim)
    dist = pd.DataFrame({'t' : [], 'run_id' : [], 'expectation' : [], 'variance' : []})
    for runID in runIDs['run_id']:
        subDist = pd.read_sql_query("SELECT t, bstrain_id, abundance, run_id FROM babundance \
                                WHERE run_id in ({0})"
                                    .format(runID), conSim)\
                    .merge(pd.read_sql_query("SELECT t, microbial_abundance, run_id FROM summary \
                                            WHERE run_id in ({0})"
                                                .format(runID), conSim), on=['t', 'run_id'])\
                    .merge(pd.read_sql_query("SELECT bstrain_id, growth_rate, run_id FROM bgrowthrates \
                                            WHERE run_id in ({0})"
                                                .format(runID), conSim), on=['bstrain_id','run_id'])
        subDist['bfreq'] = subDist['abundance']/subDist['microbial_abundance']
        subDist = subDist.drop(columns=['abundance','microbial_abundance'])
        subDist['expectation'] = subDist['growth_rate']*subDist['bfreq']
        subDist['variance'] = ((subDist['growth_rate'] - subDist['expectation'])**2)*subDist['bfreq']
        subDist = subDist.groupby(['t', 'run_id']).agg(expectation=('expectation','sum'), variance=('variance','sum'))\
                    .reset_index()
        dist = pd.concat([dist,subDist],ignore_index=True)
    ###
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
    return dist, meanDist

def computeSpacerDistribution(runIDs,immuneIDs,conSim):
    print('...computing spacer statistics...')
    diversityStats = pd.DataFrame({'run_id' : [], 't' : [],
                                'richness' : [],
                                'shannon' : [],
                                'simpson' : []})
    for runID in runIDs['run_id']:
        # print(runID)
        bfreqs = pd.read_sql_query("SELECT t, bstrain_id, abundance, run_id FROM babundance \
                                    WHERE run_id in ({0})"
                                        .format(runID), conSim)\
                        .merge(pd.read_sql_query("SELECT t, microbial_abundance, run_id FROM summary \
                                                WHERE run_id in ({0})"
                                                    .format(runID), conSim), on=['t','run_id'])
        #
        bfreqs['bfreq'] = bfreqs['abundance']/bfreqs['microbial_abundance']
        bfreqs = bfreqs.drop(columns=['microbial_abundance','abundance'])
        bfreqs = bfreqs.merge(immuneIDs,on=['bstrain_id','run_id'])\
                    .groupby(['run_id','t','immune_id']).agg(bfreq=('bfreq','sum'))\
                    .reset_index().rename(columns={'immune_id':'bstrain_id'})
        bfreqs['shannon']  = np.exp(-1*np.array(bfreqs['bfreq'])*np.log(bfreqs['bfreq']))
        bfreqs['simpson'] = np.array(bfreqs['bfreq'])**2
        # bfreqs['simpson'] = 1/bfreqs['simpson']
        diversityStats = pd.concat([diversityStats , bfreqs.groupby(['t','run_id'])\
            .agg(richness=('bfreq','size'), 
                shannon=('shannon','prod'), 
                simpson=('simpson','sum')).reset_index()], ignore_index=True)
    #
    diversityStats['simpson'] = 1/np.array(diversityStats['simpson'])
    diversityStats = diversityStats.groupby(['t'])\
        .agg(richness_mean=('richness','mean'),
                richness_std=('richness','std'), 
                shannon_mean=('shannon','mean'),
                shannon_std=('shannon','std'), 
                simpson_mean=('simpson','mean'),
                simpson_std=('simpson','std')).reset_index()
    return diversityStats

def computeCladeDiversity(runIDs,immuneIDs,cladeIDs,conSim):
    print('...computing clade diversity...')
    cladeDiversity = pd.DataFrame({'t' : [], 'clade_id' : [], 
                                'prop_richness' : [], 
                                'prop_shannon' : [], 
                                #'prop_simpson' : [],
                                'clade_richness' : [],
                                'clade_shannon' : [],
                                'clade_simpson' : [],
                                'run_id' : []})
    #
    for runID in runIDs['run_id']:
        # print(runID)
        subDF = cladeIDs.merge(\
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
        subDF = subDF.merge(immuneIDs,on=['bstrain_id','run_id'])\
                    .groupby(['run_id','t','immune_id','clade_id']).agg(bfreq=('bfreq','sum'))\
                    .reset_index().rename(columns={'immune_id':'bstrain_id'})
        subDF['shannon'] = np.exp(-1*np.array(subDF['bfreq'])*np.log(subDF['bfreq']))
        subDF['simpson'] = np.array(subDF['bfreq'])**2
        # subDF['simpson'] = 1/subDF['simpson']
        subDF = subDF.groupby(['run_id','t','clade_id']).agg(norm=('bfreq','sum')).reset_index()\
                .merge(subDF,on=['run_id','t','clade_id'])
        subDF = subDF[['run_id','t','bstrain_id','bfreq','shannon','simpson']].drop_duplicates().groupby(['run_id','t'])\
            .agg(richness=('bfreq','size'),
                shannon=('shannon','prod'),
                simpson=('simpson','sum')).reset_index()\
            .merge(subDF.drop(columns=['shannon','simpson']),on=['run_id','t'])
        subDF['prop_shannon'] = np.exp(-1*np.array(subDF['bfreq'])*np.log(subDF['bfreq']))
        # subDF['prop_simpson'] = np.array(subDF['bfreq'])**2
        subDF['clade_shannon'] = np.exp(-1*np.array(subDF['bfreq']/subDF['norm'])*np.log(subDF['bfreq']/subDF['norm']))
        subDF['clade_simpson'] = np.array(subDF['bfreq']/subDF['norm'])**2
        subDF = subDF.groupby(['t','clade_id', 'richness','shannon','simpson', 'run_id'])\
                                .agg(prop_richness=('bfreq','size'), 
                                    prop_shannon=('prop_shannon','prod'),
                                    # prop_simpson=('prop_simpson','sum'),
                                    clade_richness=('bfreq','size'),
                                    clade_shannon=('clade_shannon','prod'),
                                    clade_simpson=('clade_simpson','sum')).reset_index()
        subDF['prop_richness'] = subDF['prop_richness']/subDF['richness']
        subDF['prop_shannon'] = subDF['prop_shannon']/subDF['shannon']
        # subDF['prop_simpson'] = subDF['prop_simpson']/subDF['prop_simpson']
        cladeDiversity = pd.concat([cladeDiversity, subDF.drop(columns=['richness','shannon','simpson'])],ignore_index=True)
    ###
    ###
    cladeDiversity = cladeDiversity.groupby(['t','clade_id'])\
                        .agg(#bfreq_mean=('bfreq','mean'),
                            #bfreq_std=('bfreq','std'),
                            n=('prop_richness','size'),
                            prop_richness_mean=('prop_richness','mean'), 
                            prop_richness_std=('prop_richness','std'), 
                            prop_shannon_mean=('prop_shannon','mean'),
                            prop_shannon_std=('prop_shannon','std'),
                            # prop_simpson=('prop_simpson','sum'),
                            clade_richness_mean=('clade_richness','mean'),
                            clade_richness_std=('clade_richness','std'),
                            clade_shannon_mean=('clade_shannon','mean'),
                            clade_shannon_std=('clade_shannon','std'),
                            clade_simpson_mean=('clade_simpson','mean'),
                            clade_simpson_std=('clade_simpson','std')).reset_index()
    ###
    cladeDiversity = cladeDiversity[cladeDiversity.n >= 2]
    return cladeDiversity

def computeVstats(runIDs,cID,conSim):
    print('...computing viral statistics...')
    virus_total = pd.read_sql_query("SELECT run_id, t, viral_abundance FROM summary \
                                    WHERE run_id in ({}) AND viral_abundance > 0"
                                    .format(', '.join(map(str, runIDs[runIDs.combo_id==cID]['run_id']))), conSim)\
        .rename(columns={"viral_abundance": "vtotal"})
    vstats = virus_total.groupby(['t'])\
        .agg(exp_vtotal=('vtotal', 'mean'), std_vtotal=('vtotal', 'std'), n=('vtotal', 'size'))\
        .reset_index()
    vstats = vstats[vstats.n >= 2]
    return vstats

def make_DFs(i,combos,scale,simDir,DBSIM_PATH,conSim,conSimNoV,figxy,vthreshold,bthreshold,micMutRep,hosts_per_strain, viruses_per_strain,\
                abundthresholdM,abundthresholdV,filterStrains,spacingM,spacingV,maxticksize):
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
        speciesMullerTreeHostVirus(True,sampleRunID, DBSIM_PATH, DBTREE_PATH, maxticksize, 
                            abundthresholdM, abundthresholdV,  filterStrains, 
                            spacingM, spacingV, figxy)
    immuneIDs = findImmuneIDs(runIDs,conSim)
    immuneIDsNoV = findImmuneIDs(runIDsNoV,conSimNoV)
    cladeIDs = findCladeIDs(runIDs,conSim)
    ##########
    noVdist, meanNoVDist = computeTraitDistribution(runIDsNoV,cIDNoV,conSimNoV)
    dist, meanDist = computeTraitDistribution(runIDs,cID,conSim)
    ## DIVERSITY
    diversityStats = computeSpacerDistribution(runIDs,immuneIDs,conSim)
    diversityStatsNoV = computeSpacerDistribution(runIDsNoV,immuneIDsNoV,conSimNoV)
    ## Clade Immune Proportions
    cladeDivAllN = computeCladeDiversity(runIDs,immuneIDs,cladeIDs,conSim)
    cladeDiversity = cladeDivAllN[cladeDivAllN['n']>9]
    #
    vstats = computeVstats(runIDs,cID,conSim)
    ####
    ####
    return vstats,meanDist,meanNoVDist,cladeDiversity,diversityStats,diversityStatsNoV,cID,cIDNoV,scale,\
                speciesColorDict,strainTimes,microbeStrainsDF

def make_figure(simDir,vstats,meanDist,meanNoVDist,cladeDiversity,diversityStats,diversityStatsNoV,ti,tf,cID,cIDNoV,scale,\
                resolve,speciesColorDict,strainTimes,microbeStrainsDF,rateName,diversity,propTitle,tree_parent_strain_id,tree_strain_id):
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
                    label=' '.join([r" $\sim$","{}".format(np.round(gRate,3))]),\
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
    axes[3].legend(handles,labels,loc='upper left', fontsize=8.5, title=propTitle, title_fontsize=8.5,ncol=4)
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
    fig.savefig(os.path.join(simDir,'diversity-cID{0}_scale{1}.pdf'.format(cID,scale)),dpi=resolve)
    plt.close('all')
    ##

def make_suppDFs2(i,combos,scale,simDir,DBSIM_PATH,conSim,conSimNoV,figxy,vthreshold,bthreshold,micMutRep,hosts_per_strain, viruses_per_strain,\
                  abundthresholdM,abundthresholdV,filterStrains,spacingM,spacingV,maxticksize):
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
        speciesMullerTreeHostVirus(True,sampleRunID, DBSIM_PATH, DBTREE_PATH, maxticksize, 
                            abundthresholdM, abundthresholdV,  filterStrains, 
                            spacingM, spacingV, figxy)
    immuneIDs = findImmuneIDs(runIDs,conSim)
    immuneIDsNoV = findImmuneIDs(runIDsNoV,conSimNoV)
    cladeIDs = findCladeIDs(runIDs,conSim)
    ##########
    noVdist, meanNoVDist = computeTraitDistribution(runIDsNoV,cIDNoV,conSimNoV)
    dist, meanDist = computeTraitDistribution(runIDs,cID,conSim)
    ## DIVERSITY
    diversityStats = computeSpacerDistribution(runIDs,immuneIDs,conSim)
    diversityStatsNoV = computeSpacerDistribution(runIDsNoV,immuneIDsNoV,conSimNoV)
    ## Clade Immune Proportions
    cladeDivAllN = computeCladeDiversity(runIDs,immuneIDs,cladeIDs,conSim)
    cladeDiversity = cladeDivAllN[cladeDivAllN['n']>9]
    #
    vstats = computeVstats(runIDs,cID,conSim)
    return vstats,meanDist,meanNoVDist,cladeDiversity,diversityStats,diversityStatsNoV,cID,cIDNoV,\
                     micMutSpacer, strainTimes,microbeStrainsDF,speciesColorDict

def make_suppfigure2(i,outer_grid,fig,resolve,simDir,scale,vstats,meanDist,meanNoVDist,cladeDiversity,diversityStats,diversityStatsNoV,ti,tf,cID,cIDNoV,\
                     micMutSpacer, strainTimes,microbeStrainsDF,tree_parent_strain_id,tree_strain_id,speciesColorDict,rateName,diversity,propTitle):
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
                    label=' '.join([r" $\sim$","{}".format(np.round(gRate,3))]),\
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
    axes[3].legend(handles,labels,loc='upper center', fontsize=8.5, title=propTitle, title_fontsize=8.5,ncol=4)        
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
    ##
    return fig
