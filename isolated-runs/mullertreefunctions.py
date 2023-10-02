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

DBSIM_PATHS = []
DBSIM_PATHS.append(os.path.join(
    '/Volumes', 'Yadgah', 'sylvain-martin-collab/9_MOI3/isolates/runID1265-c7-r65/runID1265-c7-r65.sqlite'))
DBSIM_PATHS.append(os.path.join(
    '/Volumes', 'Yadgah','sylvain-martin-collab/9_MOI3/isolates/runID10804-c55-r4/runID10804-c55-r4.sqlite'))
DBTREE_PATHS = []
DBTREE_PATHS.append(os.path.join(
    '/Volumes', 'Yadgah', 'sylvain-martin-collab/9_MOI3/isolates/runID1265-c7-r65/trees_output.sqlite'))
DBTREE_PATHS.append(os.path.join(
    '/Volumes', 'Yadgah', 'sylvain-martin-collab/9_MOI3/isolates/runID10804-c55-r4/trees_output.sqlite'))
run_ids=[]
for i in range(0,len(DBSIM_PATHS)):
    conSim = sqlite3.connect(DBSIM_PATHS[i])
    curSim = conSim.cursor()
    run_id = curSim.execute('SELECT DISTINCT run_id FROM summary').fetchall()
    run_ids.append(run_id[0][0])


def get_inner_order(children_by_parent,identity):
    children_identities = children_by_parent.get(identity, [])
    if len(children_identities) == 0:
        return [identity, identity]
    inner = [get_inner_order(children_by_parent,c) for c in sorted(children_identities)]
    return [identity] + sum(inner, []) + [identity]

def speciesMullerPlotSideBySide(run_ids,species,DBSIM_PATHS,DBTREE_PATHS,treepalette,maxticksize,abundthreshold):
    if species == 'virus':
        strains = 'vstrains'
        strain = 'vstrain'
        strain_id = 'vstrain_id'
        abundanceTitle = 'Viral Strain Abundances'
        strain_abundance = "viral_abundance"
        straintotal = "vtotal"
        abundance = 'vabundance'
        tree_strain_id = "tree_vstrain_id"
        tree_parent_strain_id = "tree_parent_vstrain_id"
        parent_strain_id = "parent_vstrain_id"
    else:
        strains = 'bstrains'
        strain = 'bstrain'
        strain_id = 'bstrain_id'
        abundanceTitle = 'Host Strain\nAbundance'
        strain_abundance = "microbial_abundance"
        straintotal = "btotal"
        abundance = 'babundance'
        tree_strain_id = "tree_bstrain_id"
        tree_parent_strain_id = "tree_parent_bstrain_id"
        parent_strain_id = "parent_bstrain_id"
    fig = plt.figure()
    gs = fig.add_gridspec(2, 2, hspace=0, wspace=.15, width_ratios=[
                          1, 1], height_ratios=[1, 3])
    (ax1, ax2), (ax3, ax4) = gs.subplots(sharex='col')
    axes = [ax1, ax2, ax3, ax4]
    for i in range(0, len(run_ids)):
        run_id = run_ids[i]
        DBSIM_PATH = DBSIM_PATHS[i]
        DBTREE_PATH = DBTREE_PATHS[i]
        conSim = sqlite3.connect(DBSIM_PATH)
        curSim = conSim.cursor()
        ID = curSim.execute(
            'SELECT combo_id,replicate FROM runs WHERE run_id = {}'.format(run_id)).fetchall()
        combo_id = ID[0][0]
        replicate = ID[0][1]
        RUN_DIR = os.path.join(
            'runID{0}-c{1}-r{2}'.format(run_id, combo_id, replicate))
        species_stacked = pd.read_sql_query(
            "SELECT t,{0},abundance FROM {1} WHERE run_id = {2}".format(strain_id, abundance, run_id), conSim)
        t = max(pd.read_sql_query(
            "SELECT t FROM vabundance WHERE run_id = {}".format(run_id), conSim).t.values)
        species_stacked = species_stacked[species_stacked.t <= t]
        species_total = pd.read_sql_query("SELECT t,{0} FROM summary WHERE run_id = {1}".format(strain_abundance, run_id), conSim)\
            .rename(columns={strain_abundance: straintotal})
        species_total = species_total[species_total.t <= max(
            species_stacked.t)]
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
            FROM {0} WHERE {1} in ({2})".format(strains, strain_id, ', '.join(map(str, parentTrack))), conSim)
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
        strainTimes = pd.read_sql_query(
            "SELECT tree_{0}, t_creation, t_extinction, tree_parent_{0} \
        FROM tree_{1}_creation_extinction WHERE tree_{0} in ({2})".format(
                strain_id, strain, ', '.join(map(str, keepTreeStrains))), conTree)
        strainTimesNoExt = pd.read_sql_query(
            "SELECT {0}, t_creation, parent_{0} \
        FROM {1} WHERE {0} in ({2})".format(
                strain_id, strains, ', '.join(map(str, keepTreeStrainsDF['bstrain_id']))), conSim)\
            .merge(keepTreeStrainsDF, on=strain_id).drop(columns=[strain_id])\
            .rename(columns={'parent_{0}'.format(strain_id): strain_id})\
            .merge(keepTreeStrainsDF.rename(columns={tree_strain_id: 'tree_parent_{0}'.format(strain_id)}), on=strain_id)\
            .drop(columns=[strain_id])
        strainTimesNoExt = strainTimesNoExt[[(i not in list(strainTimes[tree_strain_id]))
                                            for i in list(strainTimesNoExt[tree_strain_id])]]
        strainTimesNoExt['t_extinction'] = len(strainTimesNoExt)*[t]
        strainTimes = pd.concat(
            [strainTimes, strainTimesNoExt]).reset_index(drop=True)
        treeAbundances = pd.read_sql_query(
            "SELECT t, tree_{0}, abundance \
        FROM tree_{1} WHERE tree_{0} in ({2})"
            .format(strain_id, abundance, ', '.join(map(str, keepTreeStrains))), conTree)
        newTreeOrder = {}
        newTreeID = 1
        for treeStrain in keepTreeStrains[::-1]:
            newTreeOrder[treeStrain] = newTreeID
            newTreeID += 1
        treeAbundances.replace({tree_strain_id: newTreeOrder}, inplace=True)
        strainTimes.replace({tree_strain_id: newTreeOrder,
                            tree_parent_strain_id: newTreeOrder}, inplace=True)
        numStrains = len(keepTreeStrains)
        #
        Vcmap = cm.get_cmap(treepalette)
        # Vcmap = sns.color_palette("icefire",as_cmap=True)
        norm = Normalize(vmin=float(1), vmax=float(max(newTreeOrder.values())))
        speciesColorDict = {}
        maxAbundanceSpecies = np.max(treeAbundances.abundance.values)
        markerIncSpecies = maxticksize/maxAbundanceSpecies
        markerColorsSpecies = []
        hlinecSpecies = []
        vlinecSpecies = []
        hcolorsSpecies = []
        vcolorsSpecies = []
        print('Compiling stacked species abundances and tree plots')
        for strainID in sorted(treeAbundances[tree_strain_id].unique()):
            # print(strainID)
            tCreate = strainTimes[strainTimes[tree_strain_id]
                                  == strainID].t_creation.values[0]
            tExtinct = strainTimes[strainTimes[tree_strain_id]
                                   == strainID].t_extinction.values[0]
            parent = strainTimes[strainTimes[tree_strain_id]
                                 == strainID][tree_parent_strain_id].values[0]
            hlinecSpecies.append([[tCreate, strainID], [tExtinct, strainID]])
            vlinecSpecies.append([[tCreate, parent], [tCreate, strainID]])
            speciesColorDict[strainID] = Vcmap(norm(strainID))
            hcolorsSpecies.append(speciesColorDict[strainID])
            vcolorsSpecies.append(speciesColorDict[strainID])
            markerColorsSpecies.append(speciesColorDict[strainID])
            axes[i+2].scatter(treeAbundances[treeAbundances[tree_strain_id] == strainID]['t'],
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
        axes[i+2].add_collection(strainLineages)
        axes[i+2].add_collection(creationLines)
        axes[i+2].set_ylim(0, np.max(strainTimes[tree_strain_id])+1)
        axes[i+2].set_xlim(0, t)
        axes[i+2].xaxis.set_minor_locator(ticker.MultipleLocator(25))
        axes[i+2].set_yticks([])
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
            .rename(columns={'tree_bstrain_id': 'Identity', 'tree_parent_bstrain_id': 'Parent'})
        check = check[check.Parent != 0]
        children_by_parent = check.groupby(
            'Parent')['Identity'].apply(lambda x: list(sorted(x)))
        order = []
        for mrca in sorted(treeAbundances[treeAbundances.t == 0][tree_strain_id].values):
            if mrca not in order:
                order += get_inner_order(children_by_parent, mrca)
        mullerAbunds[order].plot.area(
            ax=axes[i], stacked=True, legend=False, linewidth=0, color=speciesColorDict, sort_columns=True)
        # mullerAbunds[order].plot(stacked=True, ax=axes[i], legend=False, color='white',sort_columns=True,linewidth=.1)
        if species != 'virus':
            virus_total = pd.read_sql_query("SELECT t, viral_abundance FROM summary WHERE run_id = {}".format(run_id), conSim)\
                .rename(columns={'viral_abundance': 'vtotal'})
            virus_total = virus_total[virus_total.t <= t]
            # virus_total['vtotal'] = np.log10(virus_total['vtotal'])
            axes.append(axes[i].twinx())
            axes[i+4].plot(virus_total['t'], virus_total['vtotal'],
                           linewidth=0, color='grey')
            axes[i+4].fill_between(
                virus_total['t'], virus_total['vtotal'], color='grey', alpha=0.75)
            # maxO = int(np.ceil(virus_total['vtotal']))
            # axes[i+4].set_yticks([i for i in range(0,maxO+1,1)])
            # axes[i+4].set_yticklabels([r'$10^{}$'.format(i) for i in range(0,maxO+1,1)])
            axes[i+4].set_ylabel(ylabel='Viral\nAbundance',
                                 labelpad=35, fontsize=12, rotation=270)
            axes[i+4].tick_params(axis='y', labelsize=12)
            axes[i+4].yaxis.get_offset_text().set_fontsize(12)
            lim = axes[i+4].get_ylim()
            axes[i+4].set_ylim(0, lim[1])
            if i == 0:
                axes[i+4].set_ylabel(ylabel='')
        axes[i].set_ylabel(ylabel=abundanceTitle, labelpad=15, fontsize=12)
        if i == 1:
            axes[i].set_ylabel(ylabel='')
        axes[i+2].set_xlabel(xlabel='Time t', labelpad=10, fontsize=12)
        # axes[i].ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
        maxO = int(np.ceil(np.log10(max(species_total[straintotal]))))
        axes[i].set_yticks([i for i in range(0, maxO+1, 1)])
        axes[i].set_yticklabels([r'$10^{}$'.format(i)
                                for i in range(0, maxO+1, 1)])
        axes[i].xaxis.set_minor_locator(ticker.MultipleLocator(25))
        lim = axes[i].get_ylim()
        axes[i].set_ylim(0, lim[1])
        axes[i].set_xlim(0, t)
        axes[i+2].tick_params(axis='x', labelsize=12)
        axes[i].tick_params(axis='y', labelsize=12)
        axes[i+2].tick_params(axis='y', labelsize=12)
        # axes[i].yaxis.get_offset_text().set_fontsize(20)
    plt.show()
    return fig, axes


DBSIM_PATH = os.path.join(
    '/Volumes', 'Yadgah', 'sylvain-martin-collab/12_MOI3/isolates/runID4209-c15-r9/runID4209-c15-r9.sqlite')
DBTREE_PATH = os.path.join(
    '/Volumes', 'Yadgah', 'sylvain-martin-collab/12_MOI3/isolates/runID4209-c15-r9/trees_output.sqlite')
conSim = sqlite3.connect(DBSIM_PATH)
curSim = conSim.cursor()
run_id = curSim.execute('SELECT DISTINCT run_id FROM summary').fetchall()
run_id = run_id[0][0]


def speciesMullerPlotSingle(run_id, species, DBSIM_PATH, DBTREE_PATH, treepalette, maxticksize, abundthreshold):
    if species == 'virus':
        strains = 'vstrains'
        strain = 'vstrain'
        strain_id = 'vstrain_id'
        abundanceTitle = 'Viral Strain Abundances'
        strain_abundance = "viral_abundance"
        straintotal = "vtotal"
        abundance = 'vabundance'
        tree_strain_id = "tree_vstrain_id"
        tree_parent_strain_id = "tree_parent_vstrain_id"
        parent_strain_id = "parent_vstrain_id"
    else:
        strains = 'bstrains'
        strain = 'bstrain'
        strain_id = 'bstrain_id'
        abundanceTitle = 'Host Strain\nAbundance'
        strain_abundance = "microbial_abundance"
        straintotal = "btotal"
        abundance = 'babundance'
        tree_strain_id = "tree_bstrain_id"
        tree_parent_strain_id = "tree_parent_bstrain_id"
        parent_strain_id = "parent_bstrain_id"
    fig = plt.figure()
    gs = fig.add_gridspec(2, 1, hspace=0, wspace=.15, height_ratios=[1, 3])
    ax1, ax2 = gs.subplots(sharex='col')
    axes = [ax1, ax2]
    ##
    conSim = sqlite3.connect(DBSIM_PATH)
    curSim = conSim.cursor()
    ID = curSim.execute(
        'SELECT combo_id,replicate FROM runs WHERE run_id = {}'.format(run_id)).fetchall()
    combo_id = ID[0][0]
    replicate = ID[0][1]
    RUN_DIR = os.path.join(
        'runID{0}-c{1}-r{2}'.format(run_id, combo_id, replicate))
    species_stacked = pd.read_sql_query(
        "SELECT t,{0},abundance FROM {1} WHERE run_id = {2}".format(strain_id, abundance, run_id), conSim)
    t = max(pd.read_sql_query(
        "SELECT t FROM vabundance WHERE run_id = {}".format(run_id), conSim).t.values)
    species_stacked = species_stacked[species_stacked.t <= t]
    species_total = pd.read_sql_query("SELECT t,{0} FROM summary WHERE run_id = {1}".format(strain_abundance, run_id), conSim)\
        .rename(columns={strain_abundance: straintotal})
    species_total = species_total[species_total.t <= max(
        species_stacked.t)]
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
        FROM {0} WHERE {1} in ({2})".format(strains, strain_id, ', '.join(map(str, parentTrack))), conSim)
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
    strainTimes = pd.read_sql_query(
        "SELECT tree_{0}, t_creation, t_extinction, tree_parent_{0} \
    FROM tree_{1}_creation_extinction WHERE tree_{0} in ({2})".format(
            strain_id, strain, ', '.join(map(str, keepTreeStrains))), conTree)
    strainTimesNoExt = pd.read_sql_query(
        "SELECT {0}, t_creation, parent_{0} \
    FROM {1} WHERE {0} in ({2})".format(
            strain_id, strains, ', '.join(map(str, keepTreeStrainsDF[strain_id]))), conSim)\
        .merge(keepTreeStrainsDF, on=strain_id).drop(columns=[strain_id])\
        .rename(columns={'parent_{0}'.format(strain_id): strain_id})\
        .merge(keepTreeStrainsDF.rename(columns={tree_strain_id: 'tree_parent_{0}'.format(strain_id)}), on=strain_id)\
        .drop(columns=[strain_id])
    strainTimesNoExt = strainTimesNoExt[[(i not in list(strainTimes[tree_strain_id]))
                                        for i in list(strainTimesNoExt[tree_strain_id])]]
    strainTimesNoExt['t_extinction'] = len(strainTimesNoExt)*[t]
    strainTimes = pd.concat(
        [strainTimes, strainTimesNoExt]).reset_index(drop=True)
    treeAbundances = pd.read_sql_query(
        "SELECT t, tree_{0}, abundance \
    FROM tree_{1} WHERE tree_{0} in ({2})"
        .format(strain_id, abundance, ', '.join(map(str, keepTreeStrains))), conTree)
    newTreeOrder = {}
    newTreeID = 1
    for treeStrain in keepTreeStrains[::-1]:
        newTreeOrder[treeStrain] = newTreeID
        newTreeID += 1
    treeAbundances.replace({tree_strain_id: newTreeOrder}, inplace=True)
    strainTimes.replace({tree_strain_id: newTreeOrder,
                        tree_parent_strain_id: newTreeOrder}, inplace=True)
    numStrains = len(keepTreeStrains)
    #
    Vcmap = cm.get_cmap(treepalette)
    # Vcmap = sns.color_palette("icefire",as_cmap=True)
    norm = Normalize(vmin=float(1), vmax=float(max(newTreeOrder.values())))
    speciesColorDict = {}
    maxAbundanceSpecies = np.max(treeAbundances.abundance.values)
    markerIncSpecies = maxticksize/maxAbundanceSpecies
    markerColorsSpecies = []
    hlinecSpecies = []
    vlinecSpecies = []
    hcolorsSpecies = []
    vcolorsSpecies = []
    print('Compiling stacked species abundances and tree plots')
    for strainID in sorted(treeAbundances[tree_strain_id].unique()):
        # print(strainID)
        tCreate = strainTimes[strainTimes[tree_strain_id]
                                == strainID].t_creation.values[0]
        tExtinct = strainTimes[strainTimes[tree_strain_id]
                                == strainID].t_extinction.values[0]
        parent = strainTimes[strainTimes[tree_strain_id]
                                == strainID][tree_parent_strain_id].values[0]
        hlinecSpecies.append([[tCreate, strainID], [tExtinct, strainID]])
        vlinecSpecies.append([[tCreate, parent], [tCreate, strainID]])
        speciesColorDict[strainID] = Vcmap(norm(strainID))
        hcolorsSpecies.append(speciesColorDict[strainID])
        vcolorsSpecies.append(speciesColorDict[strainID])
        markerColorsSpecies.append(speciesColorDict[strainID])
        axes[1].scatter(treeAbundances[treeAbundances[tree_strain_id] == strainID]['t'],
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
    axes[1].add_collection(strainLineages)
    axes[1].add_collection(creationLines)
    axes[1].set_ylim(0, np.max(strainTimes[tree_strain_id])+1)
    axes[1].set_xlim(0, t)
    axes[1].xaxis.set_minor_locator(ticker.MultipleLocator(25))
    axes[1].set_yticks([])
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
        ax=axes[0], stacked=True, legend=False, linewidth=0, color=speciesColorDict, sort_columns=True)
    # mullerAbunds[order].plot(stacked=True, ax=axes[0], legend=False, color='white',sort_columns=True,linewidth=.1)
    if species == 'microbe':
        virus_total = pd.read_sql_query("SELECT t, viral_abundance FROM summary WHERE run_id = {}".format(run_id), conSim)\
            .rename(columns={'viral_abundance': 'vtotal'})
        virus_total = virus_total[virus_total.t <= t]
        # virus_total['vtotal'] = np.log10(virus_total['vtotal'])
        axes.append(axes[0].twinx())
        axes[2].plot(virus_total['t'], virus_total['vtotal'],
                        linewidth=0, color='grey')
        axes[2].fill_between(
            virus_total['t'], virus_total['vtotal'], color='grey', alpha=0.8)
        # maxO = int(np.ceil(virus_total['vtotal']))
        # axes[i+4].set_yticks([i for i in range(0,maxO+1,1)])
        # axes[i+4].set_yticklabels([r'$10^{}$'.format(i) for i in range(0,maxO+1,1)])
        axes[2].set_ylabel(ylabel='Viral\nAbundance',
                                labelpad=35, fontsize=12, rotation=270)
        axes[2].tick_params(axis='y', labelsize=12)
        axes[2].yaxis.get_offset_text().set_fontsize(12)
        lim = axes[2].get_ylim()
        axes[2].set_ylim(0, lim[1])
    axes[0].set_ylabel(ylabel=abundanceTitle, labelpad=15, fontsize=12)
    axes[1].set_xlabel(xlabel='Time t', labelpad=10, fontsize=12)
    # axes[0].ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
    maxO = int(np.ceil(np.log10(max(species_total[straintotal]))))
    axes[0].set_yticks([i for i in range(0, maxO+1, 1)])
    axes[0].set_yticklabels([r'$10^{}$'.format(i)
                            for i in range(0, maxO+1, 1)])
    axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(25))
    # axes[0].set_xscale('log', base=10)
    # axes[1].set_xscale('log', base=10)
    # axes[2].set_xscale('log', base=10)
    lim = axes[0].get_ylim()
    axes[0].set_ylim(0, lim[1])
    lim = axes[0].get_xlim()
    axes[0].set_xlim(lim[0], t)
    # axes[1].set_xscale('log', base=10)
    # axes[2].set_xscale('log', base=10)
    axes[1].tick_params(axis='x', labelsize=12)
    axes[0].tick_params(axis='y', labelsize=12)
    axes[1].tick_params(axis='y', labelsize=12)
    # axes[0].yaxis.get_offset_text().set_fontsize(20)
    plt.show()
    return fig, axes





DBSIM_PATH = os.path.join(
    '/Volumes', 'Yadgah', 'sylvain-martin-collab/12_MOI3/isolates/runID4209-c15-r9/runID4209-c15-r9.sqlite')
DBTREE_PATH = os.path.join(
    '/Volumes', 'Yadgah', 'sylvain-martin-collab/12_MOI3/isolates/runID4209-c15-r9/trees_output.sqlite')
conSim = sqlite3.connect(DBSIM_PATH)
curSim = conSim.cursor()
run_id = curSim.execute('SELECT DISTINCT run_id FROM summary').fetchall()
run_id = run_id[0][0]


def speciesMullerPlotBoth(run_id, species, DBSIM_PATH, DBTREE_PATH, treepalette, maxticksize, abundthreshold):
    fig = plt.figure()
    gs = fig.add_gridspec(7, 1, hspace=0, wspace=.15,
                        height_ratios=[1, 0.1, 1, 0.3, 3, .1, 3])
    ax1, ax6, ax2, ax3, ax4, ax7, ax5 = gs.subplots(sharex='col')
    ax3.remove()
    ax6.remove()
    ax7.remove()
    axes = [ax1, ax2, ax4, ax5]
    ##
    conSim = sqlite3.connect(DBSIM_PATH)
    curSim = conSim.cursor()
    ID = curSim.execute(
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
            treepalette = 'plasma_r'
            Vcmap = cm.get_cmap(treepalette)
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
            treepalette = 'viridis'
            Vcmap = cm.get_cmap(treepalette)
            i = 0
        species_stacked = pd.read_sql_query(
            "SELECT t,{0},abundance FROM {1} WHERE run_id = {2}".format(strain_id, abundance, run_id), conSim)
        t = max(pd.read_sql_query(
            "SELECT t FROM vabundance WHERE run_id = {}".format(run_id), conSim).t.values)
        species_stacked = species_stacked[species_stacked.t <= t]
        species_total = pd.read_sql_query("SELECT t,{0} FROM summary WHERE run_id = {1}".format(strain_abundance, run_id), conSim)\
            .rename(columns={strain_abundance: straintotal})
        species_total = species_total[species_total.t <= max(
            species_stacked.t)]
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
            FROM {0} WHERE {1} in ({2})".format(strains, strain_id, ', '.join(map(str, parentTrack))), conSim)
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
        strainTimes = pd.read_sql_query(
            "SELECT tree_{0}, t_creation, t_extinction, tree_parent_{0} \
        FROM tree_{1}_creation_extinction WHERE tree_{0} in ({2})".format(
                strain_id, strain, ', '.join(map(str, keepTreeStrains))), conTree)
        strainTimesNoExt = pd.read_sql_query(
            "SELECT {0}, t_creation, parent_{0} \
        FROM {1} WHERE {0} in ({2})".format(
                strain_id, strains, ', '.join(map(str, keepTreeStrainsDF[strain_id]))), conSim)\
            .merge(keepTreeStrainsDF, on=strain_id).drop(columns=[strain_id])\
            .rename(columns={'parent_{0}'.format(strain_id): strain_id})\
            .merge(keepTreeStrainsDF.rename(columns={tree_strain_id: 'tree_parent_{0}'.format(strain_id)}), on=strain_id)\
            .drop(columns=[strain_id])
        strainTimesNoExt = strainTimesNoExt[[(i not in list(strainTimes[tree_strain_id]))
                                            for i in list(strainTimesNoExt[tree_strain_id])]]
        strainTimesNoExt['t_extinction'] = len(strainTimesNoExt)*[t]
        strainTimes = pd.concat(
            [strainTimes, strainTimesNoExt]).reset_index(drop=True)
        treeAbundances = pd.read_sql_query(
            "SELECT t, tree_{0}, abundance \
        FROM tree_{1} WHERE tree_{0} in ({2})"
            .format(strain_id, abundance, ', '.join(map(str, keepTreeStrains))), conTree)
        newTreeOrder = {}
        newTreeID = 1
        for treeStrain in keepTreeStrains[::-1]:
            newTreeOrder[treeStrain] = newTreeID
            newTreeID += 1
        treeAbundances.replace({tree_strain_id: newTreeOrder}, inplace=True)
        strainTimes.replace({tree_strain_id: newTreeOrder,
                            tree_parent_strain_id: newTreeOrder}, inplace=True)
        numStrains = len(keepTreeStrains)
        #
        # Vcmap = sns.color_palette("icefire",as_cmap=True)
        norm = Normalize(vmin=float(1), vmax=float(max(newTreeOrder.values())))
        speciesColorDict = {}
        maxAbundanceSpecies = np.max(treeAbundances.abundance.values)
        markerIncSpecies = maxticksize/maxAbundanceSpecies
        markerColorsSpecies = []
        hlinecSpecies = []
        vlinecSpecies = []
        hcolorsSpecies = []
        vcolorsSpecies = []
        print('Compiling stacked species abundances and tree plots')
        for strainID in sorted(treeAbundances[tree_strain_id].unique()):
            # print(strainID)
            tCreate = strainTimes[strainTimes[tree_strain_id]
                                == strainID].t_creation.values[0]
            tExtinct = strainTimes[strainTimes[tree_strain_id]
                                == strainID].t_extinction.values[0]
            parent = strainTimes[strainTimes[tree_strain_id]
                                == strainID][tree_parent_strain_id].values[0]
            hlinecSpecies.append([[tCreate, strainID], [tExtinct, strainID]])
            vlinecSpecies.append([[tCreate, parent], [tCreate, strainID]])
            speciesColorDict[strainID] = Vcmap(norm(strainID))
            hcolorsSpecies.append(speciesColorDict[strainID])
            vcolorsSpecies.append(speciesColorDict[strainID])
            markerColorsSpecies.append(speciesColorDict[strainID])
            axes[i+2].scatter(treeAbundances[treeAbundances[tree_strain_id] == strainID]['t'],
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
        axes[i+2].add_collection(strainLineages)
        axes[i+2].add_collection(creationLines)
        axes[i+2].set_ylim(0, np.max(strainTimes[tree_strain_id])+1)
        axes[i+2].set_xlim(0, t)
        axes[i+2].xaxis.set_minor_locator(ticker.MultipleLocator(25))
        axes[i+2].set_yticks([])
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
            ax=axes[i], stacked=True, legend=False, linewidth=0, color=speciesColorDict, sort_columns=True)
        axes[i].set_ylabel(ylabel=abundanceTitle, labelpad=15, fontsize=12)
        axes[i+2].set_ylabel(ylabel=strainTitle, labelpad=15, fontsize=12)
        axes[i+2].set_xlabel(xlabel='Time t', labelpad=10, fontsize=12)
        maxO = int(np.ceil(np.log10(max(species_total[straintotal]))))
        axes[i].set_yticks([i for i in range(0, maxO+1, 1)])
        yticklabels = [r'$10^{}$'.format(i)
                       for i in range(0, maxO, 1)]
        yticklabels.append('')
        axes[i].set_yticklabels(yticklabels)
        axes[i].xaxis.set_minor_locator(ticker.MultipleLocator(25))
        lim = axes[i].get_ylim()
        axes[i].set_ylim(0, lim[1])
        lim = axes[i].get_xlim()
        axes[i].set_xlim(lim[0], t)
        axes[i+2].tick_params(axis='x', labelsize=12)
        axes[i].tick_params(axis='y', labelsize=7)
    plt.show()
    return fig, axes
