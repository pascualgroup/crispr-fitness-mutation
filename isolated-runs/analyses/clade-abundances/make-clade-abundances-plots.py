#!/usr/bin/env python3

import pandas as pd
import numpy as np
import scipy as sp
import math
import matplotlib.pyplot as plt
import sys
import os
import seaborn as sns
from scipy import stats
import sqlite3
# from scipy.interpolate import griddata
# from scipy.interpolate import interp1d
import matplotlib.ticker as ticker
import matplotlib.cm as cm
from matplotlib.colors import Normalize
import matplotlib.colors as mc
import colorsys

run_id = sys.argv[1]

resolve = 500

SCRIPT_PATH = os.path.abspath(os.path.dirname(__file__)) # cluster

DBSIM_PATH = os.path.join(SCRIPT_PATH,'..','..','..','simulation','sweep_db_gathered.sqlite') # cluster
# DBSIM_PATH = os.path.join('/Volumes','Yadgah','run_id1455_combo73_replicate15.sqlite')
# DBSIM_PATH = os.path.join('/Volumes','Yadgah','crispr-sweep-7-2-2022/isolates/runID1723-c35-r23/runID1723-c35-r23.sqlite')
DBCLADE_PATH = os.path.join(SCRIPT_PATH,'..','..','isolates',RUN_DIR,'clade-abundances_output.sqlite') # cluster
# DBCLADE_PATH = os.path.join('/Volumes','Yadgah','clade-abundances_output.sqlite') # local. run_id fixed; for testing
# DBCLADE_PATH = os.path.join('/Volumes','Yadgah','crispr-sweep-7-2-2022/isolates/runID1723-c35-r23','clade-abundances_output.sqlite') # local. run_id fixed; for testing

conSim = sqlite3.connect(DBSIM_PATH)
curSim = conSim.cursor()
conClade = sqlite3.connect(DBCLADE_PATH)
curClade = conClade.cursor()
ID = curSim.execute('SELECT combo_id,replicate FROM runs WHERE run_id = {}'.format(run_id)).fetchall()
combo_id = ID[0][0]
replicate = ID[0][1]
RUN_DIR = os.path.join('runID{0}-c{1}-r{2}'.format(run_id,combo_id,replicate))


print('SQLite Query: microbial abundance time series data')
microbeSim = pd.read_sql_query("SELECT t,bstrain_id,abundance FROM babundance WHERE run_id = {}".format(run_id), conSim)
microbe_stacked = microbeSim.pivot(index='t',columns='bstrain_id',values='abundance')
microbe_stacked.drop(2000.0,inplace=True)
print('SQLite Query: viral abundance time series data')
virusSim = pd.read_sql_query("SELECT t,vstrain_id,abundance FROM vabundance WHERE run_id = {}".format(run_id), conSim)
virus_stacked = virusSim.pivot(index='t',columns='vstrain_id',values='abundance')


print('SQLite Query: virus shannon data')
microbeClades = pd.read_sql_query("SELECT DISTINCT clade_id, bstrain_id \
FROM babundances", conClade)
microbeCladeAbundances = pd.read_sql_query("SELECT t, clade_id, abundance \
FROM clade_babundances", conClade)
microbeCladeIDs = pd.read_sql_query("SELECT DISTINCT clade_id \
FROM babundances", conClade)
virusClades = pd.read_sql_query("SELECT DISTINCT clade_id, vstrain_id \
FROM vabundances", conClade)
virusCladeIDs = pd.read_sql_query("SELECT DISTINCT clade_id \
FROM vabundances", conClade)
virusCladeAbundances = pd.read_sql_query("SELECT t, clade_id, abundance \
FROM clade_vabundances", conClade)


print('Compiling microbial clade abundance plots')
# Get a color map
Mcmap = cm.get_cmap('turbo')
# Get normalize function (takes data in range [vmin, vmax] -> [0, 1])
Mnorm = Normalize(vmin=1, vmax=len(microbeCladeIDs))
microbeColorDict = {}
cladeColorDict = {}
cladeDict = {}
colorInd = 0
for clade_id in microbeCladeIDs.values:
    # print(id[0])
    cladeColorDict[clade_id[0]] = colorInd
    cladeDict[clade_id[0]] = []
    colorInd += 1

for strain in microbe_stacked.columns.values:
    clade = microbeClades[microbeClades['bstrain_id']==strain]['clade_id'].values[0]
    cladeDict[clade] = np.append(cladeDict[clade],strain)

columnOrder = []
for clade_id in microbeCladeIDs.values:
    columnOrder = np.append(columnOrder,cladeDict[clade_id[0]])

columnOrder = columnOrder.astype(int)
microbe_stacked = microbe_stacked[columnOrder]

for strain in microbe_stacked.columns.values:
    clade = microbeClades[microbeClades['bstrain_id']==strain]['clade_id'].values[0]
    microbeColorDict[strain] = Mcmap(Mnorm(np.arange(1, len(microbeCladeIDs)+1, 1)))[cladeColorDict[clade]]



fig, ax = plt.subplots(1)
fig.suptitle('(run{0}-c{1}-r{2})'.format(run_id,combo_id,replicate))
microbe_stacked.plot.area(ax = ax,stacked=True,legend=False, linewidth=0,color=microbeColorDict,sort_columns=True)
microbe_stacked.plot(stacked=True, ax=ax, legend=False, color='white',sort_columns=True,linewidth=.3)
ax.set_ylabel(ylabel ='Microbial Strain Abundances',labelpad=15,fontsize=7)
ax.set_xlabel(xlabel = 'Time t',fontsize=7)
ax.ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(50))
lim = ax.get_ylim()
ax.set_ylim(0,lim[1])
fig.tight_layout()
# plt.show()
fig.savefig(os.path.join(SCRIPT_PATH,'..','..','isolates',RUN_DIR,'microbe-clade-abundances.png'),dpi=resolve)
# plt.show()


print('Compiling viral clade abundance plots')
# Get a color map
Vcmap = cm.get_cmap('turbo')
# Get normalize function (takes data in range [vmin, vmax] -> [0, 1])
Vnorm = Normalize(vmin=1, vmax=len(virusCladeIDs))
virusColorDict = {}
cladeColorDict = {}
shadeColorDict = {}
cladeDict = {}
colorInd = 0
for id in virusCladeIDs.values:
    # print(id[0])
    cladeColorDict[id[0]] = colorInd
    cladeDict[id[0]] = []
    colorInd += 1

for strain in virus_stacked.columns.values:
    clade = virusClades[virusClades['vstrain_id']==strain]['clade_id'].values[0]
    cladeDict[clade] = np.append(cladeDict[clade],strain)

columnOrder = []
for id in virusCladeIDs.values:
    columnOrder = np.append(columnOrder,cladeDict[id[0]])

columnOrder = columnOrder.astype(int)
virus_stacked = virus_stacked[columnOrder]

for strain in virus_stacked.columns.values:
    clade = virusClades[virusClades['vstrain_id']==strain]['clade_id'].values[0]
    print('strain: {}, clade: {}'.format(strain,clade))
    virusColorDict[strain] = Vcmap(Vnorm(np.arange(1, len(virusCladeIDs)+1, 1)))[cladeColorDict[clade]]

fig, ax = plt.subplots(1)
fig.suptitle('(run{0}-c{1}-r{2})'.format(run_id,combo_id,replicate))
virus_stacked.plot.area(ax = ax,stacked=True,legend=False, linewidth=0,color=virusColorDict,sort_columns=True)
virus_stacked.plot(stacked=True, ax=ax, legend=False, color='black',sort_columns=True,linewidth=.2)
ax.set_ylabel(ylabel ='Viral Strain Abundances',labelpad=15,fontsize=7)
ax.set_xlabel(xlabel = 'Time t',fontsize=7)
ax.ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(50))
lim = ax.get_ylim()
ax.set_ylim(0,lim[1])
fig.tight_layout()
fig.savefig(os.path.join(SCRIPT_PATH,'..','..','isolates',RUN_DIR,'virus-clade-abundances.png'),dpi=resolve)
# plt.show()


print('Compiling stacked microbial and viral clade abundance plots')
fig, ax = plt.subplots(2,sharex=True)
fig.suptitle('Strain Abundances (run{0}-c{1}-r{2})'.format(run_id,combo_id,replicate))
axes = [ax[0], ax[1]]
microbe_stacked.plot.area(ax = axes[0],stacked=True,legend=False, linewidth=0,color=microbeColorDict,sort_columns=True)
microbe_stacked.plot(stacked=True, ax=axes[0], legend=False, color='black',sort_columns=True,linewidth=.1)
axes[0].set_ylabel(ylabel ='Microbial Strain Abundances',labelpad=15,fontsize=7)
axes[0].set_xlabel(xlabel = 'Time t',fontsize=7)
axes[0].ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(50))
lim = axes[0].get_ylim()
axes[0].set_ylim(0,lim[1])
virus_stacked.plot.area(ax = axes[1],stacked=True,legend=False, linewidth=0,color=virusColorDict,sort_columns=True)
virus_stacked.plot(stacked=True, ax=axes[1], legend=False, color='black',sort_columns=True,linewidth=.1)
axes[1].set_ylabel(ylabel ='Viral Strain Abundances',labelpad=15,fontsize=7)
axes[1].set_xlabel(xlabel = 'Time t',fontsize=7)
axes[1].ticklabel_format(axis = 'y',style='sci',scilimits=(0,0))
axes[1].xaxis.set_minor_locator(ticker.MultipleLocator(50))
lim = axes[1].get_ylim()
axes[1].set_ylim(0,lim[1])
fig.tight_layout()
fig.savefig(os.path.join(SCRIPT_PATH,'..','..','isolates',RUN_DIR,'microbe-virus-clades-stacked-abundances.png'),dpi=resolve)
# plt.show()
