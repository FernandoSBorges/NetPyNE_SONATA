
"""
netParams.py

High-level specifications for M1 network model using NetPyNE

Contributors: salvadordura@gmail.com
"""

from netpyne import specs
import pickle, json
from neuron import gui
import matplotlib.pyplot as plt
import numpy as np
import os

netParams = specs.NetParams()   # object of class NetParams to store the network parameters


try:
    from __main__ import cfg  # import SimConfig object with params from parent module
except:
    from cfg import cfg
    
#------------------------------------------------------------------------------
# General connectivity parameters
#------------------------------------------------------------------------------
netParams.defaultThreshold = -10.0 # spike threshold, 10 mV is NetCon default, lower it for all cells
netParams.defaultDelay = 2.0 # default conn delay (ms)
netParams.propVelocity = 500.0 # propagation velocity (um/ms)

#------------------------------------------------------------------------------
# Cell parameters
#------------------------------------------------------------------------------
cellModels = ['HH_reduced', 'HH_full']
layer = {'Ori':[0.0, 0.1], 'Pyr':[0.1, 0.2], 'Rad':[0.2, 0.8], 'LM':[0.8, 1.0], 'Ca3':[2.0, 3.0], 'EC':[5.0, 6.0]}  # normalized layer boundaries


cellTypes = ['axoaxoniccell',  'bistratifiedcell',  'cckcell',  'ivycell',  'ngfcell',  'olmcell',  'poolosyncell',
 'pvbasketcell',  'scacell']

cellnumber=0
for cellType in cellTypes:

    if cellType == 'poolosyncell':
        cellRule = netParams.importCellParams(label='cell' + str(cellnumber), conds={'cellType':cellType, 'cellModel':'HH_full'}, 
            fileName='cell_data/class_' + cellName + '.hoc', cellName=cellType, cellInstance = True)
    else:
        cellRule = netParams.importCellParams(label='cell' + str(cellnumber), conds={'cellType':cellType, 'cellModel':'HH_reduced'}, 
            fileName='cell_data/class_' + cellName + '.hoc', cellName=cellType, cellInstance = True)

    cellnumber = cellnumber + 1
#------------------------------------------------------------------------------
# Population parameters
#------------------------------------------------------------------------------

netParams.popParams['olmcell'] = {'cellType': 'olmcell', 'ynormRange': layer['Ori'], 'cellModel': 'HH_reduced', 'numCells': 1} 

netParams.popParams['poolosyncell'] = {'cellType': 'poolosyncell', 'ynormRange': layer['Pyr'], 'cellModel': 'HH_full', 'numCells': 1} 
netParams.popParams['axoaxoniccell'] = {'cellType': 'axoaxoniccell', 'ynormRange': layer['Pyr'], 'cellModel': 'HH_reduced', 'numCells': 1} 
netParams.popParams['bistratifiedcell'] = {'cellType': 'bistratifiedcell', 'ynormRange': layer['Pyr'], 'cellModel': 'HH_reduced', 'numCells': 1} 
netParams.popParams['cckcell'] = {'cellType': 'cckcell', 'ynormRange': layer['Pyr'], 'cellModel': 'HH_reduced', 'numCells': 1} 
netParams.popParams['ivycell'] = {'cellType': 'ivycell', 'ynormRange': layer['Pyr'], 'cellModel': 'HH_reduced', 'numCells': 1} 
netParams.popParams['pvbasketcell'] = {'cellType': 'pvbasketcell', 'ynormRange': layer['Pyr'], 'cellModel': 'HH_reduced', 'numCells': 1} 

netParams.popParams['scacell'] = {'cellType': 'scacell', 'ynormRange': layer['Rad'], 'cellModel': 'HH_reduced', 'numCells': 1} 

netParams.popParams['ngfcell'] = {'cellType': 'ngfcell', 'ynormRange': layer['LM'], 'cellModel': 'HH_reduced', 'numCells': 1} 

#------------------------------------------------------------------------------
# Current inputs (IClamp)
#------------------------------------------------------------------------------
if cfg.addIClamp:
     for key in [k for k in dir(cfg) if k.startswith('IClamp')]:
        params = getattr(cfg, key, None)
        [pop,sec,loc,start,dur,amp] = [params[s] for s in ['pop','sec','loc','start','dur','amp']]

        #cfg.analysis['plotTraces']['include'].append((pop,0))  # record that pop

        # add stim source
        netParams.stimSourceParams[key] = {'type': 'IClamp', 'delay': start, 'dur': dur, 'amp': amp}
        
        # connect stim source to target
        netParams.stimTargetParams[key+'_'+pop] =  {
            'source': key, 
            'conds': {'pop': pop},
            'sec': sec, 
            'loc': loc}


