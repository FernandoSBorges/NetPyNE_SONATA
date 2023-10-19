"""
cfg.py 

Simulation configuration for A1 model (using NetPyNE)
This file has sim configs as well as specification for parameterized values in netParams.py 

Contributors: @gmail.com, @gmail.com
"""

from netpyne import specs
import pickle
import os

cfg = specs.SimConfig()  

#------------------------------------------------------------------------------
#
# SIMULATION CONFIGURATION
#
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Run parameters
#------------------------------------------------------------------------------
cfg.duration = 7.0*1e2 ## Duration of the sim, in ms  
cfg.dt = 0.05
# ~ cfg.seeds = {'conn': 4321, 'stim': 1234, 'loc': 4321} 
cfg.hParams = {'celsius': 34, 'v_init': -65}  
cfg.verbose = True
cfg.createNEURONObj = True
cfg.createPyStruct = True  
cfg.cvode_active = False
cfg.cvode_atol = 1e-6
cfg.cache_efficient = True
cfg.printRunTime = 0.5

cfg.includeParamsLabel = False
cfg.printPopAvgRates = True

cfg.checkErrors = False
#------------------------------------------------------------------------------
# Recording 
#------------------------------------------------------------------------------

allpops = []
cfg.cellsrec = 1
if cfg.cellsrec == 0:  cfg.recordCells = ['allpops'] # record all cells
elif cfg.cellsrec == 1: cfg.recordCells = [(pop,0) for pop in allpops] # record one cell of each pop
 
cfg.recordTraces = {'V_soma': {'sec':'soma_0', 'loc':0.5, 'var':'v'}}  ## Dict with traces to record
cfg.recordStim = True			
cfg.recordTime = True  		
cfg.recordStep = 0.1            



#------------------------------------------------------------------------------
# Saving
#------------------------------------------------------------------------------

cfg.simLabel = 'Ca1detailed'
cfg.saveFolder = '.'
# cfg.filename =                	## Set file output name
cfg.savePickle = False         	## Save pkl file
cfg.saveJson = True           	## Save json file
cfg.saveDataInclude = ['simConfig', 'netParams'] ## 'simData' , 'simConfig', 'netParams'
cfg.backupCfgFile = None 		##  
cfg.gatherOnlySimData = False	##  
cfg.saveCellSecs = True			##  
cfg.saveCellConns = True		##  

#------------------------------------------------------------------------------
# Analysis and plotting 
#------------------------------------------------------------------------------
cfg.analysis['plotRaster'] = {'include': allpops, 'timeRange': [0,cfg.duration], 'saveFig': True, 'showFig': False, 'popRates': True, 'orderInverse': True, 'figSize': (12,8), 'lw': 2.5, 'markerSize':5, 'marker': '.', 'dpi': 300} 
cfg.analysis['plotTraces'] = {'include': cfg.recordCells, 'timeRange': [0,cfg.duration], 'overlay': False, 'oneFigPer': 'trace', 'figSize': (12,8), 'saveFig': True, 'showFig': False} 

#------------------------------------------------------------------------------
# Cells
#------------------------------------------------------------------------------
cfg.cellmod =  {'L1_1': 'HH_full', 'L1_2': 'HH_full', 'L1_3': 'HH_full', 'L1_4': 'HH_full'}

#------------------------------------------------------------------------------
# Synapses
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Subcellular distribution
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Current inputs 
#------------------------------------------------------------------------------
cfg.addIClamp = 1

cfg.IClamp1 = {'pop': 'L1_1', 'sec': 'soma_0', 'loc': 0.5, 'start': 200, 'dur': 400, 'amp': -0.05}
cfg.IClamp2 = {'pop': 'L1_2', 'sec': 'soma_0', 'loc': 0.5, 'start': 200, 'dur': 400, 'amp': 0.15}
cfg.IClamp3 = {'pop': 'L1_2', 'sec': 'soma_0', 'loc': 0.5, 'start': 0, 'dur': 700, 'amp': -0.05}
cfg.IClamp4 = {'pop': 'L1_3', 'sec': 'soma_0', 'loc': 0.5, 'start': 200, 'dur': 400, 'amp': 0.25}
cfg.IClamp5 = {'pop': 'L1_3', 'sec': 'soma_0', 'loc': 0.5, 'start': 0, 'dur': 700, 'amp': -0.05}
cfg.IClamp6 = {'pop': 'L1_4', 'sec': 'soma_0', 'loc': 0.5, 'start': 200, 'dur': 400, 'amp': 0.35}
cfg.IClamp7 = {'pop': 'L1_4', 'sec': 'soma_0', 'loc': 0.5, 'start': 0, 'dur': 700, 'amp': -0.05}


#------------------------------------------------------------------------------
# NetStim inputs 
#------------------------------------------------------------------------------

## Attempt to add Background Noise inputs 
cfg.addNetStim = 0
