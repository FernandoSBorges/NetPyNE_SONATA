from netpyne import specs

#############################################
####		SIMULATION PARAMETERS		#####
#############################################

#SIMDUR = STARTDEL + (THETA*8)	// simulation duration (msecs)
cfg = specs.SimConfig() 

cfg.dt = 0.025                 # Internal integration timestep to use
cfg.verbose = 0
cfg.duration = 2.0e3
cfg.recordStim = True
cfg.recordStep = 0.1             # Step size in ms to save data (e.g. V traces, LFP, etc)

cfg.printRunTime = 0.1
cfg.printPopAvgRates = True

cfg.seeds = {'conn': 1, 'stim': 1, 'loc': 1} # Seeds for randomizers (connectivity, input stimulation and cell locations)
cfg.hParams['celsius'] = 34.

cfg.allpops = ['Pyramidal','OLM','BS','Basket','AA']

cfg.recordTraces = {'V_soma':{'sec':'soma','loc':0.5,'var':'v'}} #, 'V_lmT':{'sec':'lm_thick1','loc':0.5,'var':'v'}}  # Dict with traces to record
cfg.analysis['plotTraces'] = {'include': [('Pyramidal',[10,50,90]),('AA',0),('Basket',[0,1]),('OLM',0),('BS',0)],'saveFig': True, 'oneFigPer':'trace', 'overlay': False, 'figSize':(18,12)}
cfg.analysis['plotRaster'] = {'saveFig': True, 'showFig': False, 'orderInverse': True, 'figSize': (18,12), 'labels': 'legend', 'popRates': True, 'fontSize':9, 'lw': 2, 'markerSize':4, 'marker': '.', 'dpi': 300} 
cfg.analysis['plot2Dnet'] = {'saveFig': True}
#------------------------------------------------------------------------------
# Saving
#------------------------------------------------------------------------------
cfg.simLabel = 'v0_batch0'
cfg.saveFolder = '../data/'+cfg.simLabel
# cfg.filename =                	## Set file output name
cfg.savePickle = False         	## Save pkl file
cfg.saveJson = True	           	## Save json file
cfg.saveDataInclude = ['simData'] ## 'simData' , 'simConfig', 'netParams'
cfg.backupCfgFile = None 		##  
cfg.gatherOnlySimData = False	##  
cfg.saveCellSecs = False			
cfg.saveCellConns = True	

#cfg.saveJson=True
#cfg.saveMat=True

# cfg.recordLFP = [[netParams.sizeX/2, netParams.sizeY*1/4, netParams.sizeZ/2], 
# 						[netParams.sizeX/2, netParams.sizeY*2/4, netParams.sizeZ/2],
# 						[netParams.sizeX/2, netParams.sizeY*3/4, netParams.sizeZ/2]]

# cfg.recordLFP = [[x,y,z] for x in range(900, netParams.sizeX, 900) for y in range(40, netParams.sizeY, 40) for z in range(40, netParams.sizeZ, 40)]
