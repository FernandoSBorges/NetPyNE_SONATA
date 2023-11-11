from netpyne import specs

cfg = specs.SimConfig()		

#------------------------------------------------------------------------------
# Options
#------------------------------------------------------------------------------

cfg = specs.SimConfig()					            # object of class SimConfig to store simulation configuration
cfg.duration = 1131.0			            # Duration of the simulation, in ms
cfg.dt = 0.01								                # Internal integration timestep to use
cfg.verbose = False							                # Show detailed messages 
cfg.recordTraces = {'V_soma':{'sec':'soma','loc':0.5,'var':'v'}}  # Dict with traces to record
cfg.recordStep = 0.01 
cfg.printRunTime = 0.1 # in sec			

cfg.seeds = {'conn': 1333, 'stim': 1333, 'loc': 1333} 
cfg.hParams = {'celsius': 34, 'v_init': -65}  
cfg.verbose = False
cfg.createNEURONObj = True
cfg.createPyStruct = True  
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
cfg.saveCellSecs = True			
cfg.saveCellConns = True	

#------------------------------------------------------------------------------
# ploting
#------------------------------------------------------------------------------
cfg.analysis['plotTraces'] = {'include': [0,1,2,3,4,5,6,7,8,9,10,11], 'timeRange': [400,1000], 'ylim': [-90,30], 'saveFig': True, 'showFig': True, 'figSize':(12,4)} # Plot recorded traces for this list of cells
cfg.analysis['plotShape'] = {'includePre': [n for n in range(0,24,1)],'includePost': [n for n in range(0,24,1)], 'includeAxon': False, 'saveFig': True, 'showFig': True, 'figSize':(22,22)}
cfg.analysis['plot2Dnet'] = {'view':'xy','saveFig': '../data/'+cfg.simLabel[0:9]+'/'+cfg.simLabel +'_xy_.png', 'showFig': True, 'figSize':(16,16), 'fontSize': 10}

