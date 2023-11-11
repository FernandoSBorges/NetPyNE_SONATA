from netpyne import specs

cfg = specs.SimConfig()					            # object of class SimConfig to store simulation configuration
#------------------------------------------------------------------------------
#
# SIMULATION CONFIGURATION
#
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Run parameters
#------------------------------------------------------------------------------
cfg.duration = 1131.0 # Duration of the sim, in ms  
cfg.dt = 0.025 # Internal integration timestep to use
cfg.seeds = {'conn': 1333, 'stim': 1333, 'loc': 1333} 
cfg.hParams = {'celsius': 34, 'v_init': -70}  # -65
cfg.verbose = False
cfg.createNEURONObj = True
cfg.createPyStruct = True  
cfg.cvode_active = False
cfg.cvode_atol = 1e-6
cfg.cache_efficient = True

cfg.printRunTime = 0.1 # in sec	

cfg.includeParamsLabel = False
cfg.printPopAvgRates = True

cfg.checkErrors = False


#------------------------------------------------------------------------------
# Cells
#------------------------------------------------------------------------------
cfg.rootFolder = os.getcwd()

cfg.CellImportMode = 'template' # 'pkl'
cfg.celldiversity = True 

cfg.poptypeNumber = 8 # max 55 PAREI AQUI ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cfg.celltypeNumber = 12 # max 207
    
cfg.popParamLabels = ['L2e', 'L2i', 'L4e', 'L4i', 'L5e', 'L5i', 'L6e', 'L6i'] # to debug
cfg.cellParamLabels = ['L2e', 'L2i_BC', 'L2i_MC', 'L4e', 'L4i_BC', 'L4i_MC', 'L5e', 'L5i_BC', 'L5i_MC', 'L6e', 'L6i_BC', 'L6i_MC']


#------------------------------------------------------------------------------
# Options
#------------------------------------------------------------------------------

cfg.recordTraces = {'V_soma':{'sec':'soma','loc':0.5,'var':'v'}}  # Dict with traces to record
cfg.recordStep = 0.01 
cfg.seeds = {'conn': 1333, 'stim': 1333, 'loc': 1333} 
cfg.hParams = {'celsius': 34, 'v_init': -65}  
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

