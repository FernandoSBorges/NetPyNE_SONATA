"""
cfg.py 

Contributors: joao.moreira@downstate.edu, salvadordura@gmail.com
"""

from netpyne import specs
import pickle
import NetPyNE_BBP
import numpy as np
import os

#------------------------------------------------------------------------------
# Simulation options
#------------------------------------------------------------------------------
cfg = specs.SimConfig()     # object of class SimConfig to store simulation configuration

cfg.duration        = 0.001*1e3              # Duration of the simulation, in ms
cfg.dt              = 0.025                 # Internal integration timestep to use
cfg.hParams         = {'v_init': -79, 'celsius': 34} 
cfg.verbose         = False                 # Show detailed messages
cfg.printRunTime    = 0.1

cfg.coreneuron = False
cfg.random123 = True

cfg.recordTraces = {'V_soma': {'sec': 'soma_0', 'loc': 0.5, 'var': 'v'}}

# cfg.recordStep = 1              # Step size in ms to save data (eg. V traces, LFP, etc)
cfg.recordStep = cfg.dt         # Step size in ms to save data (eg. V traces, LFP, etc)

cfg.recordStim = False  
cfg.recordTime = False

#------------------------------------------------------------------------------
# Load BBP circuit
#------------------------------------------------------------------------------

cfg.base_dir = os.path.expanduser("~")

cfg.BBP_rootFolder                  = cfg.base_dir+'/Documents/BBP_thalamus_microcircuit_2'
cfg.sonataConfigFile                = cfg.BBP_rootFolder+'/sonata/circuit_sonata.json'
cfg.morphologyFolder_h5             = cfg.BBP_rootFolder+'/sonata/morphologies_h5'
cfg.morphologyFolder_asc            = cfg.BBP_rootFolder+'/sonata/morphologies/morphologies_asc'
cfg.virtual_input_spikes            = cfg.BBP_rootFolder+'/sonata/simulation_spike_files'
cfg.ct_virtual_noise                = cfg.virtual_input_spikes+'/input_spikes_ct_noise.dat'
cfg.ml_virtual_noise                = cfg.virtual_input_spikes+'/input_spikes_ml_noise.dat'

# --- Path to NetPyNE files
cfg.NetPyNE_rootFolder              = cfg.base_dir+'/Documents/thalamus_netpyne'
cfg.NetPyNE_sim                     = cfg.NetPyNE_rootFolder+'/sim'
cfg.NetPyNE_templateCells           = cfg.NetPyNE_rootFolder+'/mod'
cfg.NetPyNE_exportedCells           = cfg.NetPyNE_rootFolder+'/cells/morphologies_swc'
cfg.NetPyNE_JSON_cells              = cfg.NetPyNE_rootFolder+'/cells/netpyne_morphologies'
cfg.NetPyNE_conn                    = cfg.NetPyNE_rootFolder+'/conn'
cfg.NetPyNE_data                    = cfg.NetPyNE_rootFolder+'/data'
cfg.NetPyNE_node_pathway_edges      = cfg.NetPyNE_conn+'/node_pathway_edges'
cfg.NetPyNE_input_noise_savePath    = cfg.NetPyNE_rootFolder+'/conn/external_inputs'

#------------------------------------------------------------------------------
# Sim Properties
#------------------------------------------------------------------------------

cfg.select_microcircuit = None

# --- Convert cell models
cfg.convertCellMorphologies=False

cfg.target = 'VPL_TC'
gids = [36636]; 

iclamp_amp = 0.0061 # 0.0068 # 0.0072 # 0.008 # 0.01 # 0.03 # 0.061 # 'VPL_TC'     

cfg.select_thal_gids = NetPyNE_BBP.ConfigCfg.sampleGids(pop=cfg.target, cfg_file=cfg.sonataConfigFile, gids = gids)

# --- Cell Model
cfg.loadCellModel = False;     
cfg.convertCellModel= True;    
cfg.saveCellModel   = True; 
cfg.plotSynLocation = True; 
cfg.plotCellShape = True

# --- Connectivity Type
cfg.connType = 'original'    # adds original edge connectivity to the target cell model

# --- Connectivity
cfg.th_updateConns      = False     # (True): generates the edge connectivity to the target cell model | (False): Loads from file | Defaults to True if file doesn't exist
cfg.th_topological      = False      # Redistributes the synapses onto the dendrites according to Driver/IN/RT/CT hierarchy
cfg.th_useTableVals     = False     # Replace edges values for the values in Table S2 from the paper
cfg.saveIndividualNoise = True      # Saves individual input noise files for each postsynaptic cell for later plotting
cfg.th_connOverrides    = True      # Alters the edge properties to match the parameters that were changed before model run, but are not stored in the original model
cfg.removeChemical      = False
cfg.removeElectrical    = True      # wont work with TRN inputs as vecstims
cfg.th_singleTcPop      = False

cfg.modType          = 'Prob_original'
cfg.th_stimType         = 'VecStim'

if cfg.target == 'VPL_TC':
    cfg.th_inter_sources    = ['Rt_RC','VPL_IN']
    cfg.th_select_pathways  = [
                                'MedialLemniscus_projections',
                                'CorticoThalamic_projections', 
                                'thalamus_neurons|Rt_RC',
                                'thalamus_neurons|VPL_IN'
                                ]
elif cfg.target == 'Rt_RC':
    cfg.th_inter_sources    = ['VPL_TC','Rt_RC']
    cfg.th_select_pathways  = [
                                'CorticoThalamic_projections', 
                                'thalamus_neurons|Rt_RC',
                                # 'thalamus_neurons|VPL_IN'
                                'thalamus_neurons|VPL_TC'
                                ]
# --- Extracellular Calcium concentration

cfg.re_rescale = 1

cfg.cao_secs            = 1.2
cfg.rescaleUSE          = 0.4029343148532312 * cfg.re_rescale # From BlueConfig file
# --- Change conn weights
cfg.rescale_conn_weight = { 'MedialLemniscus_projections':  1,
                            'CorticoThalamic_projections':  1,
                            'thalamus_neurons|Rt_RC':       1,
                            'thalamus_neurons|VPL_IN':      1,
                            'thalamus_neurons|VPL_TC':      1,
                            }
                    
# --- Change biophysics
cfg.removeExtraCurrents=False
cfg.changeCurrents  = []
cfg.rescaleCurrents = []
# --- Stims
cfg.add_current_stims_noise   = False
noise_mean = 0
noise_var = 0.001
cfg.noise_string = 'uniform('+str(noise_mean)+','+str(noise_var)+')'
cfg.add_current_stims   = True
cfg.current_stim_amp    = iclamp_amp
cfg.current_stim_start    = 0
cfg.current_stim_duration = cfg.duration
cfg.th_boostFibers = False
# --- Load Spike Times
cfg.th_spikes_file = cfg.Fig_4_1_spikes


#------------------------------------------------------------------------------
# File Name
#------------------------------------------------------------------------------

fileName = 'test0'
cellFigName = 'test0.png'
folderName = 'single_cell_inputs'

#------------------------------------------------------------------------------
# Saving
#------------------------------------------------------------------------------
cfg.filename            = fileName          # Set file output name
cfg.savePickle          = False             # Save params, network and sim output to pickle file
cfg.saveJson            = True              # Save params, network and sim output to JSON file
cfg.saveDataInclude     = ['simData', 'simConfig', 'netParams']#, 'net']
cfg.saveCellConns       = True

cfg.saveFolder          = '../data/init_sims/'+folderName
#------------------------------------------------------------------------------
# Plotting
#------------------------------------------------------------------------------

# cfg.analysis['plotRaster']  = { 'orderInverse': True, 'labels': None, 'markerSize':1,'saveFig':'../data/raster_plot.png','dpi':1000}           
cfg.analysis['plotTraces']  = { 
                                'include': list(range(len(cfg.select_thal_gids))),
                                'timeRange': [0, cfg.duration],
                                'ylim': [-90, 60],
                                'overlay':True,
                                'saveFig':True,
                                }

# cfg.analysis['plotShape']   = { 'saveFig':True}
# cfg.analysis['plot2Dnet'] = {'figSize':(8, 20),'saveFig': 'model_2dnet__.png'}                                                # plot 2D cell positions and connections