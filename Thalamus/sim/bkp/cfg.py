"""
cfg.py 

Contributors: joao.moreira@downstate.edu, salvadordura@gmail.com
"""

from netpyne import specs
import pickle
import NetPyNE_BBP
import numpy as np
import os
saveFigFlag = '' # declared here, but altered afterwards

#------------------------------------------------------------------------------
# Simulation options
#------------------------------------------------------------------------------
cfg                 = specs.SimConfig()     # object of class SimConfig to store simulation configuration
# cfg.skipTime        = 0                  # pre-run simulation for this amount of time 
cfg.skipTime        = 5000                  # pre-run simulation for this amount of time 
cfg.duration        = 5.25*1e3 + cfg.skipTime              # Duration of the simulation, in ms
# cfg.duration        = 0.001*1e3              # Duration of the simulation, in ms
cfg.dt              = 0.025                 # Internal integration timestep to use
# cfg.hParams         = {'v_init': -75, 'celsius': 34} 
# cfg.hParams         = {'v_init': -75.85079, 'celsius': 34} 
cfg.hParams         = {'v_init': -79, 'celsius': 34} 
# cfg.hParams         = {'v_init': -74}       # default VPL cell RMP = -76 mV
cfg.verbose         = False                 # Show detailed messages
cfg.printRunTime    = 0.1
# cfg.connRandomSecFromList = False

################################################################################################################################################
# Adding patches from issues opened on S1 repo (
# https://github.com/BlueBrain/CoreNeuron/issues/843
# https://github.com/suny-downstate-medical-center/S1_netpyne/commit/513c6db577aa574be6d54c1e576cd28b43f4fcdf
# )
cfg.coreneuron = False
cfg.random123 = True
################################################################################################################################################

cfg.recordTraces = {    'V_soma': {'sec': 'soma_0', 'loc': 0.5, 'var': 'v'},
                        # 'V_ptr': {'var': 'ptr'},
                        }

cfg.recordCurrents=True
rec_curr = [('SK_E2','ik'),('TC_HH','ina'),('TC_HH','ik'),('TC_Nap_Et2','ina'),
            ('TC_iA','ik'),('TC_iL','ica'),('TC_iT_Des98','ica'),('TC_ih_Bud97','ih')]
if cfg.recordCurrents:
    for curr in rec_curr: cfg.recordTraces.update({'i__soma_0__'+curr[0]+'__'+curr[1]:{'sec':'soma_0','loc':0.5,'mech':curr[0],'var':curr[1]},})


# cfg.recordStep = 1              # Step size in ms to save data (eg. V traces, LFP, etc)
cfg.recordStep = cfg.dt         # Step size in ms to save data (eg. V traces, LFP, etc)

# cfg.recordStim = False  
# cfg.recordTime = False

#------------------------------------------------------------------------------
# Load BBP circuit
#------------------------------------------------------------------------------

# --- Path to BBP model files

# cfg.base_dir='/Users/joao'
cfg.base_dir=os.path.expanduser("~")

cfg.BBP_rootFolder                  = cfg.base_dir+'/Research/Models/BBP/BBP_thalamus_microcircuit_2'
cfg.sonataConfigFile                = cfg.BBP_rootFolder+'/sonata/circuit_sonata.json'
cfg.morphologyFolder_h5             = cfg.BBP_rootFolder+'/sonata/morphologies_h5'
cfg.morphologyFolder_asc            = cfg.BBP_rootFolder+'/sonata/morphologies/morphologies_asc'
cfg.virtual_input_spikes            = cfg.BBP_rootFolder+'/sonata/simulation_spike_files'
cfg.ct_virtual_noise                = cfg.virtual_input_spikes+'/input_spikes_ct_noise.dat'
cfg.ml_virtual_noise                = cfg.virtual_input_spikes+'/input_spikes_ml_noise.dat'

# --- Path to BBP data files
cfg.BBP_dataFolder                  = cfg.base_dir+'/Research/Models/BBP/BBP_paper_data'
cfg.Fig_4_1_dataFolder              = cfg.BBP_dataFolder+'/Fig4_wakefulness_spontaneous'
cfg.Fig_4_1_spikes                  = cfg.Fig_4_1_dataFolder+'/out.h5'
cfg.Fig_4_1_traces                  = cfg.Fig_4_1_dataFolder+'/soma.bbp.h5'

# --- Path to NetPyNE files
cfg.NetPyNE_rootFolder              = cfg.base_dir+'/Research/Models/BBP/thalamus_netpyne'
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

# # --- Loading Default Config settings
# from defaultCfg import *

# --- Import circuit properties
# cfg.mc_number   = 2
# cfg.mc_num      = 'mc'+str(cfg.mc_number)

cfg.select_microcircuit = None

# --- Convert cell models
cfg.convertCellMorphologies=False

# --- Select Cells
'''
import random;cfg.select_thal_gids=[]; 
import NetPyNE_BBP; circuit_dict    = NetPyNE_BBP.LoadBBPCircuit.getDataFrames ( cfg_file=cfg.sonataConfigFile, microcircuit_number=cfg.mc_number)
for i in range(100):  cfg.select_thal_gids.append(random.choice(list(circuit_dict[cfg.mc_num]['VPL_TC_dAD_ltb'].index)))
for j in range(100):  cfg.select_thal_gids.append(random.choice(list(circuit_dict[cfg.mc_num]['VPL_TC_dNAD_ltb'].index)))

print(' --- select_thal_gids --- \n',cfg.select_thal_gids)
import sys;sys.exit()
# VPL_TC_dAD_ltb:
# [37654, 38204, 35311, 34560, 38647, 35838, 38381, 40674, 38491, 38917, 36707, 39355, 34805, 35342, 40051, 36043, 38911, 37773, 36438, 39040]
# VPL_TC_dNAD_ltb
# [34183, 35877, 42452, 40705, 38808, 40811, 36837, 39026, 38613, 34934, 42180, 42412, 36861, 34578, 40026, 34650, 36844, 42321, 39430, 34441]
# import sys;sys.exit()
'''
# # VPL_TC_dAD_ltb (20) and VPL_TC_dNAD_ltb (20)
# cfg.select_thal_gids = [  35103, 36243, 37082, 33702, 37625, 41429, 35879, 41240, 41615, 34361, 37543, 37177, 41876, 34569, 36963, 41912, 39985, 37055, 36484, 35847,
#                           33774, 41571, 42273, 41996, 38098, 36368, 41395, 37033, 39864, 39123, 36611, 40153, 39451, 35662, 42357, 40624, 40363, 36612, 36499, 33806]

# # # VPL_TC_dAD_ltb (100) and VPL_TC_dNAD_ltb (100)
# cfg.select_thal_gids = [  
#                           35103, 36243, 37082, 33702, 37625, 41429, 35879, 41240, 41615, 34361, 37543, 37177, 41876, 34569, 36963, 41912, 39985, 37055, 36484, 35847, 
#                           33798, 34368, 36219, 39232, 34389, 34242, 38102, 35703, 38487, 41067, 37463, 38468, 36711, 34932, 38346, 34503, 36248, 41454, 36721, 33741, 
#                           40602, 34274, 41534, 33640, 36881, 34859, 36169, 38276, 37409, 34707, 38440, 41237, 38052, 36302, 33602, 41247, 38036, 39429, 38474, 35824, 
#                           38651, 37968, 40213, 42177, 40168, 40215, 41723, 36655, 38134, 41695, 42422, 42460, 36521, 38775, 35220, 35162, 34349, 36440, 35739, 34954, 
#                           37256, 41168, 39751, 38748, 33967, 35343, 40876, 39755, 36185, 41399, 39299, 38971, 37093, 37917, 37599, 34471, 39745, 39477, 42073, 36043, 
#                           41388, 38169, 34773, 34401, 41379, 37475, 38090, 40659, 37782, 38709, 42405, 41353, 41307, 40641, 37685, 39390, 39239, 35684, 34363, 37548, 
#                           34976, 35398, 34977, 34209, 37751, 39276, 38218, 41138, 37435, 37966, 42345, 35864, 34506, 40105, 38470, 34418, 37141, 39362, 33676, 36674, 
#                           36748, 36059, 35158, 40735, 35483, 42198, 34433, 41390, 39229, 40044, 37740, 40122, 36364, 35113, 38793, 40560, 36857, 37553, 41271, 39981, 
#                           41439, 38171, 39183, 41890, 37925, 37824, 38002, 35649, 41579, 38806, 37520, 40430, 33822, 39202, 37863, 41253, 33571, 35332, 35748, 39340, 
#                           33774, 41571, 42273, 41996, 38098, 36368, 41395, 37033, 39864, 39123, 36611, 40153, 39451, 35662, 42357, 40624, 40363, 36612, 36499, 33806
# ]

# VPL_TC_dAD_ltb
# cfg.select_thal_gids = [37654, 38204, 35311, 34560, 38647, 35838, 38381, 40674, 38491, 38917, 36707, 39355, 34805, 35342, 40051, 36043, 38911, 37773, 36438, 39040]
# VPL_TC_dNAD_ltb
# cfg.select_thal_gids = [34183, 35877, 42452, 40705, 38808, 40811, 36837, 39026, 38613, 34934, 42180, 42412, 36861, 34578, 40026, 34650, 36844, 42321, 39430, 34441]

# Rt_RC 100 random gids
# cfg.select_thal_gids = [31571]
# cfg.select_thal_gids = [29876, 29504, 31571, 29227, 31493, 31775, 29925, 29704, 32718, 30882, 30203, 33166, 32301, 29200, 30285, 29500, 30581, 32723, 30643, 29341, 30782, 28753, 33139, 32929, 30827, 28663, 31312, 32315, 32536, 29923, 31316, 29016, 30522, 32613, 28789, 29627, 31361, 29526, 31724, 30625, 33332, 33043, 29545, 32377, 28962, 29628, 29197, 32239, 29367, 28829, 32051, 31055, 31896, 30505, 29659, 32686, 29222, 31655, 29422, 31013, 33094, 30286, 33270, 29171, 30185, 29143, 33013, 31034, 29955, 33134, 29008, 32305, 30295, 30023, 33012, 33022, 29472, 29617, 28647, 32632, 29020, 30149, 29673, 32855, 31476, 33280, 31324, 29375, 30308, 32462, 32348, 33020, 29785, 32680, 30785, 30915, 29113, 30897, 33133, 30848]


# cfg.select_thal_gids = [40408]
# cfg.select_thal_gids = [39690]
# cfg.select_thal_gids = [38906] # 196 spikes

# cfg.select_thal_gids = [35080] # 63 spikes
# cfg.select_thal_gids = [35321] # 63 spikes
# cfg.select_thal_gids = [36636] # 63 spikes
# cfg.select_thal_gids = [37547] # 63 spikes
# cfg.select_thal_gids = [38768] # 63 spikes
# cfg.select_thal_gids = [35080, 35321] # 63 spikes
# cfg.select_thal_gids = [35080, 35321, 36636, 37547, 38768] # 63 spikes

# cfg.select_thal_gids = [40408, 39690]
# cfg.select_thal_gids = [40408, 39690, 35628, 40957, 41864, 37916, 41105, 38048]
# cfg.select_thal_gids = [40408, 39690, 35628, 40957, 41864, 37916, 41105, 38048, 37184, 34832]
# cfg.select_thal_gids = [40408, 39690, 35628, 40957]

# cell_properties = NetPyNE_BBP.LoadBBPCircuit.getNodeProperties(circuit_dict,int(edge_target))

cfg.target = 'VPL_TC'
gids=[36636] ; 
# cfg.hParams['v_init']=-75.85079 ; 
iclamp_amp = 0.0061 # 0.0068 # 0.0072 # 0.008 # 0.01 # 0.03 # 0.061 # 'VPL_TC'     # RMP properties:   0 mA: -76 mV | 0.01 mA: -74.28 mV | 0.02 mA: -72.7 mV | 0.05 mA: -68.5 mV | 0.075 mA: -66.03 mV | 0.08 mA: -65.47 mV | 0.0825 mA: -65.2 mV | 0.084 mA: -65.04 mV | 0.1 mA: -63 mV | 

# cfg.duration = 5100
# gids=[38906] # 'VPL_TC'     

# cfg.target = 'Rt_RC'
# gids=[31571]      ; cfg.hParams['v_init']=-71.4145 ; iclamp_amp = 0.04056 # 0.04054 # 0.040535 # 0.04053 # 0.0405 # 0.038 # 0.041 # 0.045
# # cfg.duration = 5100


# gids = 'test'
# gids = 'test_reduced'
# gids = 'test_single'
# gids = None
cfg.select_thal_gids = NetPyNE_BBP.ConfigCfg.sampleGids(pop=cfg.target, cfg_file=cfg.sonataConfigFile, gids = gids)

# --- Cell Model
cfg.loadCellModel   = False;     cfg.convertCellModel= True;    cfg.saveCellModel   = True; cfg.plotSynLocation = True; cfg.plotCellShape = True

# --- Connectivity Type
cfg.connType        = 'original'    # adds original edge connectivity to the target cell model
# cfg.connType        = 'testSyn'     # adds a single synapse (and single contact) to test EPSP/IPSP amplitude

if cfg.connType == 'original':
    # --- Connectivity
    cfg.th_updateConns      = False     # (True): generates the edge connectivity to the target cell model | (False): Loads from file | Defaults to True if file doesn't exist
    cfg.th_topological      = False      # Redistributes the synapses onto the dendrites according to Driver/IN/RT/CT hierarchy
    cfg.th_useTableVals     = False     # Replace edges values for the values in Table S2 from the paper
    cfg.saveIndividualNoise = True      # Saves individual input noise files for each postsynaptic cell for later plotting
    cfg.th_connOverrides    = True      # Alters the edge properties to match the parameters that were changed before model run, but are not stored in the original model

    cfg.removeChemical      = False
    cfg.removeElectrical    = True      # wont work with TRN inputs as vecstims

    # cfg.th_singleSource    = True
    cfg.th_singleTcPop      = False
    # cfg.th_targetSoma      = True
    # cfg.th_singleConn      = True
    # cfg.modType          = ''
    # cfg.modType          = 'Prob'
    # cfg.modType          = 'Prob_S1'
    cfg.modType          = 'Prob_original'
    # cfg.modType          = 'Det'
    # cfg.modType          = 'Prob_CORENEURON'
    # cfg.modType          = 'Prob_S1_CORENEURON'
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
    # cfg.re_rescale = 1.55
    
    # cfg.cao_secs            = None
    # cfg.rescaleUSE          = None # From BlueConfig file
    
    # cfg.cao_secs            = 1.2
    # cfg.rescaleUSE          = 0.4029343148532312 * cfg.re_rescale # From BlueConfig file
    
    cfg.cao_secs            = 1.2
    cfg.rescaleUSE          = 0.4029343148532312 * cfg.re_rescale # From BlueConfig file


    # --- Change conn weights
    cfg.rescale_conn_weight = { 'MedialLemniscus_projections':  1,
                                'CorticoThalamic_projections':  1,
                                'thalamus_neurons|Rt_RC':       1,
                                'thalamus_neurons|VPL_IN':      1,
                                'thalamus_neurons|VPL_TC':      1,
                                }
    # cfg.rescale_conn_weight = { 'MedialLemniscus_projections':  1.1,
    #                             'CorticoThalamic_projections':  0.9,
    #                             'thalamus_neurons|Rt_RC':       0.9,
    #                             'thalamus_neurons|VPL_IN':      0.9,
    #                             'thalamus_neurons|VPL_TC':      1,
    #                             }
    # cfg.rescale_conn_weight = { 'MedialLemniscus_projections':  1,
    #                             'CorticoThalamic_projections':  0.8,
    #                             'thalamus_neurons|Rt_RC':       0.8,
    #                             'thalamus_neurons|VPL_IN':      0.8,
    #                             'thalamus_neurons|VPL_TC':      1,
    #                             }
    # cfg.rescale_conn_weight = { 'MedialLemniscus_projections':  1,
    #                             'CorticoThalamic_projections':  0.5,
    #                             'thalamus_neurons|Rt_RC':       0.5,
    #                             'thalamus_neurons|VPL_IN':      0.5,
    #                             'thalamus_neurons|VPL_TC':      1,
    #                             }
                        
    # --- Change biophysics
    cfg.removeExtraCurrents=False

    cfg.changeCurrents  = []
    cfg.rescaleCurrents = []
    # cfg.rescaleCurrents = [('TC_iA','gk_max',0.5)]
    # cfg.rescaleCurrents = [('TC_iA','gk_max',0.25)]

    # cfg.rescaleCurrents = [('TC_HH','gna_max',0.017196474865155516*100)]


    # --- Stims
    cfg.add_current_stims_noise   = False
    noise_mean = 0
    noise_var = 0.001
    # noise_var = abs(cfg.hParams['v_init']*0.001)
    cfg.noise_string = 'uniform('+str(noise_mean)+','+str(noise_var)+')'

    cfg.add_current_stims   = True
    cfg.current_stim_amp    = iclamp_amp
    cfg.current_stim_start    = 0
    cfg.current_stim_duration = cfg.duration

    # cfg.current_stim_amp    = 0.061
    # cfg.current_stim_amp    = 0.084
    # cfg.current_stim_start    = 0

    cfg.th_boostFibers = False

    # --- Load Spike Times
    cfg.th_spikes_file = cfg.Fig_4_1_spikes

    if cfg.th_topological:          saveFigFlag+='_topological'
    if cfg.th_useTableVals:         saveFigFlag+='_tableVals'
    if cfg.th_connOverrides:        saveFigFlag+='_connOverr'
    if cfg.cao_secs is not None:    saveFigFlag+='_ca|'+str(cfg.cao_secs)
    if cfg.removeExtraCurrents:     saveFigFlag+='_noIH_noIA'
    if cfg.add_current_stims_noise:       saveFigFlag+='_preIClampNoise'
    if cfg.add_current_stims:       saveFigFlag+='_IClamp_'+str(cfg.current_stim_amp)
    if len(cfg.changeCurrents)>0:   saveFigFlag+='_changedCurr'
    if len(cfg.rescaleCurrents)>0:  saveFigFlag+='_rescaledCurr'
    if np.prod(list(cfg.rescale_conn_weight.values()))!=1:  saveFigFlag+='_rescaledWeights_'+'_'.join(map(str, list(cfg.rescale_conn_weight.values())))
    if cfg.th_boostFibers:          saveFigFlag+='_boostFibers'
    if cfg.th_singleTcPop:          saveFigFlag+='_singleTcPop'
    
    for pathway in cfg.th_select_pathways:
        if 'CorticoThalamic' in pathway: saveFigFlag+='_ct'
        if 'MedialLemniscus' in pathway: saveFigFlag+='_ml'
        if 'thalamus_neurons|Rt_RC' in pathway: saveFigFlag+='_rt'
        if 'thalamus_neurons|VPL_IN' in pathway: saveFigFlag+='_in'
        if 'thalamus_neurons|VPL_TC' in pathway: saveFigFlag+='_tc'
    

elif cfg.connType == 'testSyn':
    # --- Stim source
    cfg.ts_targetMType  = 'VPL_TC'

    cfg.modType      = '' # Options: ('Det','Prob_S1','Prob')
    # cfg.modType      = 'Prob_S1' # Options: ('Det','Prob_S1','Prob')
    cfg.ts_spkFirst     = 3500; cfg.ts_spkNum = 8; cfg.ts_spkInterv = 25
    cfg.ts_spkTimes     = [500]+[cfg.ts_spkFirst+i*cfg.ts_spkInterv for i in range(cfg.ts_spkNum)]        +[5000]

    '''
    Syns per Conn (calculated using <calculate_BBP_conn_properties.py>)
    CorticoThalamic_projections 	(mean, std):  (1.2090532221570451, 0.5171314834679583)
    MedialLemniscus_projections 	(mean, std):  (12.483479074117373, 4.196180573789769)
    thalamus_neurons|Rt_RC 	        (mean, std):  (4.625980301095814, 3.9618020799219265)
    '''
    # cfg.ts_edgeSource=('external__chemical__CorticoThalamic_projections__',   None,   'CorticoThalamic_projections', 1    ) # 1.2090532221570451)
    # cfg.ts_edgeSource=('external__chemical__MedialLemniscus_projections__',   None,   'MedialLemniscus_projections', 12   ) # 12.483479074117373)
    cfg.ts_edgeSource=('internal__chemical__thalamus_neurons|Rt_RC__',        None,   'thalamus_neurons|Rt_RC',      5    ) # 4.625980301095814)
    print('edge_source: ',cfg.ts_edgeSource)



# cfg.ts_edgeSource=('MedialLemniscus_projections',1221,'ML')
# cfg.ts_edgeSource=('thalamus_neurons|Rt_RC',33473,'RT')


#------------------------------------------------------------------------------
# File Name
#------------------------------------------------------------------------------

# --- Selecting type of MOD file to be used
if   cfg.modType == 'Det':          saveFigFlag+='_MOD_S1Det' # S1 Deterministic  implementation of the BBP mod files
elif cfg.modType == 'Prob_S1':      saveFigFlag+='_MOD_S1Prob' # S1 Probabilistic  implementation of the BBP mod files
elif cfg.modType == 'Prob_original':saveFigFlag+='_MOD_ThOriginal' # S1 Probabilistic  implementation of the BBP mod files
else:                               saveFigFlag+='_MOD_ThProb' # Original Thalamus implementation of the BBP mod files


if cfg.connType == 'original':

    if len(cfg.select_thal_gids)==1:
        outputName = cfg.target+'__'+str(cfg.select_thal_gids[0])+saveFigFlag
    else:
        c_name=''
        for cell in cfg.select_thal_gids:c_name+='|'+str(cell)
        outputName=cfg.target+'__'+c_name+saveFigFlag
elif cfg.connType == 'testSyn':
    outputName=cfg.target+'__'+str(len(cfg.select_thal_gids))+'_cells__'+cfg.ts_edgeSource[2]+'|'+str(cfg.ts_edgeSource[3])+saveFigFlag



fileName = outputName
cellFigName = outputName+'.png'

if   cfg.connType == 'original':    folderName = 'single_cell_inputs'
elif cfg.connType == 'testSyn':     folderName = 'single_synapse_inputs'
else:                               folderName = ''

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

# popList = ['pop__'+str(thal_gid) for thal_gid in cfg.select_thal_gids]
# popList = ['VPL_TC__pop']
# popList.sort




# cfg.analysis['plotRaster']  = { 'orderInverse': True, 'labels': None, 'markerSize':1,'saveFig':'../data/raster_plot.png','dpi':1000}           
cfg.analysis['plotTraces']  = { 
                                'include': list(range(len(cfg.select_thal_gids))),
                                # 'include': [0,1],
                                # 'include': [(pop, 0) for pop in popList],
                                # 'include': ['all'],
                                'timeRange': [0, cfg.duration],
                                'ylim': [-90, 60],
                                'overlay':True,
                                'saveFig':'../init_figs/'+cellFigName,
                                # 'saveFig':True,
                                }
# cfg.analysis['plotShape']   = { 'saveFig':True}

# cfg.analysis['plot2Dnet'] = {'figSize':(8, 20),'saveFig': 'model_2dnet__.png'}                                                # plot 2D cell positions and connections