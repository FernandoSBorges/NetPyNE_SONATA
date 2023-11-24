"""
barreloid_cfg.py 

Contributors: joao.moreira@downstate.edu, salvadordura@gmail.com
"""

from netpyne import specs
import pickle
import NetPyNE_BBP
import json
import sys
import os

saveFigFlag = '' # declared here, but altered afterwards

#------------------------------------------------------------------------------
#   Simulation options
#------------------------------------------------------------------------------
cfg       = specs.SimConfig()     # object of class SimConfig to store simulation configuration
# cfg.duration        = 0.001*1e3 # (ms)
cfg.duration        = 8000 # (ms)
# cfg.duration        = 3000 # (ms)
# cfg.duration        = 2000 # (ms)
# cfg.duration        = 1 # (ms)
# cfg.duration        = 100 # (ms)
# cfg.duration        = 500 # (ms)
# cfg.duration        = 800 # (ms)
# cfg.duration        = 1000 # (ms)
cfg.dt              = 0.025                 # Internal integration timestep to use
cfg.printRunTime    = 0.01 # (s)

cfg.hParams         = {'v_init': -60, 'celsius': 34} # changing to -60 mV to remove initial bursts in "in vivo"-like sim
# cfg.hParams         = {'v_init': -79, 'celsius': 34} # bbp params
# cfg.hParams         = {'v_init': -65} 
# cfg.hParams         = {'v_init': -74}       # default VPL cell RMP = -76 mV
cfg.verbose         = False                 # Show detailed messages
# cfg.cvode_atol      = 1e-6
# cfg.cache_efficient = True
# cfg.connRandomSecFromList = False
# cfg.random123 = True

#------------------------------------------------------------------------------
# Recording
#------------------------------------------------------------------------------
cfg.recordStep = 0.1              # Step size in ms to save data (eg. V traces, LFP, etc)
# cfg.recordStep = cfg.dt         # Step size in ms to save data (eg. V traces, LFP, etc)

cfg.recordTraces = {    'V_soma': {'sec': 'soma_0', 'loc': 0.5, 'var': 'v'},
                        # 'V_ptr': {'var': 'ptr'},
                        }
cfg.recordCurrents=False
rec_curr = [('SK_E2','ik'),('TC_HH','ina'),('TC_HH','ik'),('TC_Nap_Et2','ina'),
            ('TC_iA','ik'),('TC_iL','ica'),('TC_iT_Des98','ica'),('TC_ih_Bud97','ih')]
if cfg.recordCurrents:
    for curr in rec_curr: cfg.recordTraces.update({'i__soma_0__'+curr[0]+'__'+curr[1]:{'sec':'soma_0','loc':0.5,'mech':curr[0],'var':curr[1]},})

#------------------------------------------------------------------------------
# Path configuration
#------------------------------------------------------------------------------
cfg.convertCellMorphologies = False
cfg.loadCellModel           = False
cfg.convertCellModel        = True
cfg.saveCellModel           = True
cfg.plotSynLocation         = True

# cfg.base_dir = '/Users/joao'
cfg.base_dir=os.path.expanduser("~")

cfg.BBP_rootFolder                  = cfg.base_dir+'/Research/Models/BBP/BBP_thalamus_microcircuit_2'
cfg.sonataConfigFile                = cfg.BBP_rootFolder+'/sonata/circuit_sonata.json'
cfg.morphologyFolder_h5             = cfg.BBP_rootFolder+'/sonata/morphologies_h5'

cfg.NetPyNE_rootFolder              = cfg.base_dir+'/Research/Models/BBP/thalamus_netpyne'
cfg.NetPyNE_JSON_cells              = cfg.NetPyNE_rootFolder+'/cells/netpyne_morphologies'
cfg.NetPyNE_templateCells           = cfg.NetPyNE_rootFolder+'/mod'
cfg.NetPyNE_exportedCells           = cfg.NetPyNE_rootFolder+'/cells/morphologies_swc'
cfg.NetPyNE_L6A_JSON_cells          = cfg.NetPyNE_rootFolder+'/cells/S1_BBP_cells/'
cfg.NetPyNE_network_template        = cfg.NetPyNE_rootFolder+'/conn/barreloid_network_template/network_template.json'

cfg.prepare_HPC_mode = False # prepares the project to run in HPC environments
cfg.run_HPC_mode     = True # configure the script to run in HPC environments
if cfg.prepare_HPC_mode or cfg.run_HPC_mode: cfg.NetPyNE_JSON_cells+='_barreloid'

cfg.loadCircuitProperties           = True
cfg.saveCircuitProperties           = True
cfg.stored_circuit_path             = cfg.NetPyNE_rootFolder+'/conn/bbp_circuit_propeties/circuit_dict.pkl'

cfg.BBP_conn_properties             = cfg.NetPyNE_rootFolder+'/conn/calculate_BBP_conn_properties/BBP_conn_propeties.json'

#------------------------------------------------------------------------------
# Simulation Configuration
#------------------------------------------------------------------------------

cfg.modType                         = 'Prob_original'
# cfg.modType                         = 'Prob_S1'

# --- Loading the pre-generated network connectivity

# cfg.select_thal_gids=[37654, 38204, 35311, 34560, 38647, 35838, 38381, 40674, 38491, 38917, 36707, 39355, 34805, 35342, 40051, 36043, 38911, 37773, 36438, 39040]

# cfg.network_template = 'network_template.json'
cfg.network_template = cfg.NetPyNE_network_template
NetPyNE_BBP.Prompt.headerMsg('Sourcing network template from \t\t'+cfg.network_template)
try:    f = open(cfg.network_template)
except: sys.exit('Loading network failed')
cfg.store_network = json.load(f)

cfg.load_vpm_gids = len(cfg.store_network['pops']['VPM__pop']['cellGids'])
cfg.load_trn_gids = len(cfg.store_network['pops']['TRN__pop']['cellGids'])

# cfg.randomizeNetwork=True
cfg.randomizeNetwork= False

cfg.simplify_gids   = True # simplifying the number of morphologies added

if cfg.randomizeNetwork:
    NetPyNE_BBP.Prompt.headerMsg('Randomizing network gids')
    cfg.VPM_gids = NetPyNE_BBP.ConfigCfg.sampleGids(pop='VPL_TC', cfg_file=cfg.sonataConfigFile, gids = cfg.load_vpm_gids)
    cfg.TRN_gids = NetPyNE_BBP.ConfigCfg.sampleGids(pop='Rt_RC',  cfg_file=cfg.sonataConfigFile, gids = cfg.load_trn_gids)
else:
    NetPyNE_BBP.Prompt.headerMsg('Loading preset network gids')
    if cfg.simplify_gids:
        print('Simplifying network morphologies')
        # VPM (318) GIDs picked from: [35050]
        cfg.VPM_gids=[35050 for i in range(318)]
        # TRN (101) GIDs picked from: [32411]
        cfg.TRN_gids=[32411 for i in range(101)]

        # # VPM (318) GIDs picked from: [40352, 34066, 37812, 38737, 35050]
        # cfg.VPM_gids=[35050, 37812, 38737, 37812, 34066, 38737, 38737, 40352, 37812, 34066, 34066, 40352, 35050, 40352, 34066, 34066, 40352, 35050, 37812, 38737, 34066, 34066, 40352, 38737, 37812, 34066, 37812, 37812, 34066, 38737, 35050, 34066, 34066, 34066, 40352, 38737, 37812, 40352, 40352, 38737, 34066, 38737, 37812, 35050, 35050, 37812, 35050, 34066, 40352, 38737, 35050, 35050, 34066, 40352, 37812, 37812, 34066, 35050, 40352, 38737, 35050, 35050, 37812, 34066, 37812, 38737, 40352, 38737, 35050, 38737, 37812, 34066, 37812, 34066, 40352, 35050, 38737, 34066, 35050, 40352, 35050, 34066, 38737, 38737, 40352, 40352, 38737, 34066, 40352, 35050, 34066, 35050, 38737, 35050, 40352, 34066, 40352, 37812, 34066, 35050, 38737, 40352, 34066, 38737, 38737, 35050, 40352, 40352, 34066, 37812, 40352, 34066, 38737, 38737, 34066, 40352, 34066, 34066, 38737, 35050, 37812, 37812, 38737, 35050, 38737, 40352, 37812, 34066, 37812, 34066, 34066, 35050, 34066, 40352, 34066, 40352, 40352, 37812, 37812, 34066, 40352, 34066, 38737, 38737, 40352, 34066, 40352, 34066, 40352, 34066, 38737, 40352, 40352, 34066, 40352, 34066, 37812, 37812, 38737, 40352, 34066, 35050, 35050, 35050, 37812, 38737, 38737, 38737, 35050, 40352, 38737, 37812, 34066, 40352, 35050, 37812, 37812, 34066, 34066, 34066, 40352, 38737, 40352, 40352, 37812, 40352, 40352, 37812, 35050, 40352, 34066, 34066, 35050, 38737, 37812, 34066, 37812, 35050, 35050, 37812, 38737, 34066, 38737, 35050, 40352, 40352, 37812, 38737, 35050, 35050, 38737, 40352, 37812, 38737, 37812, 40352, 37812, 35050, 38737, 40352, 34066, 37812, 38737, 38737, 35050, 35050, 37812, 40352, 35050, 35050, 40352, 40352, 35050, 35050, 35050, 34066, 34066, 37812, 38737, 40352, 40352, 34066, 35050, 40352, 35050, 37812, 35050, 38737, 38737, 35050, 38737, 38737, 37812, 40352, 37812, 35050, 35050, 35050, 37812, 40352, 40352, 37812, 35050, 35050, 34066, 38737, 35050, 38737, 37812, 34066, 37812, 37812, 38737, 37812, 34066, 37812, 38737, 34066, 38737, 34066, 38737, 34066, 35050, 40352, 35050, 40352, 37812, 40352, 40352, 37812, 34066, 34066, 34066, 35050, 38737, 34066, 38737, 37812, 38737, 37812, 34066, 38737, 35050, 40352, 34066, 40352, 37812, 38737, 40352, 35050, 40352, 34066, 37812, 37812, 38737, 34066, 38737, 37812]
        # # TRN (101) GIDs picked from: [30404, 32411, 31136, 32952, 29223]
        # cfg.TRN_gids=[32411, 31136, 32411, 31136, 30404, 29223, 32952, 31136, 32411, 32411, 29223, 32952, 30404, 32952, 31136, 29223, 29223, 31136, 29223, 29223, 29223, 32411, 32952, 32952, 31136, 32952, 32952, 32411, 32952, 32952, 29223, 29223, 30404, 29223, 31136, 30404, 32411, 32952, 32952, 30404, 32411, 31136, 32411, 32952, 32952, 32952, 30404, 30404, 29223, 32952, 32411, 29223, 32411, 29223, 30404, 29223, 29223, 32411, 32411, 32411, 32411, 29223, 30404, 30404, 32952, 32952, 29223, 29223, 30404, 32411, 32411, 30404, 30404, 30404, 31136, 30404, 32411, 32411, 32952, 32411, 32411, 29223, 30404, 30404, 30404, 29223, 32411, 31136, 29223, 32952, 29223, 32952, 32952, 32411, 32952, 32411, 32411, 32411, 31136, 31136, 29223]
    else:
        cfg.VPM_gids=[40352, 34066, 36777, 40166, 37812, 38737, 40296, 39900, 35050, 36666, 42272, 38785, 40826, 42023, 33933, 34466, 38284, 34844, 40465, 39867, 41216, 42326, 39167, 33887, 33559, 36658, 36388, 38608, 36299, 39299, 39251, 41563, 39633, 37330, 34195, 40932, 41549, 37891, 36698, 37851, 37428, 36855, 37478, 40842, 35711, 36579, 40249, 37405, 37301, 42351, 34692, 36310, 35985, 41386, 36520, 41474, 42158, 35215, 42321, 40463, 38254, 41946, 40128, 40039, 40816, 35193, 41190, 40600, 40850, 40988, 39574, 34574, 35455, 36670, 35877, 39030, 40488, 33700, 40106, 37496, 37444, 37426, 34747, 37101, 40674, 35132, 38297, 37758, 41972, 36705, 37135, 38593, 36018, 39293, 34664, 37342, 41990, 41412, 39180, 34903, 40977, 34286, 36937, 37325, 40079, 40269, 34673, 33834, 41329, 35825, 37624, 41490, 37614, 42178, 33677, 36258, 36395, 35841, 34269, 37281, 34671, 39485, 35312, 34693, 37944, 36934, 40523, 38097, 39204, 42120, 41049, 33592, 34605, 38808, 40746, 37217, 34309, 36336, 33772, 37477, 37649, 35219, 39933, 35387, 42283, 41852, 40970, 33988, 33758, 40845, 33689, 40199, 36333, 41191, 34732, 42423, 33737, 40843, 38058, 38726, 34346, 38513, 39488, 40697, 33600, 37144, 38973, 38386, 38718, 40777, 39544, 40444, 34787, 39026, 34423, 37580, 40068, 33920, 36134, 38708, 40036, 38719, 39155, 33905, 37264, 37151, 40453, 41337, 37789, 33562, 39932, 39922, 35064, 38248, 35647, 37979, 40020, 34817, 35795, 36598, 33985, 35998, 42276, 39691, 41416, 36638, 36443, 41706, 41254, 41863, 36222, 36365, 34493, 40836, 33862, 35899, 40611, 37020, 34949, 34376, 34504, 40345, 37422, 35177, 39764, 37443, 36304, 34710, 39791, 34897, 36422, 38384, 38540, 40865, 41027, 36096, 36713, 39218, 39964, 34763, 38743, 39782, 41319, 40184, 42079, 34306, 41217, 41135, 34665, 39008, 34956, 39088, 37048, 41379, 39607, 35124, 41322, 40979, 37712, 36709, 34047, 42106, 36603, 37607, 38261, 35941, 35797, 35232, 34896, 39099, 39425, 39383, 36036, 35978, 39064, 37567, 39386, 40720, 38734, 36073, 42100, 33765, 41042, 37268, 38324, 42065, 37046, 40622, 35461, 38062, 38757, 36794, 38525, 34648, 37965, 36058, 39154, 38321, 39955, 36947, 34892, 36176, 41843, 41270, 35459, 40379, 36575, 41413, 40712, 41249, 34417, 34447, 39220, 37459, 34960, 38394, 34801, 39512]
        cfg.TRN_gids=[31364, 30106, 29095, 30404, 32411, 32689, 31136, 32952, 29223, 29226, 31158, 30716, 30874, 29694, 32244, 28976, 31159, 29987, 29060, 28963, 31917, 33413, 29656, 28828, 30883, 33262, 28657, 31312, 30553, 32956, 33496, 31127, 30686, 32101, 28629, 29205, 29588, 32390, 31050, 31497, 32746, 31493, 32515, 30601, 30161, 32655, 29546, 28867, 31388, 30041, 30876, 33378, 33405, 29303, 32878, 29152, 28751, 30860, 30265, 32556, 33156, 32386, 30801, 29146, 29527, 32576, 29159, 32174, 29217, 29736, 31557, 32851, 29401, 31530, 31839, 30846, 30879, 31170, 32483, 31756, 29935, 29172, 30493, 29094, 33437, 30013, 30440, 28812, 30307, 29266, 31156, 33485, 30093, 30987, 32245, 31512, 29027, 31093, 32087, 29125, 32921]
    if cfg.prepare_HPC_mode: cfg.loadCellModel = False
    elif cfg.run_HPC_mode:   cfg.loadCellModel = True
    else:
        cfg.loadCellModel           = True
        cfg.convertCellMorphologies = False
cfg.select_thal_gids=cfg.VPM_gids+cfg.TRN_gids
cfg.gids_lists = [cfg.VPM_gids,cfg.TRN_gids]

#------------------------------------------------------------------------------
#   Network
#------------------------------------------------------------------------------
cfg.center_point = 500

cfg.re_rescale = 1
# cfg.re_rescale = 1.5
cfg.cao_secs            = 1.2
# cfg.rescaleUSE          = None
cfg.rescaleUSE          = 0.4029343148532312 * cfg.re_rescale # From BlueConfig file

#------------------------------------------------------------------------------
#   MLe inputs
#------------------------------------------------------------------------------
from GenerateStimDataset import SampleData

cfg.deflection_events = [[2000,2150],[3000,3150],[4000,4250],[5000,5250],[6000,6350],[7000,7350]]

# cfg.samplingInterv = 25     # us
cfg.samplingInterv = 1000   # us

# cfg.target_angle = None
# cfg.target_angle = 0
# cfg.target_angle = 45
# cfg.target_angle = 90
# cfg.target_angle = 135
# cfg.target_angle = 180
# cfg.target_angle = 225
# cfg.target_angle = 270
cfg.target_angle = 315

if cfg.target_angle is not None:    
    cfg.deflection_times = cfg.deflection_events+[[10000,10020]]
    cfg.stims_string ='_stims|' + "|".join(str(sublist[0])+'@'+str(sublist[1]) for sublist in cfg.deflection_times)
else:                               
    cfg.deflection_times = [[10000,10020]]
    cfg.stims_string = '_stims|None'

cfg.stim_source_code_folder    = '/stims/Minnery_2003/fig1'
cfg.deflection_model_folder    = cfg.NetPyNE_rootFolder+cfg.stim_source_code_folder+'/deflection_model/'   # Folder to save the output deflection model and sampled raster
# --- Load deflection model and sampled raster (stored in JSON for readability - small files)
cfg.deflection_dict_path = SampleData.LoadJSON(file_path=   cfg.deflection_model_folder+'deflectionModel_'+str(cfg.target_angle)+'|deg'+'_samplingBins|'+str(cfg.samplingInterv)+'us'+cfg.stims_string+'.json')
cfg.deflection_dataset_path =                cfg.deflection_model_folder+'spike_dicts/mleSpikes_deflectionAngle|'+str(cfg.target_angle)+'_samplingBins|'+str(cfg.samplingInterv)+'us'+cfg.stims_string+'_simplifiedDataset'+'.pkl'

#------------------------------------------------------------------------------
#   Stimulation
#------------------------------------------------------------------------------

cfg.add_current_stims       = True
if cfg.add_current_stims:
    cfg.current_stim_targets    = ['VPM__pop']
    # cfg.current_stim_targets    = ['VPM__pop','TRN__pop']
    cfg.current_stim_amp        = 0.1
    cfg.current_stim_start      = 0
    cfg.current_stim_duration   = cfg.duration

cfg.add_bkg_stim = True
if cfg.add_bkg_stim:
    cfg.bkg_rate    = [40,          200,            200,            200,]
    cfg.bkg_noise   = [1,           1,              1,              1,]
    cfg.bkg_weight  = [0.001,       0.0005,         0.0005,         0.0005,]
    cfg.bkg_delay   = [0,           0,              0,              0,]
    cfg.bkg_synMech = ['exc',       'exc',          'exc',          'exc',]
    cfg.bkg_pop     = ['L6A__cell', 'L6A__cell',    'VPM__cell',    'TRN__cell',] # 'cellType'
    
    # cfg.bkg_CTX_rate    = [40,      200     ]
    # cfg.bkg_CTX_noise   = [1,       1       ]
    # cfg.bkg_CTX_weight  = [0.001,   0.0005  ]
    # cfg.bkg_CTX_delay   = [0,       0       ]
    # cfg.bkg_CTX_synMech = ['exc',   'exc'   ]

#------------------------------------------------------------------------------
#   Connectivity
#------------------------------------------------------------------------------

def get_synsPerConn(conn_data, pre_pop, post_pop, conn_type):
    try:    
        syns_per_conn = round(conn_data[post_pop][conn_type][pre_pop]['synsPerConn'][0])
        if (conn_data[post_pop][conn_type][pre_pop]['synsPerConn'][0]>0) and int(conn_data[post_pop][conn_type][pre_pop]['synsPerConn'][0])==0: syns_per_conn=1
    except: syns_per_conn = 1
    return syns_per_conn

# --- Connectivity properties
import Build_Net as BN
# --- Loads the conn data calculated from the projections in the BBP somatosensory thalamus model
cfg.conn_data = BN.NetworkConversion.convertConnProperties(filePath=cfg.BBP_conn_properties)

# --- Enabling Connectivity
cfg.conn_MLe_VPM    = True     # works
cfg.conn_TRN_VPM    = True     # fail still - works now
cfg.conn_TRN_TRN    = True     # works
cfg.conn_VPM_TRN    = True     # works

cfg.conn_TRNe_TRNe  = True     # 2023_11_1: BREAKS THE MPI SIMULATION - works now with modified gap.mod with NET_RECEIVE block

cfg.conn_L6A_VPM    = True     # fail still - works now
cfg.conn_L6A_TRN    = True     # fail still - works now
cfg.conn_VPM_L6A    = True     # works


#------------------------------------------------------------------------------
# scale_weights = 1.75    # 000008_8000ms_allWeights1.75x_270degStim_adjustedProjections_adjustedStims_noTRNIclamp_singleCellTemplate
# scale_weights = 2.5     # 000009_8000ms_allWeights2.5x_270degStim_adjustedProjections_adjustedStims_noTRNIclamp_singleCellTemplate
scale_weights = 1.25     # 000010_8000ms_allWeights1.25x_90degStim_adjustedProjections_adjustedStims_noTRNIclamp_singleCellTemplate

#------------------------------------------------------------------------------
# --- MLe => VPM
''' (Takeuchi, 2017) # of boutons / mle fiber:   40.50919230368733 '''    
cfg.MLe_VPM__syns_per_conn  = 40
cfg.MLe_VPM__synMech        = 'syn|MLe|VPM|exc'
cfg.MLe_VPM__conn_type      = 'chem'
cfg.MLe_VPM__secList        = 'inputs__proximal'        # (Jones, 2002) - driver input targets proximal sections
cfg.MLe_VPM__weight         = 10
#------------------------------------------------------------------------------
# --- L6A => VPM - distanceBasedProbability_1D_exponential
cfg.L6A_VPM__conn_prob      = 2.0827527202658684
cfg.L6A_VPM__pre_pop        = 'L6A'
cfg.L6A_VPM__post_pop       = 'VPM'
cfg.L6A_VPM__pre_pop_axis   = 'x'
cfg.L6A_VPM__post_pop_axis  = 'y'
cfg.L6A_VPM__conn_type      = 'chem'
cfg.L6A_VPM__syn_mech       = 'syn|L6A|VPM|exc'
cfg.L6A_VPM__target_secs    = 'inputs__distal'
cfg.L6A_VPM__weight         = scale_weights
cfg.L6A_VPM__syns_per_conn  = get_synsPerConn(conn_data=cfg.conn_data,pre_pop=cfg.L6A_VPM__pre_pop,post_pop=cfg.L6A_VPM__post_pop,conn_type=cfg.L6A_VPM__conn_type)

#------------------------------------------------------------------------------
# --- L6A => TRN - distanceBasedProbability_1D_exponential
cfg.L6A_TRN__conn_prob      = 0.8010587385637955
cfg.L6A_TRN__pre_pop        = 'L6A'
cfg.L6A_TRN__post_pop       = 'TRN'
cfg.L6A_TRN__pre_pop_axis   = 'x'
cfg.L6A_TRN__post_pop_axis  = 'y'   # (2023_11_04) - testing adjustment of projections
# cfg.L6A_TRN__post_pop_axis  = 'x'
cfg.L6A_TRN__conn_type      = 'chem'
cfg.L6A_TRN__syn_mech       = 'syn|L6A|TRN|exc'
cfg.L6A_TRN__target_secs    = 'inputs__distal'
cfg.L6A_TRN__weight         = scale_weights
cfg.L6A_TRN__syns_per_conn  = get_synsPerConn(conn_data=cfg.conn_data,pre_pop=cfg.L6A_TRN__pre_pop,post_pop=cfg.L6A_TRN__post_pop,conn_type=cfg.L6A_TRN__conn_type)

#------------------------------------------------------------------------------
# --- TRN => VPM
cfg.TRN_VPM__ref_pop        = 'MLe'
cfg.TRN_VPM__pre_pop        = 'TRN'
cfg.TRN_VPM__post_pop       = 'VPM'
cfg.TRN_VPM__conn_type      = 'chem'
cfg.TRN_VPM__syn_mech       = 'syn|TRN|VPM|inh'
cfg.TRN_VPM__target_secs    = 'inputs__intermediate'
cfg.TRN_VPM__weight         = scale_weights
cfg.TRN_VPM__syns_per_conn  = get_synsPerConn(conn_data=cfg.conn_data,pre_pop=cfg.TRN_VPM__pre_pop,post_pop=cfg.TRN_VPM__post_pop,conn_type=cfg.TRN_VPM__conn_type)
# --- Number of conns from MLe to VPM, used as a parameter to estimate the TRN to VPM convergence
# conns_MLe_to_VPM = 16.67924528301887 # must be updated if values change
conns_MLe_to_VPM = 55.59748427672956 # must be updated if values change

''' # Ratio of (TRN-to->VPM)/(MLe-to->VPM) '''
# # BBP ratio           = 454.1682034044653/122.97733554449974  = 3.6931049237127365
ratio_TRNtoVPM_MLetoVPM = (cfg.conn_data[cfg.TRN_VPM__post_pop][cfg.TRN_VPM__conn_type][cfg.TRN_VPM__pre_pop]['synsPerConn'][0]*cfg.conn_data[cfg.TRN_VPM__post_pop][cfg.TRN_VPM__conn_type][cfg.TRN_VPM__pre_pop]['convergence'][0])/(cfg.conn_data[cfg.TRN_VPM__post_pop][cfg.TRN_VPM__conn_type][cfg.TRN_VPM__ref_pop]['synsPerConn'][0]*cfg.conn_data[cfg.TRN_VPM__post_pop][cfg.TRN_VPM__conn_type][cfg.TRN_VPM__ref_pop]['convergence'][0])
# (çavdar, 2007) ratio  = 48/28                                 = 1.7142857142857142    # --- Testing this value instead to check if BBP (TRN->VPM) conns are not too high
# ratio_TRNtoVPM_MLetoVPM = 48/28

'''
conns = convergence_ * syns_per_conn_
convergence_TRN-to->VPM = ((TRN-to->VPM)/(MLe-to->VPM))*(conns_MLe_to_VPM/syns_per_conn_TRN_to_VPM)
'''
cfg.TRN_VPM__conn_method = 'convergence'
cfg.TRN_VPM__conn_rule   = (ratio_TRNtoVPM_MLetoVPM)*(conns_MLe_to_VPM/cfg.TRN_VPM__syns_per_conn)

#------------------------------------------------------------------------------
# --- TRN => TRN
cfg.TRN_TRN__pre_pop        = 'TRN'
cfg.TRN_TRN__post_pop       = 'TRN'
cfg.TRN_TRN__conn_type      = 'chem'
cfg.TRN_TRN__syn_mech       = 'syn|TRN|TRN|inh'
cfg.TRN_TRN__target_secs    = 'inputs__distal'
cfg.TRN_TRN__syns_per_conn  = get_synsPerConn(conn_data=cfg.conn_data,pre_pop=cfg.TRN_TRN__pre_pop,post_pop=cfg.TRN_TRN__post_pop,conn_type=cfg.TRN_TRN__conn_type)
cfg.TRN_TRN__conn_method    = 'convergence'
cfg.TRN_TRN__weight         = scale_weights
cfg.TRN_TRN__conn_rule      = cfg.conn_data[cfg.TRN_TRN__post_pop][cfg.TRN_TRN__conn_type][cfg.TRN_TRN__pre_pop]['convergence'][0]
print('TRN TRN chemical: ', cfg.TRN_TRN__conn_rule)

#------------------------------------------------------------------------------
# --- TRNe => TRNe
cfg.TRNe_TRNe__pre_pop      = 'TRN'
cfg.TRNe_TRNe__post_pop     = 'TRN'
cfg.TRNe_TRNe__conn_type    = 'elec'
cfg.TRNe_TRNe__syn_mech     = 'gap_nr' # Gap junction mechanism from https://github.com/BlueBrain/neuron_reduce/blob/08bea3520c0f535cdba27ef0c3e4d8f970d08604/tests/TestsFiles/Test_4_LBC_amsalem/mod/halfgap.mod#L4 
cfg.TRNe_TRNe__target_secs  = 'inputs__distal'
cfg.TRNe_TRNe__weight       = 0.001
cfg.TRNe_TRNe__syns_per_conn  = get_synsPerConn(conn_data=cfg.conn_data,pre_pop=cfg.TRNe_TRNe__pre_pop,post_pop=cfg.TRNe_TRNe__post_pop,conn_type=cfg.TRNe_TRNe__conn_type)
cfg.TRNe_TRNe__conn_method  = 'convergence'
rescale_gap_convergenge     = 0.1
cfg.TRNe_TRNe__conn_rule    = rescale_gap_convergenge * cfg.conn_data[cfg.TRNe_TRNe__post_pop][cfg.TRNe_TRNe__conn_type][cfg.TRNe_TRNe__pre_pop]['convergence'][0]
print('TRN TRN electrical: ', cfg.TRNe_TRNe__conn_rule)
print('     ----------- >   Using GAP junction ', cfg.TRNe_TRNe__syn_mech)

# cfg.TRNe_TRNe__syn_mech     = 'gap' # Gap junction mechanism from Iavarone, 2023
# cfg.TRNe_TRNe__syn_mech     = 'gap2' # Gap junction mechanism from https://github.com/BlueBrain/neuron_reduce/blob/08bea3520c0f535cdba27ef0c3e4d8f970d08604/tests/TestsFiles/Test_4_LBC_amsalem/mod/halfgap.mod#L4 
# cfg.TRNe_TRNe__syn_mech     = 'esyn'  # Gap junction mechanism from NetPyNE repo

#------------------------------------------------------------------------------
# --- VPM => TRN - distanceBasedProbability_1D_exponential
cfg.VPM_TRN__conn_prob      = 0.48063524313827727
cfg.VPM_TRN__pre_pop        = 'VPM'
cfg.VPM_TRN__post_pop       = 'TRN'
cfg.VPM_TRN__pre_pop_axis   = 'y'
cfg.VPM_TRN__post_pop_axis  = 'y'   # (2023_11_04) - testing adjustment of projections
# cfg.VPM_TRN__post_pop_axis  = 'x'
cfg.VPM_TRN__conn_type      = 'chem'
cfg.VPM_TRN__syn_mech       = 'syn|VPM|TRN|exc'
cfg.VPM_TRN__target_secs    = 'inputs__proximal'
cfg.VPM_TRN__weight         = scale_weights
cfg.VPM_TRN__syns_per_conn  = get_synsPerConn(conn_data=cfg.conn_data,pre_pop=cfg.VPM_TRN__pre_pop,post_pop=cfg.VPM_TRN__post_pop,conn_type=cfg.VPM_TRN__conn_type)

#------------------------------------------------------------------------------
# --- VPM => L6A - distanceBasedProbability_1D_exponential
cfg.simplifyL6A=False

# cfg.VPM_L6A__conn_prob      = 0.075
# cfg.VPM_L6A__coef           = [-185.45820847,  228.16359477] # from (Meyer,2010) m and b terms from y=(mx+b) x=depth/y=#_of_boutons

cfg.VPM_L6A__conn_prob      = 0.2 # --- TEST: 2023_11_09 - removing linear scaling of conn prob because the axis don't match - using distance-based conn instead

cfg.VPM_L6A__pre_pop        = 'VPM'
cfg.VPM_L6A__post_pop       = 'L6A'
cfg.VPM_L6A__pre_pop_axis   = 'y'
cfg.VPM_L6A__post_pop_axis  = 'x'   # (2023_11_04) - testing adjustment of projections
# cfg.VPM_L6A__post_pop_axis  = 'y'
cfg.VPM_L6A__conn_type      = 'chem'
cfg.VPM_L6A__syn_mech       = 'syn|VPM|L6A|exc'
cfg.VPM_L6A__weight         = scale_weights
cfg.VPM_L6A__TC_syns_per_conn = 9        # from S1-netpyne model

CT_cells            = 666
CT_bouton_per_cell  = 133      # from (Meyer, 2010) 
TC_cells            = 318
TC_CT_conns = (CT_cells*CT_bouton_per_cell)/(TC_cells*cfg.VPM_L6A__TC_syns_per_conn)
print('TC cell target convergence onto CT neurons: ',round(TC_CT_conns))

#------------------------------------------------------------------------------
# Saving
#------------------------------------------------------------------------------

# cfg.sim_tag = '000008_8000ms_allWeights1.75x_noStim_adjustedProjections_adjustedStims_noTRNIclamp_singleCellTemplate'

# cfg.sim_tag = '000009_8000ms_allWeights2.5x_90degStim_adjustedProjections_adjustedStims_noTRNIclamp_singleCellTemplate'

# cfg.sim_tag = '000010_8000ms_allWeights1.25x_270degStim_adjustedProjections_adjustedStims_noTRNIclamp_singleCellTemplate'
# cfg.sim_tag = '000010_8000ms_allWeights1.25x_noStim_adjustedProjections_adjustedStims_noTRNIclamp_singleCellTemplate'

# cfg.sim_tag = '000011_8000ms_allWeights1.25x_noStim_adjustedProjections_adjustedStims_noTRNIclamp_singleCellTemplate'
# cfg.sim_tag = '000011_8000ms_allWeights1.25x_225degStim_adjustedProjections_adjustedStims_noTRNIclamp_singleCellTemplate'

# cfg.sim_tag = '000012_8000ms_allWeights1.25x_noStim_adjustedProjections_adjustedStims_noTRNIclamp_singleCellTemplate'
# cfg.sim_tag = '000012_8000ms_allWeights1.25x_315degStim_adjustedProjections_adjustedStims_noTRNIclamp_singleCellTemplate'

# cfg.sim_tag = '000013_8000ms_noStim'
cfg.sim_tag = '000013_8000ms_315degStim'

if cfg.simplify_gids:           cfg.sim_tag+='_simplifiedGids'
if cfg.rescaleUSE is not None:  cfg.sim_tag+='_USE_'+str(cfg.rescaleUSE)

cfg.filename            = 'barr_net_'+cfg.sim_tag          # Set file output name
cfg.savePickle          = True             # Save params, network and sim output to pickle file
cfg.saveJson            = False              # Save params, network and sim output to JSON file
cfg.saveDataInclude     = ['simData', 'simConfig', 'netParams', 'net']
cfg.saveCellConns       = False
NetPyNE_BBP.Prompt.headerMsg('THIS IS A DEVELOPMENT/DEBUG SESSION: CELL CONNS ARE NOT BEING SAVED - for final runs set cfg.saveCellConns = True')
folderName              = 'barreloid_network'

cfg.saveFolder          = '../data/barreloid_sims/'+folderName

#------------------------------------------------------------------------------
# --- Plotting
#------------------------------------------------------------------------------
cfg.saveFigPath = cfg.NetPyNE_rootFolder+'/figs/barreloid_figs'
# cfg.analysis['plot2Dnet']   = {'figSize':(8, 20),
#                                'saveFig': cfg.saveFigPath+'/'+cfg.filename+'_2Dnet'+'.png'} # plot 2D cell positions and connections

# cfg.analysis['plotConn']    = {'figSize':(10, 10), 
#                                'includePre':['VPM__pop', 'TRN__pop', 'L6A__pop'],
#                                'includePost':['VPM__pop', 'TRN__pop', 'L6A__pop'],
#                                'saveFig': cfg.saveFigPath+'/'+cfg.filename+'_connMatrix'+'.png'} # plot connectivity matrix

cfg.analysis['plotRaster']  = {'figSize':(25, 20), 
                            #    'orderBy': 'y',
                               'saveFig': cfg.saveFigPath+'/'+cfg.filename+'raster'+'.png'} # Plot a raster

pops = ['VPM__pop', 'TRN__pop', 'L6A__pop']

import numpy as np
cfg.analysis['plotTraces']  = {
                                'include': [(pop,list(np.arange(0,10))) for pop in pops], 
                                # 'include': ['all'], 
                                'ylim':[-90,60],
                                'figSize':[24,15],
                                'saveFig': cfg.saveFigPath+'/'+cfg.filename+'traces'+'.png'}

# (includePre=['VPM__pop', 'TRN__pop', 'L6A__pop'], includePost=['VPM__pop', 'TRN__pop', 'L6A__pop'], figSize=(10,10), saveFig='figs_analysis/fullConn_topological_2023_09_21___connMatrix.png')
# feature='strength', orderBy='gid', groupBy='pop', groupByIntervalPre=None, groupByIntervalPost=None, graphType='matrix', removeWeightNorm=False, synOrConn='syn', synMech=None, connsFile=None, tagsFile=None, clim=None, figSize=(8, 8), fontSize=12, saveData=None, 

