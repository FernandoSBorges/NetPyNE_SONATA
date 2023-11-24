'''
netParams.py 

Contributors: joao.moreira@downstate.edu, salvadordura@gmail.com
'''

#------------------------------------------------------------------------------
# --- LIBRARY IMPORTS
#------------------------------------------------------------------------------
import NetPyNE_BBP
import numpy as np
from netpyne import specs
import pickle, json
import sys
import math

import pandas as pd
import os

netParams = specs.NetParams()   # object of class NetParams to store the network parameters
cfg       = specs.SimConfig()     # object of class SimConfig to store simulation configuration

# cfg.duration        = 100 # (ms)
cfg.duration        = 0.001*1e3 # (ms)
cfg.printRunTime    = 0.01 # (s)
cfg.verbose         = False

cfg.convertCellMorphologies=False
cfg.BBP_rootFolder                  = '/Users/joao/Research/Models/BBP/BBP_thalamus_microcircuit_2'
cfg.sonataConfigFile                = cfg.BBP_rootFolder+'/sonata/circuit_sonata.json'
cfg.morphologyFolder_h5             = cfg.BBP_rootFolder+'/sonata/morphologies_h5'
cfg.loadCellModel   = False    
cfg.convertCellModel= True
cfg.saveCellModel   = True
cfg.plotSynLocation = True
cfg.NetPyNE_rootFolder              = '/Users/joao/Research/Models/BBP/thalamus_netpyne'
cfg.NetPyNE_JSON_cells              = cfg.NetPyNE_rootFolder+'/cells/netpyne_morphologies'
cfg.NetPyNE_templateCells           = cfg.NetPyNE_rootFolder+'/mod'
cfg.NetPyNE_exportedCells           = cfg.NetPyNE_rootFolder+'/cells/morphologies_swc'
cfg.NetPyNE_L6A_JSON_cells          = cfg.NetPyNE_rootFolder+'/cells/S1_BBP_cells/'
cfg.modType                         = 'Prob_S1'

cfg.BBP_conn_properties             = cfg.NetPyNE_rootFolder+'/conn/calculate_BBP_conn_properties/BBP_conn_propeties.json'

# cfg.select_thal_gids=[37654, 38204, 35311, 34560, 38647, 35838, 38381, 40674, 38491, 38917, 36707, 39355, 34805, 35342, 40051, 36043, 38911, 37773, 36438, 39040]

# --- Loading the pre-generated network connectivity
try:    f = open('network_template.json')
except: sys.exit('Loading network failed')
store_network = json.load(f)

vpm_cellNumber = len(store_network['pops']['VPM__pop']['cellGids'])
trn_cellNumber = len(store_network['pops']['TRN__pop']['cellGids'])

randomizeNetwork=False
if randomizeNetwork:
    VPM_gids = NetPyNE_BBP.ConfigCfg.sampleGids(pop='VPL_TC', cfg_file=cfg.sonataConfigFile, gids = vpm_cellNumber)
    TRN_gids = NetPyNE_BBP.ConfigCfg.sampleGids(pop='Rt_RC',  cfg_file=cfg.sonataConfigFile, gids = trn_cellNumber)
else:
    VPM_gids=[40352, 34066, 36777, 40166, 37812, 38737, 40296, 39900, 35050, 36666, 42272, 38785, 40826, 42023, 33933, 34466, 38284, 34844, 40465, 39867, 41216, 42326, 39167, 33887, 33559, 36658, 36388, 38608, 36299, 39299, 39251, 41563, 39633, 37330, 34195, 40932, 41549, 37891, 36698, 37851, 37428, 36855, 37478, 40842, 35711, 36579, 40249, 37405, 37301, 42351, 34692, 36310, 35985, 41386, 36520, 41474, 42158, 35215, 42321, 40463, 38254, 41946, 40128, 40039, 40816, 35193, 41190, 40600, 40850, 40988, 39574, 34574, 35455, 36670, 35877, 39030, 40488, 33700, 40106, 37496, 37444, 37426, 34747, 37101, 40674, 35132, 38297, 37758, 41972, 36705, 37135, 38593, 36018, 39293, 34664, 37342, 41990, 41412, 39180, 34903, 40977, 34286, 36937, 37325, 40079, 40269, 34673, 33834, 41329, 35825, 37624, 41490, 37614, 42178, 33677, 36258, 36395, 35841, 34269, 37281, 34671, 39485, 35312, 34693, 37944, 36934, 40523, 38097, 39204, 42120, 41049, 33592, 34605, 38808, 40746, 37217, 34309, 36336, 33772, 37477, 37649, 35219, 39933, 35387, 42283, 41852, 40970, 33988, 33758, 40845, 33689, 40199, 36333, 41191, 34732, 42423, 33737, 40843, 38058, 38726, 34346, 38513, 39488, 40697, 33600, 37144, 38973, 38386, 38718, 40777, 39544, 40444, 34787, 39026, 34423, 37580, 40068, 33920, 36134, 38708, 40036, 38719, 39155, 33905, 37264, 37151, 40453, 41337, 37789, 33562, 39932, 39922, 35064, 38248, 35647, 37979, 40020, 34817, 35795, 36598, 33985, 35998, 42276, 39691, 41416, 36638, 36443, 41706, 41254, 41863, 36222, 36365, 34493, 40836, 33862, 35899, 40611, 37020, 34949, 34376, 34504, 40345, 37422, 35177, 39764, 37443, 36304, 34710, 39791, 34897, 36422, 38384, 38540, 40865, 41027, 36096, 36713, 39218, 39964, 34763, 38743, 39782, 41319, 40184, 42079, 34306, 41217, 41135, 34665, 39008, 34956, 39088, 37048, 41379, 39607, 35124, 41322, 40979, 37712, 36709, 34047, 42106, 36603, 37607, 38261, 35941, 35797, 35232, 34896, 39099, 39425, 39383, 36036, 35978, 39064, 37567, 39386, 40720, 38734, 36073, 42100, 33765, 41042, 37268, 38324, 42065, 37046, 40622, 35461, 38062, 38757, 36794, 38525, 34648, 37965, 36058, 39154, 38321, 39955, 36947, 34892, 36176, 41843, 41270, 35459, 40379, 36575, 41413, 40712, 41249, 34417, 34447, 39220, 37459, 34960, 38394, 34801, 39512]
    TRN_gids=[31364, 30106, 29095, 30404, 32411, 32689, 31136, 32952, 29223, 29226, 31158, 30716, 30874, 29694, 32244, 28976, 31159, 29987, 29060, 28963, 31917, 33413, 29656, 28828, 30883, 33262, 28657, 31312, 30553, 32956, 33496, 31127, 30686, 32101, 28629, 29205, 29588, 32390, 31050, 31497, 32746, 31493, 32515, 30601, 30161, 32655, 29546, 28867, 31388, 30041, 30876, 33378, 33405, 29303, 32878, 29152, 28751, 30860, 30265, 32556, 33156, 32386, 30801, 29146, 29527, 32576, 29159, 32174, 29217, 29736, 31557, 32851, 29401, 31530, 31839, 30846, 30879, 31170, 32483, 31756, 29935, 29172, 30493, 29094, 33437, 30013, 30440, 28812, 30307, 29266, 31156, 33485, 30093, 30987, 32245, 31512, 29027, 31093, 32087, 29125, 32921]
    cfg.loadCellModel           = True
    cfg.convertCellMorphologies = False
cfg.select_thal_gids=VPM_gids+TRN_gids

# sys.exit()

#------------------------------------------------------------------------------
# --- VERSION 
#------------------------------------------------------------------------------
netParams.version = 'thalamus_v00'

#------------------------------------------------------------------------------
#
# --- NETWORK PARAMETERS
#
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# --- General connectivity parameters
#------------------------------------------------------------------------------
netParams.scaleConnWeight = 1.0 # Connection weight scale factor (default if no model specified)
netParams.scaleConnWeightNetStims = 1.0 #0.5  # scale conn weight factor for NetStims
netParams.defaultThreshold = 0.0 # spike threshold, 10 mV is NetCon default, lower it for all cells

### reevaluate these values
netParams.defaultDelay = 2.0 # default conn delay (ms) # DEFAULT
netParams.propVelocity = 500.0 # propagation velocity (um/ms)
netParams.probLambda = 100.0  # length constant (lambda) for connection probability decay (um)

### maybe add the edge effect parameter to compensate for the error form distance dependant conn
# netParams.correctBorder = {‘threshold’: [200, 40, 200]}

netParams.defineCellShapes = True # JV 2021-02-23 - Added to fix the lack of the pt3d term in the cells, which make it unable to record i_membrane_

#------------------------------------------------------------------------------
# --- Load BBP circuit
#------------------------------------------------------------------------------
if cfg.convertCellMorphologies: NetPyNE_BBP.Conversion.convertCellMorphology(inputFolder_h5=cfg.morphologyFolder_h5,outputFolder_swc=cfg.NetPyNE_exportedCells,inputFolder_asc=cfg.morphologyFolder_asc)

# --- Load dictionary with thalamic circuit properties
# circuit_dict = NetPyNE_BBP.LoadBBPCircuit.getDataFrames ( cfg_file=cfg.sonataConfigFile, microcircuit_number=cfg.mc_number)
circuit_dict = NetPyNE_BBP.LoadBBPCircuit.getDataFrames_new( cfg_file=cfg.sonataConfigFile)

thal_gids = cfg.select_thal_gids

loading_failed=[]
if cfg.loadCellModel:
    NetPyNE_BBP.Prompt.headerMsg('Loading stored cells')
    for thal_gid in thal_gids:
        cfg.convertCellModel=False
        cell_properties = NetPyNE_BBP.LoadBBPCircuit.getNodeProperties(circuit_dict,thal_gid)
        # NetPyNE_BBP.Prompt.headerMsg('Loading cell '+str(thal_gid))
        try:        netParams.loadCellParams(label=cell_properties['mtype']+'__'+str(thal_gid),fileName=cfg.NetPyNE_JSON_cells+'/netpyne_'+str(thal_gid)+'.json')
        except:     loading_failed.append(thal_gid);print('Loading failed - Cell '+str(thal_gid)+' will be created in the next step')
if len(loading_failed)>0:thal_gids=loading_failed;cfg.convertCellModel=True

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# --- Cell Params
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
if cfg.convertCellModel:
    NetPyNE_BBP.Prompt.headerMsg('Processing cells')
    for thal_gid_ind,thal_gid in enumerate(thal_gids):
        print('\t>>\t',str(thal_gid),'\t|\t',str(len(thal_gids)-thal_gid_ind),' cells left')
        cell_properties = NetPyNE_BBP.LoadBBPCircuit.getNodeProperties(circuit_dict,thal_gid)
        thal_gid_conds={}
        # thal_gid_conds = df_thalamus_neurons.loc[thal_gid].to_dict()      # storing cell properties in NetPyNE model
        thal_gid_conds.update({'cell_index':thal_gid})
        netParams.importCellParams(
            label           = cell_properties['mtype']+'__'+str(thal_gid), 
            fileName        = cfg.NetPyNE_templateCells+'/'+cell_properties['etype']+'.hoc', 
            cellName        = cell_properties['etype'], 
            conds           = {'cellType':str(thal_gid)},
            # conds         = {'cellType':str(thal_gid),'mtype':cell_properties['mtype']},
            # conds         = {'cellType':cell_properties['mtype']},
            cellArgs        = [thal_gid,cfg.NetPyNE_exportedCells,cell_properties['morphology']+'.swc'], 
            importSynMechs  = True, 
            somaAtOrigin    = False, 
            cellInstance    = False,
        )
        # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------

        # --- Fixing bug in sections that don't have a pt3d list (or have an empty dictionary instead)
        for sec in netParams.cellParams[cell_properties['mtype']+'__'+str(thal_gid)]['secs'].keys():
            if type(netParams.cellParams[cell_properties['mtype']+'__'+str(thal_gid)]['secs'][sec]['geom']['pt3d']) is not list:
                netParams.cellParams[cell_properties['mtype']+'__'+str(thal_gid)]['secs'][sec]['geom']['pt3d']=[]
        
        # --- Creating secLists based on the path distance to soma during network setup
        for sec in netParams.cellParams[cell_properties['mtype']+'__'+str(thal_gid)]['secs'].keys():
            cell_dict = netParams.cellParams[cell_properties['mtype']+'__'+str(thal_gid)]
            soma_pathDist = NetPyNE_BBP.Utils.pathDistance(cell_dict=cell_dict,root_sec='soma_0',ignoreAxon=True)
            basalDends    = NetPyNE_BBP.Utils.findBasalDends(cell_dict,root_sec='soma_0',ignoreAxon=True)
            secLists_dict = NetPyNE_BBP.Utils.secListFromPathDistance(soma_pathDist,basalDendrites=basalDends)
            for secList in secLists_dict.keys(): netParams.cellParams[cell_properties['mtype']+'__'+str(thal_gid)]['secLists'].update({secList:secLists_dict[secList]})
        
        # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        if cfg.saveCellModel: netParams.saveCellParams(label=cell_properties['mtype']+'__'+str(thal_gid),fileName=cfg.NetPyNE_JSON_cells+'/netpyne_'+str(thal_gid)+'.json')
        
        # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------

for gids_list in [VPM_gids,TRN_gids]:
    # --- Adding cell diversity rule
    for thal_gid in gids_list:
        cell_properties = circuit_dict['thalamus_neurons'].loc[thal_gid]
        # --- Adds each GID morphology into the cell diversity dictionary
        netParams.cellParams[cell_properties['mtype']+'__'+str(thal_gid)].update({'diversityFraction':1/len(gids_list)})
        # --- Changes cellType to the same value, so that the different morphologies are added to the same pop can be identified by connParams
        netParams.cellParams[cell_properties['mtype']+'__'+str(thal_gid)]['conds']['cellType']=cell_properties['mtype']

# --- Creates a single pop for all cells and defines popParams with diversity rule so that each cell has a different morphology, based on its GIDs. This aproach makes it easier to target stims and run parallel sims
clone_pops = ['VPM__pop','TRN__pop']
pop_vars   = ['xRange','yRange','zRange','xnormRange','ynormRange','znormRange','numCells']

netParams.popParams['VPM__pop']={'cellType':'VPL_TC','diversity': True}
netParams.popParams['TRN__pop']={'cellType':'Rt_RC', 'diversity': True}

for pop_name in clone_pops:
    for pop_var in pop_vars: netParams.popParams[pop_name].update({pop_var:store_network['pops'][pop_name]['tags'][pop_var]})

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# --- L6A cells and pop
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
import Build_Net as BN
simplifyL6A=False
if simplifyL6A: 
    netParams.cellParams.update(BN.BuildNetwork.getCellTemplate(    template='izhi',pops=['L6A'],cell_flag='__cell'))
    netParams.popParams.update( BN.BuildNetwork.getPopTemplate(     pops=['L6A'],center_point=500,pop_flag='__pop',cell_flag='__cell'))
else:
    netParams.cellParams.update(BN.BuildNetwork.getL6ACellTemplate( cellsFolder=cfg.NetPyNE_L6A_JSON_cells))
    netParams.popParams.update( BN.BuildNetwork.getPopTemplate(     pops=['L6A'],center_point=500,pop_flag='__pop',cell_flag='__cell',diversity=True,volume_shape='cube'))
    print('\t>>\tWarning: Add code to remove AXON segment from L6A cells')

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# --- MLe cells and pops
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
align_cells=False
useVecStim = True
if useVecStim:

    target_angle = 90
    # samplingInterv = 25     # us
    samplingInterv = 1000   # us

    NetPyNE_rootFolder              = '/Users/joao/Research/Models/BBP/thalamus_netpyne'
    stim_source_code_folder         = '/stims/Minnery_2003/fig1'
    save_deflection_model_folder    = NetPyNE_rootFolder+stim_source_code_folder+'/deflection_model/'   # Folder to save the output deflection model and sampled raster
    from GenerateStimDataset import SampleData
    deflection_dict = SampleData.LoadJSON(file_path=save_deflection_model_folder+'deflectionModel_'+str(target_angle)+'|deg.json')
    file_path = save_deflection_model_folder+'spike_dicts/mleSpikes_deflectionAngle|'+str(target_angle)+'_samplingBins|'+str(samplingInterv)+'us.json'
    try:
        print('Loading spikes from \n',
              file_path)
        spikes_dict = SampleData.LoadJSON(file_path=file_path)
    except:
        file_path = save_deflection_model_folder+'spike_dicts/mleSpikes_deflectionAngle|0_samplingBins|1000us.json'
        print('Raster Loading failed - Loading sample raster: deflection 0 degrees | sampling interval 1000 us')
        spikes_dict = SampleData.LoadJSON(file_path=file_path)


    netParams.popParams.update( BN.BuildNetwork.getMLePopTemplate_VecStim(spkts_dict=spikes_dict['spkts']))
else:
    netParams.cellParams.update(BN.BuildNetwork.getCellTemplate(        template='izhi',pops=['MLe'],cell_flag='__cell'))
    netParams.popParams.update( BN.BuildNetwork.getMLePopTemplate(      align_cells=align_cells))
# import sys;sys.exit()
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# --- Synaptic Mechanisms
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# --- Dummy exc and inh synMechs for testing
netParams.synMechParams = BN.BuildNetwork.getSynMechParams()

# --- Selecting type of MOD file to be used
if   cfg.modType == 'Det':       modAMPANMDA = 'DetAMPANMDA';            modGABA     = 'DetGABAAB'              # S1 Deterministic  implementation of the BBP mod files
elif cfg.modType == 'Prob_S1':   modAMPANMDA = 'ProbAMPANMDA_EMS_S1';    modGABA     = 'ProbGABAAB_EMS_S1'      # S1 Probabilistic  implementation of the BBP mod files
else:                            modAMPANMDA = 'ProbAMPANMDA_EMS';       modGABA     = 'ProbGABAAB_EMS'         # Original Thalamus implementation of the BBP mod files
print('\n\t>>\tMOD template\tAMPA: ', modAMPANMDA, '\tGABA: ', modGABA)

TableS1_barreloidThalamus   = NetPyNE_BBP.StoreParameters.getTableS1_barreloidThalamus()
for pre_pop in TableS1_barreloidThalamus.keys():
    if   pre_pop == 'L6A':  mod_file = modAMPANMDA; mech_flag='exc'
    elif pre_pop == 'MLe':  mod_file = modAMPANMDA; mech_flag='exc'
    elif pre_pop == 'VPM':  mod_file = modAMPANMDA; mech_flag='exc'
    else:                   mod_file = modGABA;     mech_flag='inh'
    for post_pop in TableS1_barreloidThalamus[pre_pop].keys():
        syn_values = TableS1_barreloidThalamus[pre_pop][post_pop]
        if cfg.modType == 'Prob':syn_values['paper_reference_values'].update({'n_rrp_vesicles':1}) # only present in the Prob MOD, not in the Det MOD
        edge_dict = NetPyNE_BBP.CreateNetPyNE.modTemplate(syn_values['paper_reference_values'],mod_file)
        netParams.synMechParams.update({'syn|'+pre_pop+'|'+post_pop+'|'+mech_flag:edge_dict})

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# --- Connectivity
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# --- Loads the conn data calculated from the projections in the BBP somatosensory thalamus model
conn_data = BN.NetworkConversion.convertConnProperties(filePath=cfg.BBP_conn_properties)

# --- Connect the cells based on the loaded network
#########       TESTING CONNECTIONS USING THE buildNet.py IMPLEMENTATION TO SEE IF THEY MATCH HERE      ########

# --- Connect using cilinder with exponential decay or uniform distribution

if      simplifyL6A:    secTarget = 'soma_0'    # single section
else:                   secTarget = 'basal'     # secList

# ########## NOT USING ANYMORE - START ##########

# # prePop,postPop,sec/listOfSections/secList
# connect_pops={
#                 'cilinder_exp':[
#                                 # ('VPM','TRN','inputs__proximal'),       # verify - assumed proximal because input is driver - Excitatory Depressing (Fox, 2008)
#                                 # ('VPM','L6A',secTarget),                # soma_0 for point process and basal dendrites for detailed morphology (driver input)
#                                 # # ('L6A','VPM','inputs__distal'),         # (Jones, 2002)
#                                 # # ('L6A','TRN','inputs__distal'),         # verify - assumed distal because input is modulator
#                                 ],
#                 'cilinder_uni':[
#                                 # ('L6A','VPM','inputs__distal'),         # (Jones, 2002) / # See figs from (Deschênes, 1998)
#                                 # ('L6A','TRN','inputs__distal'),         # verify - assumed distal because input is modulator
#                                 # ('TRN','VPM','inputs__intermediate'),   # (Jones, 2002)
#                                 # ('TRN','TRN','inputs__distal'),         # verify - assumed distal (dendro-dendritic synapses)
#                                 ],
#                 }

# for conn_type in connect_pops.keys():
#     for (pre_pop,post_pop,secList) in connect_pops[conn_type]:
#         if   conn_type=='cilinder_exp':(conn_method,conn_rule)=BN.BuildNetwork.cilinderProjection_expDecay(pre_pop,post_pop)
#         elif conn_type=='cilinder_uni':(conn_method,conn_rule)=BN.BuildNetwork.cilinderProjection_uniform( pre_pop,post_pop)

#         # Adding a flag with the words that NetPyNE recognizes to color the connection blue in plot2Dnet
#         if 'TRN' in pre_pop:mech_flag='inh'
#         else:               mech_flag='exc'

#         if 'syn|'+pre_pop+'|'+post_pop+'|'+mech_flag in netParams.synMechParams.keys(): syn_mech='syn|'+pre_pop+'|'+post_pop+'|'+mech_flag
#         else:
#             if pre_pop=='TRN':syn_mech='inh'
#             else:             syn_mech='exc'
#             print('\t>>\tUsing a temporary <',syn_mech,'> mechanism to represent pathway ',pre_pop,'->',post_pop ,' - verify which one is the proper syn mech to add')

#         # try:syn_mech='syn|'+pre_pop+'|'+post_pop+'|'+mech_flag
#         # except:        
#         #     if pre_pop=='TRN':syn_mech='inh'
#         #     else:             syn_mech='exc'

#         # syns per conn - obs: conn_data is organized as conn_data[post_pop]['chem'][pre_pop][synsPerConn/convergence][MEAN,STD]
#         try:    
#             syns_per_conn = round(conn_data[post_pop]['chem'][pre_pop]['synsPerConn'][0])
#             if (conn_data[post_pop]['chem'][pre_pop]['synsPerConn'][0]>0) and int(conn_data[post_pop]['chem'][pre_pop]['synsPerConn'][0])==0: syns_per_conn=1
#         except: syns_per_conn = 1

#         conn_dict = BN.BuildNetwork.getConnDict(pre_pop,post_pop,conn_method,conn_rule,syn_mech,syns_per_conn,target_secs=secList)
#         netParams.connParams.update(conn_dict)

# ########## NOT USING ANYMORE - END ##########

conn_MLe_VPM=False
conn_L6A_VPM=True
conn_L6A_TRN=False
conn_TRN_VPM=False
conn_TRN_TRN=False
conn_TRN_TRN_electrical=False
conn_VPM_TRN=False
conn_VPM_L6A=False

# --- MLe to VPM connections
if conn_MLe_VPM:

    # # syns per conn - obs: conn_data is organized as conn_data[post_pop]['chem'][pre_pop][synsPerConn/convergence][MEAN,STD]
    # try:    
    #     syns_per_conn = round(conn_data['VPM']['MLe']['synsPerConn'][0])
    #     if (conn_data['VPM']['MLe']['synsPerConn'][0]>0) and int(conn_data['VPM']['MLe']['synsPerConn'][0])==0: syns_per_conn=1
    #     # print('MLe synsPerConn = ', syns_per_conn)
    # except: syns_per_conn = 1

    '''
    (Takeuchi, 2017)
    # of boutons / mle fiber:   40.50919230368733
    '''
    syns_per_conn = 40

    synMech = 'syn|MLe|VPM|exc'
    conn_type = 'chem'

    theta_pops = [pop_name for pop_name in netParams.popParams.keys() if 'MLe' in pop_name]
    
    if useVecStim:
        print('\t>>\tAdding VecStim MLe connections - One rule per (cell/pop)')
        for pre_pop in theta_pops:
            mle_pop = pre_pop.split('__')[0]
            mle_degree = int(mle_pop.split('@')[1])
            print(pre_pop,'\t',mle_pop,'\t',mle_degree,'\t')

            secList = 'inputs__proximal'        # (Jones, 2002) - driver input targets proximal sections
            
            # # --- Changing a few projections from 'exc' to 'inh' to see if the topology is correct using plot2dnet
            # if mle_degree<20: synMech='inh'
            # else: synMech='syn|MLe|VPM|exc'
            # print('\t>>\t',synMech)

            (mle_conn_method, mle_conn_rule)=BN.BuildNetwork.laminarProjection_VecStim(pre_angle=mle_degree)
            conn_dict = BN.BuildNetwork.getConnDict(pre_pop=mle_pop,post_pop='VPM',conn_method=mle_conn_method,conn_rule=mle_conn_rule,syn_mech=synMech,syns_per_conn=syns_per_conn,conn_type=conn_type)
            netParams.connParams.update(conn_dict)

    else:
        (mle_conn_method, mle_conn_rule)=BN.BuildNetwork.laminarProjection()

        for pre_pop in theta_pops:
            mle_pop = pre_pop.split('__')[0]
            secList = 'inputs__proximal'        # (Jones, 2002) - driver input targets proximal sections


            conn_dict = BN.BuildNetwork.getConnDict(pre_pop=mle_pop,post_pop='VPM',conn_method=mle_conn_method,conn_rule=mle_conn_rule,syn_mech=synMech,syns_per_conn=syns_per_conn,conn_type=conn_type)
            netParams.connParams.update(conn_dict)

# scale_prob = 40/12.483479074117373

if conn_L6A_VPM:
    pre_pop  = 'L6A'
    post_pop = 'VPM'
    conn_type = 'chem'
    try:    
        syns_per_conn = round(conn_data[post_pop]['chem'][pre_pop]['synsPerConn'][0])
        if (conn_data[post_pop]['chem'][pre_pop]['synsPerConn'][0]>0) and int(conn_data[post_pop]['chem'][pre_pop]['synsPerConn'][0])==0: syns_per_conn=1
    except: syns_per_conn = 1

    conn_method_1D, conn_rule_1D = BN.BuildNetwork.distanceBasedProbability_1D_exponential(pre_pop=pre_pop,post_pop=post_pop,
                                                                                           pre_pop_axis='x',post_pop_axis='y',
                                                                                           conn_prob=2.0827527202658684,center_point=500)
    conn_dict_1D = BN.BuildNetwork.getConnDict( pre_pop         = pre_pop,
                                                post_pop        = post_pop,
                                                conn_method     = conn_method_1D,
                                                conn_rule       = conn_rule_1D,
                                                syn_mech        = 'syn|L6A|VPM|exc',
                                                syns_per_conn   = syns_per_conn,
                                                conn_type       = conn_type,
                                                target_secs     = 'inputs__distal',
                                                )
    netParams.connParams.update(conn_dict_1D)

if conn_L6A_TRN:
    pre_pop  = 'L6A'
    post_pop = 'TRN'
    conn_type = 'chem'
    try:    
        syns_per_conn = round(conn_data[post_pop]['chem'][pre_pop]['synsPerConn'][0])
        if (conn_data[post_pop]['chem'][pre_pop]['synsPerConn'][0]>0) and int(conn_data[post_pop]['chem'][pre_pop]['synsPerConn'][0])==0: syns_per_conn=1
    except: syns_per_conn = 1

    conn_method_1D, conn_rule_1D = BN.BuildNetwork.distanceBasedProbability_1D_exponential(pre_pop=pre_pop,post_pop=post_pop,
                                                                                           pre_pop_axis='x',post_pop_axis='x',
                                                                                           conn_prob=0.8010587385637955,center_point=500)
    conn_dict_1D = BN.BuildNetwork.getConnDict( pre_pop         = pre_pop,
                                                post_pop        = post_pop,
                                                conn_method     = conn_method_1D,
                                                conn_rule       = conn_rule_1D,
                                                syn_mech        = 'syn|L6A|TRN|exc',
                                                syns_per_conn   = syns_per_conn,
                                                conn_type       = conn_type,
                                                target_secs     = 'inputs__distal')
    netParams.connParams.update(conn_dict_1D)

if conn_TRN_VPM:
    pre_pop  = 'TRN'
    post_pop = 'VPM'
    conn_type = 'chem'
    try:    
        syns_per_conn = round(conn_data[post_pop]['chem'][pre_pop]['synsPerConn'][0])
        if (conn_data[post_pop]['chem'][pre_pop]['synsPerConn'][0]>0) and int(conn_data[post_pop]['chem'][pre_pop]['synsPerConn'][0])==0: syns_per_conn=1
    except: syns_per_conn = 1

    # --- Number of conns from MLe to VPM, used as a parameter to estimate the TRN to VPM convergence
    # conns_MLe_to_VPM = 16.67924528301887 # must be updated if values change
    conns_MLe_to_VPM = 55.59748427672956 # must be updated if values change
    
    '''
    # Ratio of (TRN-to->VPM)/(MLe-to->VPM)
    '''
    # # BBP ratio           = 454.1682034044653/122.97733554449974  = 3.6931049237127365
    ratio_TRNtoVPM_MLetoVPM = (conn_data['VPM']['TRN']['synsPerConn'][0]*conn_data['VPM']['TRN']['convergence'][0])/(conn_data['VPM']['MLe']['synsPerConn'][0]*conn_data['VPM']['MLe']['convergence'][0])
    # (çavdar, 2007) ratio  = 48/28                                 = 1.7142857142857142    # --- Testing this value instead to check if BBP (TRN->VPM) conns are not too high
    # ratio_TRNtoVPM_MLetoVPM = 48/28

    convergence_ = (ratio_TRNtoVPM_MLetoVPM)*(conns_MLe_to_VPM/syns_per_conn)
    
    '''
    conns = convergence_ * syns_per_conn_

    convergence_TRN-to->VPM = ((TRN-to->VPM)/(MLe-to->VPM))*(conns_MLe_to_VPM/syns_per_conn_TRN_to_VPM)
    '''

    conn_method_ = 'convergence'
    conn_rule_   = convergence_

    conn_dict_1D = BN.BuildNetwork.getConnDict( pre_pop         = pre_pop,
                                                post_pop        = post_pop,
                                                conn_method     = conn_method_,
                                                conn_rule       = conn_rule_,
                                                syn_mech        = 'syn|TRN|VPM|inh',
                                                syns_per_conn   = syns_per_conn,
                                                conn_type       = conn_type,
                                                target_secs     = 'inputs__intermediate')
    netParams.connParams.update(conn_dict_1D)

if conn_TRN_TRN:
    pre_pop  = 'TRN'
    post_pop = 'TRN'
    conn_type = 'chem'
    try:    
        syns_per_conn = round(conn_data[post_pop]['chem'][pre_pop]['synsPerConn'][0])
        if (conn_data[post_pop]['chem'][pre_pop]['synsPerConn'][0]>0) and int(conn_data[post_pop]['chem'][pre_pop]['synsPerConn'][0])==0: syns_per_conn=1
    except: syns_per_conn = 1

    convergence_ = conn_data[post_pop]['chem'][pre_pop]['convergence'][0]
    
    conn_method_ = 'convergence'
    conn_rule_   = convergence_
    print('TRN TRN chemical: ', convergence_)

    conn_dict_1D = BN.BuildNetwork.getConnDict( pre_pop         = pre_pop,
                                                post_pop        = post_pop,
                                                conn_method     = conn_method_,
                                                conn_rule       = conn_rule_,
                                                syn_mech        = 'syn|TRN|TRN|inh',
                                                syns_per_conn   = syns_per_conn,
                                                conn_type       = conn_type,
                                                target_secs     = 'inputs__distal')
    netParams.connParams.update(conn_dict_1D)

if conn_TRN_TRN_electrical:
    pre_pop  = 'TRN'
    post_pop = 'TRN'
    conn_type = 'elec'
    try:    
        syns_per_conn = round(conn_data[post_pop]['elec'][pre_pop]['synsPerConn'][0])
        if (conn_data[post_pop]['elec'][pre_pop]['synsPerConn'][0]>0) and int(conn_data[post_pop]['elec'][pre_pop]['synsPerConn'][0])==0: syns_per_conn=1
    except: syns_per_conn = 1

    convergence_ = conn_data[post_pop]['elec'][pre_pop]['convergence'][0]
    
    conn_method_ = 'convergence'
    conn_rule_   = convergence_
    print('TRN TRN electrical: ', convergence_)

    conn_dict_1D = BN.BuildNetwork.getConnDict( pre_pop         = pre_pop,
                                                post_pop        = post_pop,
                                                conn_method     = conn_method_,
                                                conn_rule       = conn_rule_,
                                                # syn_mech        = 'Gap',
                                                syn_mech        = 'esyn',
                                                # syn_mech        = 'syn|TRN|TRN|inh',
                                                syns_per_conn   = syns_per_conn,
                                                conn_type       = conn_type,
                                                target_secs     = 'inputs__distal')
    netParams.connParams.update(conn_dict_1D)

if conn_VPM_TRN:
    pre_pop  = 'VPM'
    post_pop = 'TRN'
    conn_type = 'chem'
    try:    
        syns_per_conn = round(conn_data[post_pop]['chem'][pre_pop]['synsPerConn'][0])
        if (conn_data[post_pop]['chem'][pre_pop]['synsPerConn'][0]>0) and int(conn_data[post_pop]['chem'][pre_pop]['synsPerConn'][0])==0: syns_per_conn=1
    except: syns_per_conn = 1

    conn_method_1D, conn_rule_1D = BN.BuildNetwork.distanceBasedProbability_1D_exponential(pre_pop=pre_pop,post_pop=post_pop,
                                                                                           pre_pop_axis='y',post_pop_axis='x',
                                                                                           conn_prob=0.48063524313827727,center_point=500)
    conn_dict_1D = BN.BuildNetwork.getConnDict( pre_pop         = pre_pop,
                                                post_pop        = post_pop,
                                                conn_method     = conn_method_1D,
                                                conn_rule       = conn_rule_1D,
                                                syn_mech        = 'syn|VPM|TRN|exc',
                                                syns_per_conn   = syns_per_conn,
                                                conn_type       = conn_type,
                                                target_secs     = 'inputs__proximal')
    netParams.connParams.update(conn_dict_1D)

if conn_VPM_L6A:
    print()
    CT_cells            = 666
    CT_bouton_per_cell  = 133      # from (Meyer, 2010) 
    TC_cells            = 318
    TC_syns_per_conn    = 9        # from S1-netpyne model
    
    CT_boutons  = CT_cells*CT_bouton_per_cell
    TC_syns     = TC_cells*TC_syns_per_conn

    TC_CT_conns = CT_boutons/TC_syns
    print('TC cell target convergence onto CT neurons: ',round(TC_CT_conns))

    pre_pop  = 'VPM'
    post_pop = 'L6A'
    conn_type = 'chem'

    conn_method_, conn_rule_ = BN.BuildNetwork.positionBasedProbability_1D_linear(  pre_pop=pre_pop,post_pop=post_pop,pre_pop_axis='y',post_pop_axis='y',
                                                                                    conn_prob=0.075,
                                                                                    coef=[-185.45820847,  228.16359477], # from (Meyer,2010) m and b terms from y=(mx+b) x=depth/y=#_of_boutons
                                                                                    center_point=500,
                                                                                    )
    
    conn_dict_ = BN.BuildNetwork.getConnDict(   pre_pop         = pre_pop,
                                                post_pop        = post_pop,
                                                conn_method     = conn_method_,
                                                conn_rule       = conn_rule_,
                                                syn_mech        = 'syn|VPM|L6A|exc',
                                                # syn_mech        = 'exc', # temporary - until I implement the syn|VPM|L6A|exc mechanism
                                                syns_per_conn   = TC_syns_per_conn,
                                                conn_type       = conn_type,
                                                target_secs     = secTarget)
    netParams.connParams.update(conn_dict_)

    
    # 133*666/
    # (318*9)

    # see (Meyer 2010)



# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# --- Stimulation
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------











# ############################################################################################################
# # --- Redistribute synapses into postsynaptic cell soma/dendrites according to literature (Jones, 2002)

# for conn in netParams.connParams.keys():
#     redistribute_conns = {
#                             ('MLe','VPM'):'inputs__proximal',
#                             ('TRN','VPM'):'inputs__intermediate',
#                             ('L6A','VPM'):'inputs__distal',
#                             }
#     conn_=conn.split('|')
#     if '@' in conn_[1]: pre_pop  = conn_[1].split('@')[0]   # removes the cell index in MLe pops
#     else:               pre_pop  = conn_[1]
#     post_pop = conn_[2]
#     if (pre_pop,post_pop) in redistribute_conns.keys():




#         netParams.subConnParams[conn] = {   'preConds':  {'pop': conn_[1]}, 
#                                             'postConds': {'pop': conn_[2]},
#                                             'sec': ,   # probability of connection
#                                             'density': {'type':'1Dmap','gridY':gridY,'gridValues':gridValues,'fixedSomaY':0}
#             } 

############################################################################################################
# --- Plotting
# cfg.analysis['plotConn']  = {'includePre':['VPM__pop', 'TRN__pop', 'L6A__pop'],'includePost':['VPM__pop', 'TRN__pop', 'L6A__pop'],'figSize':(10, 10),'saveFig': 'figs_analysis/fullConn_topological_2023_09_21___connMatrix.png'}                                                 # plot connectivity matrix
# cfg.analysis['plot2Dnet'] = {'figSize':(8, 20), 'saveFig': 'figs_analysis/fullConn_topological_2023_10_05___mleConns_aligned_testInh_.png'} # plot 2D cell positions and connections
# cfg.analysis['plot2Dnet'] = {'figSize':(8, 20), 'saveFig': 'figs_analysis/fullConn_topological_2023_09_21___2Dnetwork.png'} # plot 2D cell positions and connections

# (includePre=['VPM__pop', 'TRN__pop', 'L6A__pop'], includePost=['VPM__pop', 'TRN__pop', 'L6A__pop'], figSize=(10,10), saveFig='figs_analysis/fullConn_topological_2023_09_21___connMatrix.png')
# feature='strength', orderBy='gid', groupBy='pop', groupByIntervalPre=None, groupByIntervalPost=None, graphType='matrix', removeWeightNorm=False, synOrConn='syn', synMech=None, connsFile=None, tagsFile=None, clim=None, figSize=(8, 8), fontSize=12, saveData=None, 

# cfg.analysis['plotRaster'] = {'figSize':(25, 20), 'saveFig': 'figs_analysis/raster_MLe_input.png'} # Plot a raster

############################################################################################################
# --- Running simulation
from netpyne import sim
sim.createSimulateAnalyze(netParams = netParams, simConfig = cfg)
# sim.create(netParams = netParams, simConfig = cfg)

############################################################################################################
runPostAnalysis=False
if runPostAnalysis:
    # --- Verify conns
    verifyConns=False
    if verifyConns:
        for cell_ind in range(len(sim.net.cells)):
            print('cell ', cell_ind)
            secList_proximal    = sim.net.cells[cell_ind].secLists['inputs__proximal']
            secList_intermediate= sim.net.cells[cell_ind].secLists['inputs__intermediate']
            secList_distal      = sim.net.cells[cell_ind].secLists['inputs__distal']
            if len(sim.net.cells[cell_ind].conns)==0:
                print('No conns')
                continue
            for conn_ind in range(len(sim.net.cells[cell_ind].conns)):
                conn_sec = sim.net.cells[cell_ind].conns[conn_ind]['sec']
                if   conn_sec in secList_proximal:      print('0 - proximal conn')
                elif conn_sec in secList_intermediate:  print('1 - intermediate conn')
                elif conn_sec in secList_distal:        print('2 - distal conn')

    ############################################################################################################
    # === Evaluate connectivity
    print('\t>>\tRunning network evaluation - Divergence and Convergence values')
    from EvaluateModel import Evaluate as Eval

    pop_gids  = Eval.getPopGids(sim)
    pop_names = list(pop_gids.keys())

    pop_convergence  = Eval.getPopConvergence(sim)
    pop_divergence   = Eval.getPopDivergence(sim)

    cell_convergence = Eval.getCellConvergence(sim)
    cell_divergence  = Eval.getCellDivergence(sim)

    # pop_gids,pop_names,pop_convergence,pop_divergence,cell_convergence,cell_divergence = Eval.runAll(sim)

    # --- Pathway ratios for cell convergence/divergence
    convergence_ratio = Eval.getConvergenceRatio(sim)
    divergence_ratio  = Eval.getDivergenceRatio(sim)

    print('\n------------------------------------------------------------------------------------------------------------------------------\n')

    print('\t>>\t--- Convergence Data ---')
    print('\t>>\tConvergence Absolute [(number of conns * syns per conn) for each cell]:  ')
    [print('\t\t',key,cell_convergence[key]) for key in cell_convergence.keys()]
    print('\t>>\tRatio - Convergence Absolute [(number of conns * syns per conn) for each cell]:  ')
    [print('\t\t',key,convergence_ratio[key]) for key in convergence_ratio.keys()]

    print('\n\t\t - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n')
    print('\t>>\tPathway ratios:  ')
    try: print('\t>>\t (L6A->VPM)/(L6A->MLe) ratio: ', cell_convergence['VPM']['L6A']/cell_convergence['VPM']['MLe'], '\t-- target - BBP: ', 581.4493258129902/122.97733554449974)
    except: 0
    try: print('\t>>\t (TRN->VPM)/(MLe->VPM) ratio: ', cell_convergence['VPM']['TRN']/cell_convergence['VPM']['MLe'], '\t-- target - BBP: ', 454.1682034044653/122.97733554449974)
    except: 0

    try: print('\t>>\t (L6A->VPM)/(L6A->TRN) ratio: ', cell_convergence['VPM']['L6A']/cell_convergence['TRN']['L6A'], '\t-- target - BBP: ', 581.4493258129902/255.23635294982)
    except: 0
    try: print('\t>>\t (L6A->TRN)/(VPM->TRN) ratio: ', cell_convergence['TRN']['L6A']/cell_convergence['TRN']['VPM'], '\t-- target - BBP: ', 255.23635294982/161.43545983636517)
    except: 0

    print('\n\t\t - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n')

    print('\t>>\t--- Divergence Data ---')
    print('\t>>\tDivergence Absolute [(number of conns * syns per conn) for each cell]:  ')
    [print('\t\t',key,cell_divergence[key]) for key in cell_divergence.keys()]
    print('\t>>\tRatio - Divergence Absolute [(number of conns * syns per conn) for each cell]:  ')
    [print('\t\t',key,divergence_ratio[key]) for key in divergence_ratio.keys()]

    print('\n------------------------------------------------------------------------------------------------------------------------------\n')

    # --- BBP model stats
    print('\t>>\t--- BBP Convergence Data ---')
    print('\t>>\tBBP Convergence Absolue:  ')
    [print('\t\t',key2,'-to->',key,': ',conn_data[key][key2]['synsPerConn'][0]*conn_data[key][key2]['convergence'][0]) for key in conn_data.keys() for key2 in conn_data[key].keys()]

    print('\n\t\t - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n')
    print('\t>>\tPathway ratios BBP:  ')
    print('\t>>\t (L6A->VPM)/(L6A->MLe) ratio: ', 581.4493258129902/122.97733554449974)
    print('\t>>\t (TRN->VPM)/(MLe->VPM) ratio: ', 454.1682034044653/122.97733554449974)
    print('\t>>\t (L6A->VPM)/(L6A->TRN) ratio: ', 581.4493258129902/255.23635294982)
    print('\t>>\t (L6A->TRN)/(VPM->TRN) ratio: ', 255.23635294982/161.43545983636517)

    print('\n------------------------------------------------------------------------------------------------------------------------------\n')

    # --- Ratio of synapses into VB thalamus from (Van Horn and Sherman, 2007)
    driver_RL   = 28
    modul_RS    = 146
    trn_F       = 48
    sum_all = driver_RL+modul_RS+trn_F
    print('\t>>\t--- Convergence Data (Van Horn and Sherman, 2007) ---')
    print('\t>>\tConvergence Proportional:  ')
    print('\t\t MLe-to->VPM: ', driver_RL)
    print('\t\t L6A-to->VPM: ', modul_RS)
    print('\t\t TRN-to->VPM: ', trn_F)

    print('\t>>\tConvergence Ratio:  ')
    print('\t\t MLe-to->VPM: ', driver_RL/sum_all)
    print('\t\t L6A-to->VPM: ', modul_RS/sum_all)
    print('\t\t TRN-to->VPM: ', trn_F/sum_all)

    print('\n\t\t - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n')
    print('\t>>\tPathway ratios:  ')
    print('\t>>\t (L6A->VPM)/(MLe->VPM) ratio: ', modul_RS/driver_RL)
    print('\t>>\t (TRN->VPM)/(MLe->VPM) ratio: ', trn_F/driver_RL)
    # print('\t>>\tBBP Convergence ratio:  ')
    # [print('\t\t',key,conn_data[key]) for key in conn_data.keys()]

    print('\n------------------------------------------------------------------------------------------------------------------------------\n')
    # print('\n\t\t - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n')

    '''
    (çavdar, 2011) - VB terminals and synapses - Fig 8
                    driver - RL             modulator - RS          trn - F
    terminals       4.91                    86.81818181818181       6.092307692307692
    synapses        12.248803827751196      79.6969696969697        7.569230769230771       (USE THIS VALUE)
    '''
    driver_RL   = 12.248803827751196
    modul_RS    = 79.6969696969697
    trn_F       = 7.569230769230771
    sum_all = driver_RL+modul_RS+trn_F

    print('\t>>\t--- Convergence Data (çavdar, 2007) ---')
    print('\t>>\tConvergence Proportional:  ')
    print('\t\t MLe-to->VPM: ', driver_RL)
    print('\t\t L6A-to->VPM: ', modul_RS)
    print('\t\t TRN-to->VPM: ', trn_F)

    print('\t>>\tConvergence Ratio:  ')
    print('\t\t MLe-to->VPM: ', driver_RL/sum_all)
    print('\t\t L6A-to->VPM: ', modul_RS/sum_all)
    print('\t\t TRN-to->VPM: ', trn_F/sum_all)

    print('\n\t\t - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n')
    print('\t>>\tPathway ratios:  ')
    print('\t>>\t (L6A->VPM)/(MLe->VPM) ratio: ', modul_RS/driver_RL)
    print('\t>>\t (TRN->VPM)/(MLe->VPM) ratio: ', trn_F/driver_RL)

    print('\n------------------------------------------------------------------------------------------------------------------------------\n')


    # --- Bar Plot with conn values
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 36})
    plt.figure( figsize=(20,20))

    pathway_ratios={
                    '(toVPM)_L6A/MLe':  {'BBP':581.4493258129902/122.97733554449974,    'model':cell_convergence['VPM']['L6A']/cell_convergence['VPM']['MLe']},
                    '(toVPM)_TRN/MLe':  {'BBP':454.1682034044653/122.97733554449974,    'model':cell_convergence['VPM']['TRN']/cell_convergence['VPM']['MLe']},
                    '(fromL6A)_VPM/TRN':{'BBP':581.4493258129902/255.23635294982,       'model':cell_convergence['VPM']['L6A']/cell_convergence['TRN']['L6A']},
                    '(toTRN)_L6A/VPM':  {'BBP':255.23635294982/161.43545983636517,      'model':cell_convergence['TRN']['L6A']/cell_convergence['TRN']['VPM']},
                    }

    for pathway_ratio_ind,pathway_ratio in enumerate(pathway_ratios.keys()):
        for source in pathway_ratios[pathway_ratio].keys():
            width = 0.25
            if source == 'model':
                c   = 'green'
                loc =  pathway_ratio_ind+(width/2)
            elif source == 'BBP':
                c   = 'purple'
                loc =  pathway_ratio_ind-(width/2)

            plt.bar(loc,pathway_ratios[pathway_ratio][source],width=width,color=c)

    plt.xticks([0,1,2,3],['L6A/MLe ->VPM','TRN/MLe ->VPM','L6A-> VPM/TRN','L6A/VPM ->TRN'],rotation=10)
    plt.xlabel('Pathway')
    plt.ylabel('Conn Ratio')
    plt.legend(['Thal_NB model + literature','NetPyNE model'])

    plt.savefig('test_bar_plot.png',dpi=500)

    ############################################################################################################

    from    matplotlib  import  pyplot  as plt
    # print("Matplotlib backend (default): %s" %plt.get_backend())
    # modules = []
    # for module in sys.modules:
    #     if module.startswith('matplotlib'):
    #         modules.append(module)
    # for module in modules:
    #     sys.modules.pop(module)
    # import matplotlib
    # matplotlib.use("MacOSX")
    # from    matplotlib  import  pyplot  as plt
    # print("Matplotlib backend (dynamic): %s" %plt.get_backend())

    import matplotlib; matplotlib.use('agg')  # to avoid graphics error in servers

    plotConns=True
    selectGids=[419,619,800,850,1071]
    if plotConns:
        for cell_gid in selectGids: 
            fig = plt.figure()
            ax = fig.add_subplot(projection='3d')
            # --- Plot Soma center - SWC files are loaded with cells centered at (0,0,0)
            ax.scatter(0,0,0, marker='o',color='cyan',linewidths=0, alpha=1)
            
            # --- Plot 3d locations of cell morphology sections
            x3d=[];     y3d=[];     z3d=[]
            for sec in sim.net.cells[cell_gid].secs.keys():
                if 'pt3d' not in sim.net.cells[cell_gid].secs[sec]['geom'].keys():continue
                if len(sim.net.cells[cell_gid].secs[sec]['geom']['pt3d'])>0 and type(sim.net.cells[cell_gid].secs[sec]['geom']['pt3d'])==list:
                    for pt3d in sim.net.cells[cell_gid].secs[sec]['geom']['pt3d']:
                        x3d.append(pt3d[0]);    y3d.append(pt3d[1]);    z3d.append(pt3d[2])
            # # --- sec/pt3d points
            ax.scatter(x3d, y3d, z3d, marker='o',color='k',s=0.25,linewidths=0, alpha=0.1)
            
            
            # syn_keys = list(sim.net.params.synMechParams.keys())
            # for syn_key in sim.net.params.synMechParams.keys():
            
            conn_dict = {}
            colors = ['b','r','g','lightblue','orange','magenta','limegreen','grey']
            for ind,conn_key in enumerate(['conn|TRN|TRN','conn|TRN|VPM','conn|VPM|TRN','conn|MLe|VPM','conn|L6A|TRN','conn|L6A|VPM','conn|VPM|L6A']):
                # try:    conn_key = 'conn|'+syn_key.split('|')[1]+'|'+syn_key.split('|')[2]
                # except: continue
                conn_dict.update({conn_key:{'c':colors[ind],'syn':[]}})

            for conn_ind in range(len(sim.net.cells[cell_gid].conns)):
                
                postSec = sim.net.cells[cell_gid].conns[conn_ind]['sec']
                postLoc = sim.net.cells[cell_gid].conns[conn_ind]['loc']
                preSyn  = sim.net.cells[cell_gid].conns[conn_ind]['label']
                
                
                p0=sim.net.cells[cell_gid].secs[postSec]['geom']['pt3d'][0]
                p1=sim.net.cells[cell_gid].secs[postSec]['geom']['pt3d'][-1]
                p_x = (p1[0]-p0[0])*postLoc
                p_y = (p1[1]-p0[1])*postLoc
                p_z = (p1[2]-p0[2])*postLoc
                
                conn_dict[preSyn]['syn'].append((p_x,p_y,p_z))

            for syn_ind,syn_key in enumerate(conn_dict.keys()):
                try:
                    x3d = [x[0] for x in conn_dict[syn_key]['syn']]
                    y3d = [x[1] for x in conn_dict[syn_key]['syn']]
                    z3d = [x[2] for x in conn_dict[syn_key]['syn']]
                except:continue

                ax.scatter(x3d, y3d, z3d, marker='o',color=conn_dict[syn_key]['c'],s=1,linewidths=0, alpha=1)

            plt.savefig('./figs/cell_projections_'+str(cell_gid)+'.png',dpi=1000)

    plt.show()