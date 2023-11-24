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

try:
    from __main__ import cfg  # import SimConfig object with params from parent module
except:
    from barreloid_cfg import cfg

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

if cfg.loadCircuitProperties:
    NetPyNE_BBP.Prompt.headerMsg('Loading converted circuit from pickle file in \t\t'+cfg.stored_circuit_path)
    # with open(cfg.stored_circuit_path, 'rb') as f: circuit_dict = pickle.load(f)
    with open(cfg.stored_circuit_path, 'rb') as f: circuit_dict = pd.read_pickle(f)
else:
    NetPyNE_BBP.Prompt.headerMsg('Loading circuit from original project')
    circuit_dict = NetPyNE_BBP.LoadBBPCircuit.getDataFrames_new( cfg_file=cfg.sonataConfigFile)
    if cfg.saveCircuitProperties: 
        with open(cfg.stored_circuit_path, 'wb') as f: pd.to_pickle(circuit_dict, f)
        # with open(cfg.stored_circuit_path, 'wb') as f: pickle.dump(circuit_dict, f)

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# --- Cell Params
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
# --- Load cell models from converted NetPyNE morpohologies
loading_failed=[]
if cfg.loadCellModel:
    NetPyNE_BBP.Prompt.headerMsg('Loading stored cells from \t\t'+cfg.NetPyNE_JSON_cells)
    for thal_gid in cfg.select_thal_gids:
        cfg.convertCellModel=False
        cell_properties = NetPyNE_BBP.LoadBBPCircuit.getNodeProperties(circuit_dict,thal_gid)
        # NetPyNE_BBP.Prompt.headerMsg('Loading cell '+str(thal_gid))
        try:        netParams.loadCellParams(label=cell_properties['mtype']+'__'+str(thal_gid),fileName=cfg.NetPyNE_JSON_cells+'/netpyne_'+str(thal_gid)+'.json')
        except:     loading_failed.append(thal_gid);print('Loading failed - Cell '+str(thal_gid)+' will be created in the next step')
if len(loading_failed)>0:cfg.select_thal_gids=loading_failed;cfg.convertCellModel=True

# --- Load cell models from morpohology templates
if cfg.convertCellModel:
    NetPyNE_BBP.Prompt.headerMsg('Processing cells')
    for thal_gid_ind,thal_gid in enumerate(cfg.select_thal_gids):
        print('\t>>\t',str(thal_gid),'\t|\t',str(len(cfg.select_thal_gids)-thal_gid_ind),' cells left')
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

for gids_list in cfg.gids_lists:
    # --- Adding cell diversity rule
    for thal_gid in gids_list:
        cell_properties = circuit_dict['thalamus_neurons'].loc[thal_gid]
        # --- Adds each GID morphology into the cell diversity dictionary
        # netParams.cellParams[cell_properties['mtype']+'__'+str(thal_gid)].update({'diversityFraction':1/len(gids_list)})
        netParams.cellParams[cell_properties['mtype']+'__'+str(thal_gid)].update({'diversityFraction':1/len(list(set(gids_list)))})
        # --- Changes cellType to the same value, so that the different morphologies are added to the same pop can be identified by connParams
        netParams.cellParams[cell_properties['mtype']+'__'+str(thal_gid)]['conds']['cellType']=cell_properties['mtype']

# --- Creates a single pop for all cells and defines popParams with diversity rule so that each cell has a different morphology, based on its GIDs. This aproach makes it easier to target stims and run parallel sims
clone_pops = ['VPM__pop','TRN__pop']
pop_vars   = ['xRange','yRange','zRange','xnormRange','ynormRange','znormRange','numCells']

netParams.popParams['VPM__pop']={'cellType':'VPL_TC','diversity': True}
netParams.popParams['TRN__pop']={'cellType':'Rt_RC', 'diversity': True}

for pop_name in clone_pops:
    for pop_var in pop_vars: netParams.popParams[pop_name].update({pop_var:cfg.store_network['pops'][pop_name]['tags'][pop_var]})

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# --- L6A cells and pop
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
import Build_Net as BN

if cfg.simplifyL6A: 
    netParams.cellParams.update(BN.BuildNetwork.getCellTemplate(    template='izhi',pops=['L6A'],cell_flag='__cell'))
    netParams.popParams.update( BN.BuildNetwork.getPopTemplate(     pops=['L6A'],center_point=cfg.center_point,pop_flag='__pop',cell_flag='__cell'))
    secTarget = 'soma_0'    # single section
else:
    netParams.cellParams.update(BN.BuildNetwork.getL6ACellTemplate( cellsFolder=cfg.NetPyNE_L6A_JSON_cells,file_format='pkl'))
    # netParams.cellParams.update(BN.BuildNetwork.getL6ACellTemplate( cellsFolder=cfg.NetPyNE_L6A_JSON_cells))
    netParams.popParams.update( BN.BuildNetwork.getPopTemplate(     pops=['L6A'],center_point=cfg.center_point,pop_flag='__pop',cell_flag='__cell',diversity=True,volume_shape='cube'))
    secTarget = 'basal'     # secList
    print('\t>>\tWarning: Add code to remove AXON segment from L6A cells')

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# --- MLe cells and pops
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
align_cells=False
useVecStim = True
if useVecStim:
    from GenerateStimDataset import SampleData
    try:
        print('Loading spikes from \n',cfg.deflection_dataset_path)
        # spikes_dict = SampleData.LoadJSON(file_path=cfg.deflection_dataset_path)
        # --- Load spikes dictionary (stored in Pickle for stora - large files)
        spikes_dict = SampleData.LoadPickle(file_path=cfg.deflection_dataset_path)

    except:
        cfg.deflection_dataset_path = cfg.deflection_model_folder+'spike_dicts/mleSpikes_deflectionAngle|0_samplingBins|1000us_simplifiedDataset.pkl'
        print('Raster Loading failed - Loading sample raster: deflection 0 degrees | sampling interval 1000 us')
        # spikes_dict = SampleData.LoadJSON(file_path=cfg.deflection_dataset_path)
        # --- Load spikes dictionary (stored in Pickle for stora - large files)
        spikes_dict = SampleData.LoadPickle(file_path=cfg.deflection_dataset_path)


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
if   cfg.modType == 'Det':              modAMPANMDA = 'DetAMPANMDA';                modGABA     = 'DetGABAAB'              # S1 Deterministic  implementation of the BBP mod files
elif cfg.modType == 'Prob_S1':          modAMPANMDA = 'ProbAMPANMDA_EMS_S1';        modGABA     = 'ProbGABAAB_EMS_S1'      # S1 Probabilistic  implementation of the BBP mod files
elif cfg.modType == 'Prob_original':    modAMPANMDA = 'ProbAMPANMDA_EMS_original';  modGABA     = 'ProbGABAAB_EMS_original' # original MOD from BBP model
else:                                   modAMPANMDA = 'ProbAMPANMDA_EMS';           modGABA     = 'ProbGABAAB_EMS'         # Original Thalamus implementation of the BBP mod files
print('\n\t>>\tMOD template\tAMPA: ',   modAMPANMDA, '\tGABA: ', modGABA)

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
# --- Model adjustments for in vitro condition
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# --- Changes extracellular Ca2+ concentration for all sections in the biophysical cell models
if cfg.cao_secs is not None:
    print('\t>>\tChanging extracellular Ca2+ concentration to ', str(cfg.cao_secs))
    for biophys_cell in netParams.cellParams.keys():
        for sec in netParams.cellParams[biophys_cell]['secs'].keys():
            if 'ions' in netParams.cellParams[biophys_cell]['secs'][sec].keys():
                if 'ca' in netParams.cellParams[biophys_cell]['secs'][sec]['ions'].keys(): netParams.cellParams[biophys_cell]['secs'][sec]['ions']['ca']['o'] = cfg.cao_secs

# --- Rescale USE parameter (probability of synapse activation)
if cfg.rescaleUSE is not None:
    print('\t>>\tRescaling synaptic USE to ', str(cfg.rescaleUSE))
    for mech in netParams.synMechParams.keys():
        try:    netParams.synMechParams[mech]['Use']*=cfg.rescaleUSE
        except: continue

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# --- Connectivity
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# --- MLe to VPM connections

#       Obs: Projections are organized so that:
#               the top    of the MLe projects to the top    of VPM (absolute top = 0, pia surface)
#               the bottom of the MLe projects to the bottom of VPM
#            Based on the data from (Timofeeva et al., 2003)

if cfg.conn_MLe_VPM:
    # syns per conn - obs: conn_data is organized as conn_data[post_pop]['chem'][pre_pop][synsPerConn/convergence][MEAN,STD]
    theta_pops = [pop_name for pop_name in netParams.popParams.keys() if 'MLe' in pop_name]
    if useVecStim:
        print('\t>>\tAdding VecStim MLe connections - One rule per (cell/pop)')
        for pre_pop in theta_pops:
            mle_pop = pre_pop.split('__')[0]
            mle_degree = int(mle_pop.split('@')[1])
            # --- Conditions to test the organization of MLe to VPM projections
            # if mle_degree>20:continue
            # if mle_degree<340:continue
            (mle_conn_method, mle_conn_rule)=BN.BuildNetwork.laminarProjection_VecStim(pre_angle=mle_degree)
            conn_dict = BN.BuildNetwork.getConnDict(pre_pop         = mle_pop,
                                                    post_pop        = 'VPM',
                                                    conn_method     = mle_conn_method,
                                                    conn_rule       = mle_conn_rule,
                                                    syn_mech        = cfg.MLe_VPM__synMech,
                                                    syns_per_conn   = cfg.MLe_VPM__syns_per_conn,
                                                    conn_type       = cfg.MLe_VPM__conn_type,
                                                    weight          = cfg.MLe_VPM__weight,
                                                    target_secs     = cfg.MLe_VPM__secList)
            netParams.connParams.update(conn_dict)
    else:
        (mle_conn_method, mle_conn_rule)=BN.BuildNetwork.laminarProjection()
        for pre_pop in theta_pops:
            mle_pop = pre_pop.split('__')[0]
            conn_dict = BN.BuildNetwork.getConnDict(pre_pop         = mle_pop,
                                                    post_pop        = 'VPM',
                                                    conn_method     = mle_conn_method,
                                                    conn_rule       = mle_conn_rule,
                                                    syn_mech        = cfg.MLe_VPM__synMech,
                                                    syns_per_conn   = cfg.MLe_VPM__syns_per_conn,
                                                    conn_type       = cfg.MLe_VPM__conn_type,
                                                    weight          = cfg.MLe_VPM__weight,
                                                    target_secs     = cfg.MLe_VPM__secList)
            netParams.connParams.update(conn_dict)

# scale_prob = 40/12.483479074117373

if cfg.conn_L6A_VPM:
    conn_method_1D, conn_rule_1D = BN.BuildNetwork.distanceBasedProbability_1D_exponential(pre_pop      = cfg.L6A_VPM__pre_pop,     post_pop     = cfg.L6A_VPM__post_pop,
                                                                                           pre_pop_axis = cfg.L6A_VPM__pre_pop_axis,post_pop_axis= cfg.L6A_VPM__post_pop_axis,
                                                                                           conn_prob    = cfg.L6A_VPM__conn_prob,   center_point = cfg.center_point)
    conn_dict_1D = BN.BuildNetwork.getConnDict( pre_pop         = cfg.L6A_VPM__pre_pop,
                                                post_pop        = cfg.L6A_VPM__post_pop,
                                                conn_method     = conn_method_1D,
                                                conn_rule       = conn_rule_1D,
                                                syn_mech        = cfg.L6A_VPM__syn_mech,
                                                syns_per_conn   = cfg.L6A_VPM__syns_per_conn,
                                                conn_type       = cfg.L6A_VPM__conn_type,
                                                weight          = cfg.L6A_VPM__weight,
                                                target_secs     = cfg.L6A_VPM__target_secs)
    netParams.connParams.update(conn_dict_1D)

if cfg.conn_L6A_TRN:
    conn_method_1D, conn_rule_1D = BN.BuildNetwork.distanceBasedProbability_1D_exponential(pre_pop      = cfg.L6A_TRN__pre_pop,     post_pop     = cfg.L6A_TRN__post_pop,
                                                                                           pre_pop_axis = cfg.L6A_TRN__pre_pop_axis,post_pop_axis= cfg.L6A_TRN__post_pop_axis,
                                                                                           conn_prob    = cfg.L6A_TRN__conn_prob,   center_point = cfg.center_point)
    conn_dict_1D = BN.BuildNetwork.getConnDict( pre_pop         = cfg.L6A_TRN__pre_pop,
                                                post_pop        = cfg.L6A_TRN__post_pop,
                                                conn_method     = conn_method_1D,
                                                conn_rule       = conn_rule_1D,
                                                syn_mech        = cfg.L6A_TRN__syn_mech,
                                                syns_per_conn   = cfg.L6A_TRN__syns_per_conn,
                                                conn_type       = cfg.L6A_TRN__conn_type,
                                                weight          = cfg.L6A_TRN__weight,
                                                target_secs     = cfg.L6A_TRN__target_secs)
    netParams.connParams.update(conn_dict_1D)

if cfg.conn_TRN_VPM:
    conn_dict_1D = BN.BuildNetwork.getConnDict( pre_pop         = cfg.TRN_VPM__pre_pop,
                                                post_pop        = cfg.TRN_VPM__post_pop,
                                                conn_method     = cfg.TRN_VPM__conn_method,
                                                conn_rule       = cfg.TRN_VPM__conn_rule,
                                                syn_mech        = cfg.TRN_VPM__syn_mech,
                                                syns_per_conn   = cfg.TRN_VPM__syns_per_conn,
                                                conn_type       = cfg.TRN_VPM__conn_type,
                                                weight          = cfg.TRN_VPM__weight,
                                                target_secs     = cfg.TRN_VPM__target_secs)
    netParams.connParams.update(conn_dict_1D)

if cfg.conn_TRN_TRN:
    conn_dict_1D = BN.BuildNetwork.getConnDict( pre_pop         = cfg.TRN_TRN__pre_pop,
                                                post_pop        = cfg.TRN_TRN__post_pop,
                                                conn_method     = cfg.TRN_TRN__conn_method,
                                                conn_rule       = cfg.TRN_TRN__conn_rule,
                                                syn_mech        = cfg.TRN_TRN__syn_mech,
                                                syns_per_conn   = cfg.TRN_TRN__syns_per_conn,
                                                conn_type       = cfg.TRN_TRN__conn_type,
                                                weight          = cfg.TRN_TRN__weight,
                                                target_secs     = cfg.TRN_TRN__target_secs)
    netParams.connParams.update(conn_dict_1D)

if cfg.conn_TRNe_TRNe:
    conn_dict_1D = BN.BuildNetwork.getConnDict( pre_pop         = cfg.TRNe_TRNe__pre_pop,
                                                post_pop        = cfg.TRNe_TRNe__post_pop,
                                                conn_method     = cfg.TRNe_TRNe__conn_method,
                                                conn_rule       = cfg.TRNe_TRNe__conn_rule,
                                                syn_mech        = cfg.TRNe_TRNe__syn_mech,
                                                syns_per_conn   = cfg.TRNe_TRNe__syns_per_conn,
                                                conn_type       = cfg.TRNe_TRNe__conn_type,
                                                weight          = cfg.TRNe_TRNe__weight,
                                                target_secs     = cfg.TRNe_TRNe__target_secs)
    netParams.connParams.update(conn_dict_1D)

if cfg.conn_VPM_TRN:
    conn_method_1D, conn_rule_1D = BN.BuildNetwork.distanceBasedProbability_1D_exponential(pre_pop      = cfg.VPM_TRN__pre_pop,     post_pop     = cfg.VPM_TRN__post_pop,
                                                                                           pre_pop_axis = cfg.VPM_TRN__pre_pop_axis,post_pop_axis= cfg.VPM_TRN__post_pop_axis,
                                                                                           conn_prob    = cfg.VPM_TRN__conn_prob,   center_point = cfg.center_point)
    conn_dict_1D = BN.BuildNetwork.getConnDict( pre_pop         = cfg.VPM_TRN__pre_pop,
                                                post_pop        = cfg.VPM_TRN__post_pop,
                                                conn_method     = conn_method_1D,
                                                conn_rule       = conn_rule_1D,
                                                syn_mech        = cfg.VPM_TRN__syn_mech,
                                                syns_per_conn   = cfg.VPM_TRN__syns_per_conn,
                                                conn_type       = cfg.VPM_TRN__conn_type,
                                                weight          = cfg.VPM_TRN__weight,
                                                target_secs     = cfg.VPM_TRN__target_secs)
    netParams.connParams.update(conn_dict_1D)

if cfg.conn_VPM_L6A:
    # conn_method_, conn_rule_ = BN.BuildNetwork.positionBasedProbability_1D_linear(  pre_pop     = cfg.VPM_L6A__pre_pop,     post_pop     = cfg.VPM_L6A__post_pop,
    #                                                                                 pre_pop_axis= cfg.VPM_L6A__pre_pop_axis,post_pop_axis= cfg.VPM_L6A__post_pop_axis,
    #                                                                                 conn_prob   = cfg.VPM_L6A__conn_prob,   coef         = cfg.VPM_L6A__coef, 
    #                                                                                 center_point= cfg.center_point)
    
    # conn_dict_ = BN.BuildNetwork.getConnDict(   pre_pop         = cfg.VPM_L6A__pre_pop,
    #                                             post_pop        = cfg.VPM_L6A__post_pop,
    #                                             conn_method     = conn_method_,
    #                                             conn_rule       = conn_rule_,
    #                                             syn_mech        = cfg.VPM_L6A__syn_mech,
    #                                             syns_per_conn   = cfg.VPM_L6A__TC_syns_per_conn,
    #                                             conn_type       = cfg.VPM_L6A__conn_type,
    #                                             weight          = cfg.VPM_L6A__weight,
    #                                             target_secs     = secTarget)
    # netParams.connParams.update(conn_dict_)
    
    # --- TEST: 2023_11_09 - removing linear scaling of conn prob because the axis don't match - using distance-based conn instead
    conn_method_1D, conn_rule_1D = BN.BuildNetwork.distanceBasedProbability_1D_exponential(pre_pop      = cfg.VPM_L6A__pre_pop,     post_pop     = cfg.VPM_L6A__post_pop,
                                                                                           pre_pop_axis = cfg.VPM_L6A__pre_pop_axis,post_pop_axis= cfg.VPM_L6A__post_pop_axis,
                                                                                           conn_prob    = cfg.VPM_L6A__conn_prob,   center_point = cfg.center_point)
    conn_dict_1D = BN.BuildNetwork.getConnDict( pre_pop         = cfg.VPM_L6A__pre_pop,
                                                post_pop        = cfg.VPM_L6A__post_pop,
                                                conn_method     = conn_method_1D,
                                                conn_rule       = conn_rule_1D,
                                                syn_mech        = cfg.VPM_L6A__syn_mech,
                                                syns_per_conn   = cfg.VPM_L6A__TC_syns_per_conn,
                                                conn_type       = cfg.VPM_L6A__conn_type,
                                                weight          = cfg.VPM_L6A__weight,
                                                target_secs     = secTarget)
    netParams.connParams.update(conn_dict_1D)

    
    # 133*666/
    # (318*9)

    # see (Meyer 2010)



# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# --- Stimulation
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# --- Adds a current stim to thalamic populations
if cfg.add_current_stims:
    NetPyNE_BBP.Prompt.headerMsg('Adding Current stim (IClamp)')
    
    for pop_ind, pop in enumerate(netParams.popParams.keys()):
        if pop in cfg.current_stim_targets:
            netParams.stimSourceParams['IClamp_'+str(pop_ind)] = {'type': 'IClamp', 'del': cfg.current_stim_start, 'dur': cfg.duration, 'amp': cfg.current_stim_amp}
            netParams.stimTargetParams['IClamp_'+str(pop_ind)+'__'+pop] = {'source': 'IClamp_'+str(pop_ind), 'sec':'soma_0', 'loc': 0.5, 'conds': {'pop':pop}}

if cfg.add_bkg_stim:
    for ind in range(len(cfg.bkg_rate)):
        netParams.stimSourceParams['bkg_'+str(ind)] = { 
            'type':     'NetStim', 
            'rate':     cfg.bkg_rate[ind], 
            'noise':    cfg.bkg_noise[ind],
            }
        netParams.stimTargetParams['bkg_'+str(ind)+'|'+cfg.bkg_pop[ind].split('__')[0]+'__pop'] = { 
            'source':   'bkg_'+str(ind), 
            'conds':    {'cellType': cfg.bkg_pop[ind]}, 
            'weight':   cfg.bkg_weight[ind], 
            'delay':    cfg.bkg_delay[ind], 
            'synMech':  cfg.bkg_synMech[ind],
            }











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
