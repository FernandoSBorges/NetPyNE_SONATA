'''
netParams.py 

Contributors: joao.moreira@downstate.edu, salvadordura@gmail.com
'''

import NetPyNE_BBP
import numpy as np
from netpyne import specs
import pickle, json
import pandas as pd
import os


netParams = specs.NetParams()   # object of class NetParams to store the network parameters

try:
    from __main__ import cfg  # import SimConfig object with params from parent module
except:
    from cfg import cfg

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
netParams.defaultThreshold = -25.0  # spike threshold, 10 mV is NetCon default, lower it for all cells

### reevaluate these values
netParams.defaultDelay = 2.0 # default conn delay (ms) # DEFAULT
netParams.propVelocity = 500.0 # propagation velocity (um/ms)

netParams.defineCellShapes = True # JV 2021-02-23 - Added to fix the lack of the pt3d term in the cells, which make it unable to record i_membrane_

#------------------------------------------------------------------------------
# --- Load BBP circuit
#------------------------------------------------------------------------------
if cfg.convertCellMorphologies: 
    NetPyNE_BBP.Conversion.convertCellMorphology(inputFolder_h5=cfg.morphologyFolder_h5,
                                                 outputFolder_swc=cfg.NetPyNE_exportedCells,inputFolder_asc=cfg.morphologyFolder_asc)

circuit_dict    = NetPyNE_BBP.LoadBBPCircuit.getDataFrames_new(cfg_file=cfg.sonataConfigFile)

thal_gids = cfg.select_thal_gids

loading_failed=[]
for thal_gid in thal_gids:
    cell_properties = NetPyNE_BBP.LoadBBPCircuit.getNodeProperties(circuit_dict,thal_gid)
    if cfg.loadCellModel:
        cfg.convertCellModel=False
        NetPyNE_BBP.Prompt.headerMsg('Loading cell '+str(thal_gid))
        try:        
            netParams.loadCellParams(label=cell_properties['mtype']+'__'+str(thal_gid),
                                     fileName=cfg.NetPyNE_JSON_cells+'/netpyne_'+str(thal_gid)+'.json')
        except:     
            loading_failed.append(thal_gid);
            print('Loading failed - Cell '+str(thal_gid)+' will be created in the next step')

if len(loading_failed)>0:
    thal_gids=loading_failed;
    cfg.convertCellModel=True

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# --- Cell Params
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
NetPyNE_BBP.Prompt.headerMsg('Processing cells')

for thal_gid_ind,thal_gid in enumerate(thal_gids):

    print('\t>>\t',str(thal_gid),'\t|\t',str(len(thal_gids)-thal_gid_ind),' cells left')

    cell_properties = NetPyNE_BBP.LoadBBPCircuit.getNodeProperties(circuit_dict,thal_gid)

    if cfg.convertCellModel:
        thal_gid_conds={}

        thal_gid_conds.update({'cell_index':thal_gid})

        netParams.importCellParams(
            label           = cell_properties['mtype']+'__'+str(thal_gid), 
            fileName        = cfg.NetPyNE_templateCells+'/'+cell_properties['etype']+'.hoc', 
            cellName        = cell_properties['etype'], 
            conds           = {'cellType':str(thal_gid)},
            cellArgs        = [thal_gid,cfg.NetPyNE_exportedCells,cell_properties['morphology']+'.swc'], 
            importSynMechs  = True, 
            somaAtOrigin    = False, 
            cellInstance    = False,
        )
        if cfg.saveCellModel: netParams.saveCellParams(label=cell_properties['mtype']+'__'+str(thal_gid),fileName=cfg.NetPyNE_JSON_cells+'/netpyne_'+str(thal_gid)+'.json')
        
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
# --- Pop Params
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
NetPyNE_BBP.Prompt.headerMsg('Processing pops and conns')

if cfg.connType == 'original':
    # --- Declaring variables that store information at the network level
    store_edges_dict={}
    # --- Defining variables to point to the noise files
    noise_filePaths=[('CorticoThalamic_projections',cfg.ct_virtual_noise),('MedialLemniscus_projections',cfg.ml_virtual_noise)]
    
    # --- Adding cell diversity rule
    for thal_gid in thal_gids:
        cell_properties = circuit_dict['thalamus_neurons'].loc[thal_gid]
        # --- Adds each GID morphology into the cell diversity dictionary
        # --- Creates a separate pop for each cell
        # if not cfg.th_singleTcPop: netParams.popParams[cell_properties['mtype']+'|'+str(thal_gid)+'__pop']={'cellType':cell_properties['mtype'], 'numCells': 1}
        if cfg.th_singleTcPop:  netParams.cellParams[cell_properties['mtype']+'__pop'].update({'diversityFraction':1/len(thal_gids)})
        else:                   netParams.popParams[cell_properties['mtype']+'__'+str(thal_gid)+'__pop']={'cellType':str(thal_gid), 'numCells': 1}

    # --- Creates a single pop for all cells and defines popParams with diversity rule so that each cell has a different morphology, based on its GIDs. This aproach makes it easier to target stims and run parallel sims
    if cfg.th_singleTcPop: netParams.popParams[cell_properties['mtype']+'__pop']={'cellType':cell_properties['mtype'], 'numCells': len(thal_gids),'diversity': True}

    # --- Changes extracellular Ca2+ concentration for all sections in the biophysical cell models
    if cfg.cao_secs is not None:
        print('\t>>\tChanging extracellular Ca2+ concentration to ', str(cfg.cao_secs))
        for biophys_cell in netParams.cellParams.keys():
            for sec in netParams.cellParams[biophys_cell]['secs'].keys():
                if 'ions' in netParams.cellParams[biophys_cell]['secs'][sec].keys():
                    if 'ca' in netParams.cellParams[biophys_cell]['secs'][sec]['ions'].keys(): netParams.cellParams[biophys_cell]['secs'][sec]['ions']['ca']['o'] = cfg.cao_secs

    # --- Fixes Rt_RC cell biophysics by removing ih and ia currents
    if (cfg.removeExtraCurrents) and (cfg.target == 'Rt_RC'):
        remove_currents=[('TC_iA','gk_max'),('TC_ih_Bud97','gh_max')]
        for biophys_cell in netParams.cellParams.keys():
            for sec in netParams.cellParams[biophys_cell]['secs'].keys():
                for remove_current in remove_currents:
                    if remove_current[0] in netParams.cellParams[biophys_cell]['secs'][sec]['mechs'].keys(): netParams.cellParams[biophys_cell]['secs'][sec]['mechs'][remove_current[0]][remove_current[1]]=0

    # --- Change the parameters of a selected current/channel to a selected value
    if len(cfg.changeCurrents)>0:
        for biophys_cell in netParams.cellParams.keys():
            for sec in netParams.cellParams[biophys_cell]['secs'].keys():
                # --- Change_current template: cfg.changeCurrents = (current_name,current_variable,new_value)
                for change_current in cfg.changeCurrents:
                    if change_current[0] in netParams.cellParams[biophys_cell]['secs'][sec]['mechs'].keys(): netParams.cellParams[biophys_cell]['secs'][sec]['mechs'][change_current[0]][change_current[1]]=change_current[2]
    
    # --- Rescale the parameters of a selected current/channel to a multiplier value (e.g.: 0.5, 2.0, 4.0)
    if len(cfg.rescaleCurrents)>0:
        for biophys_cell in netParams.cellParams.keys():
            for sec in netParams.cellParams[biophys_cell]['secs'].keys():
                # --- Rescale current template: cfg.rescaleCurrents = (current_name,current_variable,new_value)
                for rescale_current in cfg.rescaleCurrents:
                    if rescale_current[0] in netParams.cellParams[biophys_cell]['secs'][sec]['mechs'].keys(): netParams.cellParams[biophys_cell]['secs'][sec]['mechs'][rescale_current[0]][rescale_current[1]]*=rescale_current[2]

    # --- Dictionary to store unique references to the presynaptic cells
    edge_sources_dict_temp={}
    for pathway in cfg.th_select_pathways: edge_sources_dict_temp.update({pathway:[]})

    NetPyNE_BBP.Prompt.headerMsg('Loading conns')
    for thal_gid_ind,thal_gid in enumerate(thal_gids):
        print('\n\t>>\t',str(thal_gid),'\t|\t',str(len(thal_gids)-thal_gid_ind),' cells left')

        cell_properties = NetPyNE_BBP.LoadBBPCircuit.getNodeProperties(circuit_dict,thal_gid)
        
        # --- Loads the cell directly from the SWC file to recreate the indexing done using 'afferent_section_pos'
        morphologyFullPath = cfg.NetPyNE_exportedCells+'/'+cell_properties['morphology']+'.swc'
        full_ids_list = NetPyNE_BBP.Conversion.getAfferentSectionIds(morphologyFullPath)
        
        if   cfg.select_microcircuit is None:   cell_conns_file = str(thal_gid)+'_edges_full.json'
        elif cfg.select_microcircuit == int:    cell_conns_file = str(thal_gid)+'_edges_mc'+str(cfg.select_microcircuit)+'.json'
            
        edges_fileName  = cfg.NetPyNE_node_pathway_edges+'/'+cell_conns_file
        
        # --- Overwrites cfg.th_updateConns if conn file does not exits
        if (cfg.th_updateConns) or (cell_conns_file not in os.listdir(cfg.NetPyNE_node_pathway_edges)):
            
            if (cell_conns_file not in os.listdir(cfg.NetPyNE_node_pathway_edges)): print('\t\t--\tFile does not exist... Creating conns for cell '+str(thal_gid))
            else: print('\t\t--\tUpdating conns for cell '+str(thal_gid))

            # --- Generates a Dictionary of Dataframes with the EDGES (or conns) information for each pathway targeting the selected NODE
            #     Read mode about edges: https://github.com/AllenInstitute/sonata/blob/master/docs/SONATA_DEVELOPER_GUIDE.md
            node_pathway_edges = NetPyNE_BBP.LoadBBPCircuit.getNodeEdges(       cfg.sonataConfigFile,
                                                                                thal_gid,
                                                                                special_conditions={'inter_sources':cfg.th_inter_sources,'inter_targets':[],'select_microcircuit':cfg.select_microcircuit}
                                                                                )
            # --- Adds section information in NetPyNE format to the edge properties
            node_pathway_edges = NetPyNE_BBP.ConvertSynapses.getTargetSections( cell                = netParams.cellParams[cell_properties['mtype']+'__'+str(thal_gid)],
                                                                                node_pathway_edges  = node_pathway_edges,
                                                                                full_ids_list       = full_ids_list
                                                                                )
            # --- Adds a remapped version of the 3d positions of each synapse, once the soma is set to (0,0,0) when a SWC morphology is loaded
            node_pathway_edges = NetPyNE_BBP.ConvertSynapses.convert3DLocation( cell                = netParams.cellParams[cell_properties['mtype']+'__'+str(thal_gid)],
                                                                                cell_properties     = cell_properties,
                                                                                node_pathway_edges  = node_pathway_edges,
                                                                                )
            # --- Converting conns from DICT of DFs to DICT of DICTs and writing to JSON
            node_pathway_edges_dict = NetPyNE_BBP.LoadBBPCircuit.saveNodeEdges(node_pathway_edges,edges_fileName)
        else:
            print('\t\t--\tLoading conns/edges from stored dataset')
            # --- Loads the connectivity information previously saved
            node_pathway_edges = NetPyNE_BBP.LoadBBPCircuit.loadNodeEdges(edges_fileName)

        if cfg.plotSynLocation:
            # --- Plots the location each synapse in the cell model
            NetPyNE_BBP.ConvertSynapses.plotSynapseLocation(    cell                = netParams.cellParams[cell_properties['mtype']+'__'+str(thal_gid)],
                                                                cell_properties     = cell_properties,
                                                                node_pathway_edges  = node_pathway_edges,
                                                                full_ids_list       = full_ids_list,
                                                                extra_flag          = '',
                                                                plotSimplified      = False
                                                                )
        if cfg.plotCellShape:
            NetPyNE_BBP.ConvertSynapses.plotCellShape(          cell=netParams.cellParams[cell_properties['mtype']+'__'+str(thal_gid)],
                                                                extra_flag='_cellShape'
                                                                )


        # --- Creates alternative edges and plots the original and reordered locations each synapse in the cell model
        if cfg.th_topological:
            node_pathway_edges_r = \
            NetPyNE_BBP.ConvertSynapses.reorderSynapseLocation( cell                = netParams.cellParams[cell_properties['mtype']+'__'+str(thal_gid)],
                                                                cell_properties     = cell_properties,
                                                                node_pathway_edges  = node_pathway_edges,
                                                                full_ids_list       = full_ids_list,
                                                                plotFigure          = True
                                                                )
        
        # --- Selects the dataset of edges distribution to be used (INESPECIFIC vs TOPOLOGICAL)
        if cfg.th_topological:  edges = node_pathway_edges_r
        else:                   edges = node_pathway_edges

        # --- Selects which pathways should be included from the original model (e.g.: ct, ml, rt, in)
        if cfg.th_select_pathways is not None: edges = NetPyNE_BBP.LoadBBPCircuit.selectPathways(cfg.th_select_pathways,edges)

        # --- Switches the Synaptic parameters by the average values reported in the Table S2 in the original paper instead of the model stored values
        if cfg.th_useTableVals: edges = NetPyNE_BBP.ConvertSynapses.replaceSynParams(cell_properties,edges)

        # --- Creates a dictionary to store the presynaptic cells that project to the evaluated cell
        for edge_name in edges.keys():
            a,b,edge_source_pop,c = edge_name.split('__')
            edge_sources = list(set(edges[edge_name]['@source_node']));edge_sources.sort()
            edge_sources_dict_temp[edge_source_pop]+=edge_sources
            edge_sources_dict_temp[edge_source_pop].sort()

        # --- Saves a dictionary of noise input to each postsynaptic cell
        if cfg.saveIndividualNoise:
            for (key,noise_filePath) in noise_filePaths: 
                if key in edge_sources_dict_temp: 
                    cell_source_gids = list(set(edge_sources_dict_temp[key])); cell_source_gids.sort()
                    NetPyNE_BBP.LoadSimResults.getVirtualInputSpikes(noise_filePath, filterGids=cell_source_gids, saveName= cfg.NetPyNE_input_noise_savePath+'/post_cell_inputs/'+'tgt|'+str(thal_gid)+'__'+'src|'+key, loadFromSourceFile=False)


        # --- Stores the edge datasets so that it can be re-iterated over without overwriting
        for edge_name in edges.keys():
            store_edges_dict.update({edge_name:edges[edge_name]})

    # --- Combines the presynaptic cells into a single dictionary to create the Vectstim populations only with connected cells
    edge_sources_dict={}
    for edge_source_pop in edge_sources_dict_temp.keys():
        source_gids=list(set(edge_sources_dict_temp[edge_source_pop]))
        source_gids.sort()
        edge_sources_dict.update({edge_source_pop:source_gids})

    NetPyNE_BBP.Prompt.headerMsg('Loading internal inputs')
    # --- Loading dictionary of spike times
    th_spk_dict = NetPyNE_BBP.LoadSimResults.getSpikesFromSim(filePath=cfg.th_spikes_file, cfg_file=cfg.sonataConfigFile, microcircuit_number=cfg.select_microcircuit, showFig=False) # --- intrathalamic cells (biophysiscal)
    # --- Loading dictionary of spike times for connected cells only
    th_connected_cells_dict = NetPyNE_BBP.LoadSimResults.getSpikesFromConnectedCells(th_spk_dict, thal_gids, edge_sources_dict, showFig=False) # --- intrathalamic cells (biophysiscal)

    # # --- Loading and plotting voltage reports for the BBP model mc2 validation
    # th_soma_voltage_dict = NetPyNE_BBP.LoadSimResults.loadVoltageReport(filePath=cfg.Fig_4_1_traces, cfg_file=cfg.sonataConfigFile, gids=cfg.select_thal_gids, dt = 0.1, timeWindow = [0,4000], microcircuit_number=cfg.select_microcircuit, showFig=True)

    # for node_id in list(set(store_edges_dict['internal__chemical__thalamus_neurons|Rt_RC__38906']['@source_node'])):
    #     c_p = NetPyNE_BBP.LoadBBPCircuit.getNodeProperties(circuit_dict,node_id)
    #     print(c_p['region'],'\t',node_id, th_connected_cells_dict['Rt_RC'][node_id])
    

    # spk_nums={}; spk_count_dict={}
    # for t_cell in th_spk_dict[cell_properties['mtype']].keys(): spk_nums.update({t_cell:len(th_spk_dict[cell_properties['mtype']][t_cell])})
    # for spk_count in list(set(spk_nums.values())): spk_count_dict.update({spk_count:[]})
    # for t_cell in th_spk_dict[cell_properties['mtype']].keys(): spk_count_dict[len(th_spk_dict[cell_properties['mtype']][t_cell])].append(t_cell)
    # import matplotlib.pyplot as plt
    # plt.figure();plt.hist(spk_nums.values(),bins=50,log=True);plt.savefig('test_hist.png')

    NetPyNE_BBP.Prompt.headerMsg('Loading external inputs')
    # --- Loading dictionaries of input noise
    noise_dict={}
    for (key,noise_filePath) in noise_filePaths: 
        if key in edge_sources_dict: noise_dict.update({key:NetPyNE_BBP.LoadSimResults.getVirtualInputSpikes(noise_filePath, filterGids=edge_sources_dict[key], saveName= cfg.NetPyNE_input_noise_savePath+'/'+key, loadFromSourceFile=False)})

    # --- Creating VecStims pops --- one cell per pop (pops based on preCell bbp_gid)
    for edge_source_pop in edge_sources_dict.keys():
        source_nodes = edge_sources_dict[edge_source_pop]
        if ('CorticoThalamic' in edge_source_pop) or ('MedialLemniscus' in edge_source_pop): 
            spkids = list(noise_dict[edge_source_pop].keys())
            spkts  = list(noise_dict[edge_source_pop].values())
        elif 'thalamus_neurons' in edge_source_pop:
            edge_source_pop_=edge_source_pop.split('|')
            spkids = list(th_connected_cells_dict[edge_source_pop_[1]].keys())
            spkts  = list(th_connected_cells_dict[edge_source_pop_[1]].values())
            # spkids = list(th_spk_dict[edge_source_pop_[1]].keys())
            # spkts  = list(th_spk_dict[edge_source_pop_[1]].values())


        if cfg.skipTime is not None: 
            for spkid_ind, spkid in enumerate(spkids):spkts[spkid_ind] = [spkt+cfg.skipTime for spkt in spkts[spkid_ind]]


        # --- Creates NetPyNE Vecstim Pops --- one cell per pop
        for spkid_ind, spkid in enumerate(spkids):
            if int(spkid) in edge_sources_dict[edge_source_pop]:
                try:
                    if len(spkts[spkid_ind])<1: spkts[spkid_ind]=[cfg.duration+1000]
                except: 
                    if type(spkts[spkid_ind])==None: spkts[spkid_ind]=[cfg.duration+1000] # adds a dummy spike at t>simDuration to prevent code from breaking because there were no spikes at the cell

                netParams.popParams.update({ edge_source_pop+'__'+str(spkid)+'__vecstim_pop':{  'cellModel':'VecStim',
                                                                                                'numCells':1,
                                                                                                'spkTimes': spkts[spkid_ind],
                                                                                                'pulses':[]}})
    #             # --- test - add minis using <netstim_inhpoisson.mod>
    #             netParams.popParams.update({ edge_source_pop+'__'+str(spkid)+'__minis_pop':{    'cellModel': 'NetStim', 
    #                                                                                             'rate': 0.01, 
    #                                                                                             'noise': 0.5, 
    #                                                                                             'numCells': 1}})
    # # --- test - add minis using <netstim_inhpoisson.mod>
    # netParams.synMechParams.update({'minis_mech':{'mod': 'InhPoissonStim'}})

                
    
    NetPyNE_BBP.Prompt.headerMsg('Processing edges to create conns')

    # --- Selecting type of MOD file to be used
    if   cfg.modType == 'Det':                  modAMPANMDA = 'DetAMPANMDA';                       modGABA     = 'DetGABAAB'              # S1 Deterministic  implementation of the BBP mod files
    elif cfg.modType == 'Prob_S1':              modAMPANMDA = 'ProbAMPANMDA_EMS_S1';               modGABA     = 'ProbGABAAB_EMS_S1'      # S1 Probabilistic  implementation of the BBP mod files (used in Fernando S1 model)
    elif cfg.modType == 'Prob_original':        modAMPANMDA = 'ProbAMPANMDA_EMS_original';         modGABA     = 'ProbGABAAB_EMS_original' # original MOD from BBP model
    elif cfg.modType == 'Prob_CORENEURON':      modAMPANMDA = 'ProbAMPANMDA_EMS_CORENEURON';       modGABA     = 'ProbGABAAB_EMS_CORENEURON'         # 
    elif cfg.modType == 'Prob_S1_CORENEURON':   modAMPANMDA = 'ProbAMPANMDA_EMS_S1_CORENEURON';    modGABA     = 'ProbGABAAB_EMS_S1_CORENEURON'      # S1 Probabilistic  implementation of the BBP mod files (used in Fernando S1 model)
    else:                                       modAMPANMDA = 'ProbAMPANMDA_EMS';                  modGABA     = 'ProbGABAAB_EMS'         # Original Thalamus implementation of the BBP mod files
    print('\n\t>>\tMOD template\tAMPA: ', modAMPANMDA, '\tGABA: ', modGABA)
        
    # --- Processing edges to create conns
    edge_modif_dict = {}
    for edge_name in store_edges_dict.keys():
        # print('\t>>\tProcessing ', edge_name)
        edge_type,edge_mech,edge_source_pop,edge_target = edge_name.split('__')
        print('\n\t>>\tProcessing ', edge_target)
        print('\t\t--\t',edge_source_pop,'\t|\t',edge_type,'\t|\t',edge_mech)
        
        if (cfg.removeChemical)     and edge_mech == 'chemical':            continue
        if (cfg.removeElectrical)   and edge_mech == 'electrical_synapse':  continue

        selected_edge_properties_list=[ '@target_node','@source_node','conductance','sec','afferent_section_pos',
                                        'efferent_section_pos',
                                        'u_syn','depression_time','facilitation_time','decay_time','n_rrp_vesicles','delay',
                                        #'NMDA_ratio', # updated later
                                       ]
        
        # - Obs: (synsPerConn == 1) bacause each contact is represented as an individual edge, therefore, they are already accounted for in the dataset
        if       edge_mech       == 'electrical_synapse':           mod_file = 'Gap'
        else:
            if   edge_source_pop == 'CorticoThalamic_projections':  mod_file = modAMPANMDA
            elif edge_source_pop == 'MedialLemniscus_projections':  mod_file = modAMPANMDA
            elif edge_source_pop == 'thalamus_neurons|VPL_TC':      mod_file = modAMPANMDA
            else:                                                   mod_file = modGABA

            TableS1         = NetPyNE_BBP.StoreParameters.getTableS1()
            cell_properties = NetPyNE_BBP.LoadBBPCircuit.getNodeProperties(circuit_dict,int(edge_target))
            try:
                syn_values      = TableS1[edge_source_pop][cell_properties['mtype']]
            except:
                print('Skipping invalid pathway:\t',edge_source_pop,' -> ', cell_properties['mtype'])
                continue
            syn_values['paper_reference_values'].update({'n_rrp_vesicles':1}) # only present in the Prob MOD, not in the Det MOD

        # --- Creating a dictionary with information referenced by edge gids
        bbp_mechs_dict={}

        for edge_index in store_edges_dict[edge_name].index:
            edge_data = store_edges_dict[edge_name].loc[edge_index]
            bbp_mechs_dict.update({edge_data.name:{}})
            edge_properties = edge_data.index
            edge_vals       = edge_data.values
            for edge_property_ind,edge_property in enumerate(edge_properties):
                if edge_property in selected_edge_properties_list: bbp_mechs_dict[edge_data.name].update({edge_property:edge_vals[edge_property_ind]})
            # --- Updates the NMDA_ratio value based on data from Table S2 on the paper
            if edge_mech == 'chemical': bbp_mechs_dict[edge_data.name].update({'NMDA_ratio':syn_values['paper_reference_values']['NMDA_ratio']}) # --- set to None in GABA mechanisms

        # # --- Connection Overrides - alters the edge properties to match the parameters that were changed before model run, but are not stored in the original model
        if cfg.th_connOverrides: edge_modif_dict.update({edge_name:NetPyNE_BBP.CreateNetPyNE.modifySyn(circuit_dict, edge_name)})
        # if cfg.th_connOverrides: syn_mod = NetPyNE_BBP.CreateNetPyNE.modifySyn(circuit_dict, edge_name)

        # --- Creating NetPyNE syn mechs
        edge_dicts={}
        for mech_id in bbp_mechs_dict.keys():
            edge_dict = NetPyNE_BBP.CreateNetPyNE.modTemplate(bbp_mechs_dict[mech_id],mod_file)
            # --- Connection Overrides - alters the edge properties to match the parameters that were changed before model run, but are not stored in the original model
            if cfg.th_connOverrides: edge_dict.update(edge_modif_dict[edge_name]) # --- updates the edge dict with the modifications from BlueConfig file
            edge_dicts.update({mech_id:edge_dict})

        for mech_id in bbp_mechs_dict.keys():
            netParams.synMechParams.update({edge_source_pop+'__'+str(mech_id)+'__vecstim_mech': edge_dicts[mech_id]})
        
        # --- Creating NetPyNE conns
        for mech_id in bbp_mechs_dict.keys():
            # print('conn: \t',edge_source_pop+'__'+str(bbp_mechs_dict[mech_id]['@source_node']),' --> ',bbp_mechs_dict[mech_id]['@target_node'])
            
            netParams.connParams.update({edge_source_pop+'__'+str(mech_id)+'__vecstim_conn': {  'preConds':  {'pop': edge_source_pop+'__'+str(bbp_mechs_dict[mech_id]['@source_node'])+'__vecstim_pop'}, 
                                                                                                'postConds': {'pop': cell_properties['mtype']+'__'+str(bbp_mechs_dict[mech_id]['@target_node'])+'__pop'},
                                                                                                'probability':  1,
                                                                                                'weight':       1,
                                                                                                'synsPerConn':  1, # - Obs: (synsPerConn = 1) bacause each contact is represented as an individual edge, therefore, they are already accounted for in the dataset
                                                                                                'synMech':      edge_source_pop+'__'+str(mech_id)+'__vecstim_mech',
                                                                                                'sec':          bbp_mechs_dict[mech_id]['sec'],
                                                                                                'loc':          bbp_mechs_dict[mech_id]['afferent_section_pos'],
                                                                                                }})
            if   edge_mech == 'chemical':                                   netParams.connParams[edge_source_pop+'__'+str(mech_id)+'__vecstim_conn'].update({'delay':       bbp_mechs_dict[mech_id]['delay']})
            # elif edge_mech == 'electrical_synapse':                         netParams.connParams[edge_source_pop+'__'+str(mech_id)+'__vecstim_conn'].update({'gapJunction': True})
            if   'efferent_section_pos' in bbp_mechs_dict[mech_id].keys():  netParams.connParams[edge_source_pop+'__'+str(mech_id)+'__vecstim_conn'].update({'preLoc':      bbp_mechs_dict[mech_id]['efferent_section_pos']})
            
            # # --- test - add minis using <netstim_inhpoisson.mod>
            # netParams.connParams.update({edge_source_pop+'__'+str(mech_id)+'__minis_conn': {    'preConds':  {'pop': edge_source_pop+'__'+str(bbp_mechs_dict[mech_id]['@source_node'])+'__minis_pop'}, 
            #                                                                                     'postConds': {'pop': cell_properties['mtype']+'__'+str(bbp_mechs_dict[mech_id]['@target_node'])+'__pop'},
            #                                                                                     'probability':  1,
            #                                                                                     'weight':       bbp_mechs_dict[mech_id]['conductance'],
            #                                                                                     'synsPerConn':  1, # - Obs: (synsPerConn = 1) bacause each contact is represented as an individual edge, therefore, they are already accounted for in the dataset
            #                                                                                     'synMech':      edge_source_pop+'__'+str(mech_id)+'__vecstim_mech',
            #                                                                                     # 'synMech':      'minis_mech',
            #                                                                                     # 'delay':        bbp_mechs_dict[mech_id]['delay'],
            #                                                                                     'sec':          bbp_mechs_dict[mech_id]['sec'],
            #                                                                                     'loc':          bbp_mechs_dict[mech_id]['afferent_section_pos'],
            #                                                                                     }})
            # if   edge_mech == 'chemical':                                   netParams.connParams[edge_source_pop+'__'+str(mech_id)+'__minis_conn'].update({'delay':       bbp_mechs_dict[mech_id]['delay']})
            # elif edge_mech == 'electrical_synapse':                         netParams.connParams[edge_source_pop+'__'+str(mech_id)+'__minis_conn'].update({'gapJunction': True})
            # if   'efferent_section_pos' in bbp_mechs_dict[mech_id].keys():  netParams.connParams[edge_source_pop+'__'+str(mech_id)+'__minis_conn'].update({'preLoc':      bbp_mechs_dict[mech_id]['efferent_section_pos']})

    # --- Rescale conn weights
    for conn_name in netParams.connParams.keys(): netParams.connParams[conn_name]['weight']*=cfg.rescale_conn_weight[conn_name.split('__')[0]]

    # --- Rescale USE parameter (probability of synapse activation)
    if cfg.rescaleUSE is not None:
        for mech in netParams.synMechParams.keys():
            try:    netParams.synMechParams[mech]['Use']*=cfg.rescaleUSE
            except: continue

    # --- Identifies which thalamic populations are present in the instance of the model
    pops_dict={}
    for pop_mtype in list(set(circuit_dict['thalamus_neurons']['mtype'])):
        pops_dict.update({pop_mtype:[]})
        pops_dict[pop_mtype]=[pop_ for pop_ in netParams.popParams.keys() if (pop_mtype in pop_) and ('vecstim' not in pop_)]

    ########################################################################################################################################
    # source_pops = [edge_name.split('__')[2] for edge_name in store_edges_dict.keys()]
    # store_syn_drive = NetPyNE_BBP.Utils.calculateSynapticDrive(netParams, source_pops,cfg.duration,plotFig=True)
    ########################################################################################################################################

    # if cfg.addMinis:


    # --- Adds a current stim to thalamic populations
    if cfg.add_current_stims:
        NetPyNE_BBP.Prompt.headerMsg('Adding Current stim (IClamp)')
        for ind_th_pop,th_pop in enumerate(list(set(circuit_dict['thalamus_neurons']['mtype']))):
        # for ind_th_pop,th_pop in enumerate(['VPL_TC']):
            if len(pops_dict[th_pop])>0:
                netParams.stimSourceParams['IClamp_'+str(ind_th_pop)] = {'type': 'IClamp', 'del': cfg.current_stim_start, 'dur': cfg.duration, 'amp': cfg.current_stim_amp}
                netParams.stimTargetParams['IClamp_'+str(ind_th_pop)+'__'+th_pop] = {'source': 'IClamp_'+str(ind_th_pop), 'sec':'soma_0', 'loc': 0.5, 'conds': {'pop':pops_dict[th_pop]}}
    
    # --- Adds a noise current stim to thalamic populations
    if cfg.add_current_stims_noise:
        
        # netParams.stimSourceParams['IClamp2'] = {'type': 'IClamp', 'del': 0, 'dur': cfg.skipTime, 'amp': cfg.noise_string}
        netParams.stimSourceParams['IClamp2'] = {'type': 'IClamp', 'del': cfg.current_stim_start, 'dur': cfg.current_stim_duration, 'amp': 'normal(0,0.001)'}
        
        for cell in netParams.cellParams.keys():
            netParams.stimTargetParams['IClamp2__'+cell+'__'+sec] = {'source': 'IClamp2', 'sec':'soma_0', 'loc': 0.5, 'conds': {'cellType':cell.split('__')[1]}}
            
            # for sec in netParams.cellParams[cell]['secs'].keys():
            #     if 'axon' in sec:continue
            #     netParams.stimTargetParams['IClamp2__'+cell+'__'+sec] = {'source': 'IClamp2', 'sec':sec, 'loc': 0.5, 'conds': {'cellType':cell.split('__')[1]}}

                # print(sec)



        # NetPyNE_BBP.Prompt.headerMsg('Adding Current stim (IClamp)')
        # for ind_th_pop,th_pop in enumerate(list(set(circuit_dict['thalamus_neurons']['mtype']))):
        # # for ind_th_pop,th_pop in enumerate(['VPL_TC']):
        #     if len(pops_dict[th_pop])>0:
        #         netParams.stimSourceParams['IClamp2_'+str(ind_th_pop)] = {'type': 'IClamp', 'del': 0, 'dur': cfg.duration, 'amp': 'normal(0,0.001)'}
        #         # netParams.stimTargetParams['IClamp2_'+str(ind_th_pop)+'__'+th_pop] = {'source': 'IClamp2_'+str(ind_th_pop), 'sec':'soma_0', 'loc': 0.5, 'conds': {'pop':pops_dict[th_pop]}}
        #         netParams.stimTargetParams['IClamp2_'+str(ind_th_pop)+'__'+th_pop] = {'source': 'IClamp2_'+str(ind_th_pop), 'sec':'all', 'loc': 0.5, 'conds': {'pop':pops_dict[th_pop]}}


    #################

    '''
    # --- Retrieving information about the presynaptic inputs
    for thal_gid in thal_gids:
        # --- Generates the cell connectivity

        # define pre cells
        if cfg.th_spikes_file is not None:

            # --- Loading dictionary of spike times
            spk_dict = NetPyNE_BBP.LoadSimResults.getSpikesFromSim(filePath=cfg.th_spikes_file, showFig=False)

            rt_indexes=list(circuit_dict['mc2']['Rt_RC'].index)
            rt_spk_dict={}
            ii=0
            for spk_source in spk_dict.keys():
                if spk_source in rt_indexes:
                    print('spk source in rt',ii)
                    ii+=1
                    rt_spk_dict.update({spk_source:spk_dict[spk_source]})

            spk_dicts = {'thalamus_neurons':spk_dict}

            noise_filePaths=[('CorticoThalamic_projections',cfg.ct_virtual_noise),('MedialLemniscus_projections',cfg.ml_virtual_noise)]
            noise_dict={}
            for (key,noise_filePath) in noise_filePaths: noise_dict.update({key:NetPyNE_BBP.LoadSimResults.getVirtualInputSpikes(noise_filePath)})



        
        # retrieve spikes
        # create mechs
        # create pops
        # create conns


            NetPyNE_BBP.CreateNetPyNE.createConnsFromSpkFile(edges,thal_gid,cell_properties,stimType=cfg.th_stimType,numCells=1,mod_type=cfg.modType)

        else: 
            mechs_dict, preCells_dict, prePops_dict, conns_dict = NetPyNE_BBP.CreateNetPyNE.createConns(edges,thal_gid,cell_properties,stimType=cfg.th_stimType,numCells=1,mod_type=cfg.modType)
        
        if cfg.th_boostFibers:
            edge_name = 'internal__chemical__thalamus_neurons|Rt_RC__'+str(thal_gid)
            print('Changing parameters for ', edge_name)
            from random import sample
            boost_fibers_percent=100
            boost_fibers = int((boost_fibers_percent/100)*len(edges[edge_name]))
            list_prePops=[k for k in prePops_dict[edge_name].keys() if cfg.th_stimType in k]
            boosted_fibers = sample(list_prePops,boost_fibers)
            for fiber in boosted_fibers: prePops_dict[edge_name][fiber].update({'pulses': [{'start': 2000, 'end': 3000, 'rate': 10, 'noise':0}]})

        # --- Add paramenters into netpyne
        for edge_name in edges.keys():
            netParams.synMechParams.update( mechs_dict[edge_name])
            netParams.cellParams.update(    preCells_dict[edge_name])
            netParams.popParams.update(     prePops_dict[edge_name])
            netParams.connParams.update(    conns_dict[edge_name])
    '''

elif cfg.connType == 'testSyn':
    
    # # debug cell
    # ## Cell types
    # PYRcell = {'secs': {}}

    # PYRcell['secs']['soma'] = {'geom': {}, 'mechs': {}}
    # PYRcell['secs']['soma']['geom'] = {'diam': 18.8, 'L': 18.8, 'Ra': 123.0}
    # PYRcell['secs']['soma']['mechs']['hh'] = {'gnabar': 0.12, 'gkbar': 0.036, 'gl': 0.003, 'el': -70}

    # PYRcell['secs']['dend'] = {'geom': {}, 'topol': {}, 'mechs': {}}
    # PYRcell['secs']['dend']['geom'] = {'diam': 5.0, 'L': 150.0, 'Ra': 150.0, 'cm': 1}
    # PYRcell['secs']['dend']['topol'] = {'parentSec': 'soma', 'parentX': 1.0, 'childX': 0}
    # PYRcell['secs']['dend']['mechs']['pas'] = {'g': 0.0000357, 'e': -70}

    # netParams.cellParams['PYR'] = PYRcell
    # netParams.popParams['S'] = {'cellType': 'PYR', 'numCells': 20}

    # --- Adding cell diversity rule
    for thal_gid in thal_gids:
        cell_properties = circuit_dict['thalamus_neurons'].loc[thal_gid]
        # --- Adds each GID morphology into the cell diversity dictionary
        netParams.cellParams[cell_properties['mtype']+'__'+str(thal_gid)].update({'diversityFraction':1/len(thal_gids)})
        # --- Changes cellType to the same value, so that the different morphologies are added to the same pop can be identified by connParams
        netParams.cellParams[cell_properties['mtype']+'__'+str(thal_gid)]['conds']['cellType']=cell_properties['mtype']

    # --- Creates a single pop for all cells and defines popParams with diversity rule so that each cell has a different morphology, based on its GIDs. This aproach makes it easier to target stims and run parallel sims
    netParams.popParams[cfg.ts_targetMType+'__pop']={'cellType':'VPL_TC', 'numCells': len(thal_gids),'diversity': True}
    
    # --- Adjusting threshold to plot PSPs in the raster plot
    netParams.defaultThreshold = -73.75

    # --- Creates a single stim for all cells in the target pop
    edge_name = cfg.ts_edgeSource[0]+''
    edge_type,edge_mech,edge_source,edge_target = edge_name.split('__')

    if   cfg.modType == 'Det':                  modAMPANMDA = 'DetAMPANMDA';                       modGABA     = 'DetGABAAB'
    elif cfg.modType == 'Prob_S1':              modAMPANMDA = 'ProbAMPANMDA_EMS_S1';               modGABA     = 'ProbGABAAB_EMS_S1'
    elif cfg.modType == 'Prob_original':        modAMPANMDA = 'ProbAMPANMDA_EMS_original';         modGABA     = 'ProbGABAAB_EMS_original' # original MOD from BBP model
    elif cfg.modType == 'Prob_CORENEURON':      modAMPANMDA = 'ProbAMPANMDA_EMS_CORENEURON';       modGABA     = 'ProbGABAAB_EMS_CORENEURON'         # 
    elif cfg.modType == 'Prob_S1_CORENEURON':   modAMPANMDA = 'ProbAMPANMDA_EMS_S1_CORENEURON';    modGABA     = 'ProbGABAAB_EMS_S1_CORENEURON'
    else:                                       modAMPANMDA = 'ProbAMPANMDA_EMS';                  modGABA     = 'ProbGABAAB_EMS'

    if edge_mech == 'chemical':
        if      'Rt_RC'  in edge_source:    mod_file=modGABA
        elif    'VPL_IN' in edge_source:    mod_file=modGABA
        else:                               mod_file=modAMPANMDA
    elif edge_mech == 'electrical_synapse': mod_file='Gap'
    
    TableS1 = NetPyNE_BBP.StoreParameters.getTableS1()
    syn_values = TableS1[cfg.ts_edgeSource[2]][cfg.ts_targetMType]
    if cfg.modType == 'Prob':syn_values['paper_reference_values'].update({'n_rrp_vesicles':1}) # only present in the Prob MOD, not in the Det MOD

    edge_dict = NetPyNE_BBP.CreateNetPyNE.modTemplate(syn_values['paper_reference_values'],mod_file)

    if   cfg.ts_edgeSource[2] == 'CorticoThalamic_projections': syn_delay = 'normal(2, 0.25)'
    elif cfg.ts_edgeSource[2] == 'MedialLemniscus_projections': syn_delay = 'normal(2, 0.25)'
    else:                                                       syn_delay = 'normal(3.2, 1.0)'

    # --- NetPyNE
    # --- Creating NetPyNE syn mech
    netParams.synMechParams.update( { 'test_mech':          edge_dict})
    # --- Creates NetPyNE Vecstim Pop
    netParams.popParams.update(     { 'vecstim_pop':    {   'cellModel':    'VecStim', 'numCells': 1, 'spkTimes': cfg.ts_spkTimes}})
    # --- Creating NetPyNE conns
    netParams.connParams.update(    { 'vecstim_conn':   {   'preConds':     {'pop': 'vecstim_pop'}, 
                                                            # 'postConds':    {'pop': cfg.ts_targetMType+'__pop','mtype': 'VPL_TC'},
                                                            'postConds':    {'pop': cfg.ts_targetMType+'__pop'},
                                                            'probability':  1,
                                                            'weight':       1,
                                                            'synsPerConn':  round(cfg.ts_edgeSource[3]),
                                                            'synMech':      'test_mech',
                                                            'delay':        syn_delay, # setting delay to 1 (instead of 0) prevents MPI error: <"usable mindelay is 0 (or less than dt for fixed step method)">
                                                            # 'delay':        1, # setting delay to 1 (instead of 0) prevents MPI error: <"usable mindelay is 0 (or less than dt for fixed step method)">
                                                            'sec':          'soma_0',
                                                            'loc':          0.5,
                                                            }})