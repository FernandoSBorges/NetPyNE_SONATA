"""
Script to build the 3D representation of the thalamus network using point processes
"""
import pandas as pd
import numpy as np
import json
import math
import os

import Build_Net as BN

from netpyne import specs, sim

# base_dir = '/Users/joao'
base_dir = os.path.expanduser("~")

NetPyNE_rootFolder              = base_dir+'/Research/Models/BBP/thalamus_netpyne'
NetPyNE_network_template        = NetPyNE_rootFolder+'/conn/barreloid_network_template/network_template.json'
BBP_conn_properties             = NetPyNE_rootFolder+'/conn/calculate_BBP_conn_properties/BBP_conn_propeties.json'

netParams = specs.NetParams()   # object of class NetParams to store the network parameters
simConfig = specs.SimConfig()   # object of class SimConfig to store the simulation configuration

simConfig.duration        = 0.001*1e3

###############################################################################
# NETWORK PARAMETERS
###############################################################################

# --- Get templates for cell, pops and syn mechs
netParams.cellParams = BN.BuildNetwork.getCellTemplate()
netParams.popParams  = BN.BuildNetwork.getPopTemplate()
align_cells=False
netParams.popParams.update(BN.BuildNetwork.getMLePopTemplate(align_cells=align_cells))
netParams.synMechParams = BN.BuildNetwork.getSynMechParams()

# --- Loads the conn data calculated from the projections in the BBP somatosensory thalamus model
conn_data = BN.NetworkConversion.convertConnProperties(filePath=BBP_conn_properties)

# --- Connect using cilinder with exponential decay or uniform distribution
# prePop,postPop,sec/listOfSections/secList
connect_pops={
                'cilinder_exp':[
                                ('VPM','TRN','soma_0'),
                                ('VPM','L6A','soma_0'),
                                ('L6A','VPM','soma_0'),
                                ('L6A','TRN','soma_0'),
                                ],
                'cilinder_uni':[
                                ('TRN','VPM','soma_0'),
                                ('TRN','TRN','soma_0'),
                                ],
                }

for conn_type in connect_pops.keys():
    for (pre_pop,post_pop,secList) in connect_pops[conn_type]:
        if   conn_type=='cilinder_exp':(conn_method,conn_rule)=BN.BuildNetwork.cilinderProjection_expDecay(pre_pop,post_pop)
        elif conn_type=='cilinder_uni':(conn_method,conn_rule)=BN.BuildNetwork.cilinderProjection_uniform( pre_pop,post_pop)

        if pre_pop=='TRN':syn_mech='inh'
        else:             syn_mech='exc'
        conn_type='chem' # 2023_10_30 - Adding compatibility to 'chem' and 'elec' synapses
        # syns per conn - obs: conn_data is organized as conn_data[post_pop][pre_pop][synsPerConn/convergence][MEAN,STD]
        try:    
            syns_per_conn = round(conn_data[post_pop][pre_pop]['synsPerConn'][0])
            if (conn_data[post_pop][pre_pop]['synsPerConn'][0]>0) and int(conn_data[post_pop][pre_pop]['synsPerConn'][0])==0: syns_per_conn=1
        except: syns_per_conn = 1

        conn_dict = BN.BuildNetwork.getConnDict(pre_pop,post_pop,conn_method,conn_rule,syn_mech,syns_per_conn,conn_type,target_secs=secList)
        netParams.connParams.update(conn_dict)

# --- MLe to VPM connections
theta_pops = [pop_name for pop_name in netParams.popParams.keys() if 'MLe' in pop_name]
(mle_conn_method, mle_conn_rule)=BN.BuildNetwork.laminarProjection()


for pre_pop in theta_pops:
    mle_pop = pre_pop.split('__')[0]
    secList = 'soma_0'
    conn_type='chem' # 2023_10_30 - Adding compatibility to 'chem' and 'elec' synapses
    # syns per conn - obs: conn_data is organized as conn_data[post_pop][pre_pop][synsPerConn/convergence][MEAN,STD]
    try:    
        syns_per_conn = round(conn_data['VPM']['MLe'][0])
        if (conn_data['VPM']['MLe'][0]>0) and int(conn_data['VPM']['MLe'][0])==0: syns_per_conn=1
    except: syns_per_conn = 1

    conn_dict = BN.BuildNetwork.getConnDict(pre_pop=mle_pop,post_pop='VPM',conn_method=mle_conn_method,conn_rule=mle_conn_rule,syn_mech='exc',syns_per_conn=syns_per_conn,conn_type=conn_type,target_secs=secList)
    netParams.connParams.update(conn_dict)

############################################################################################################
# --- Plotting
simConfig.analysis['plot2Dnet'] = {'figSize':(8, 20),'saveFig': '___model_2dnet__.png'}                                                # plot 2D cell positions and connections

############################################################################################################
# --- Running simulation
from netpyne import sim
sim.createSimulateAnalyze(netParams = netParams, simConfig = simConfig)

############################################################################################################
exportNet=False
if exportNet:
    store_pops={}
    for pop_name in sim.net.pops.keys():        
        store_pops.update({pop_name:{}})
        pop_tags = sim.net.pops[pop_name].tags
        pop_gids = sim.net.pops[pop_name].cellGids
        store_pops[pop_name].update({ 'tags':pop_tags, 
                                      'cellGids':pop_gids,
                                    })
    
    store_cells={}
    for cell_gid in range(len(sim.net.cells)):  
        store_cells.update({cell_gid:{}})
        cell_tags  = sim.net.cells[cell_gid].tags
        cell_conns = sim.net.cells[cell_gid].conns
        for conn in cell_conns:
            if 'hObj' in conn.keys(): del conn['hObj']

        store_cells[cell_gid].update({  'gid':  cell_gid,
                                        'tags': cell_tags,
                                        'conns':cell_conns
                                     })

    store_network = {'pops': store_pops,'cells':store_cells}
    # with open('network_template.json', 'w') as outfile:json.dump(store_network, outfile, indent = 4)
    with open(NetPyNE_network_template, 'w') as outfile:json.dump(store_network, outfile, indent = 4)