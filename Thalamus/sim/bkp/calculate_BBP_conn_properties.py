'''
Code renamed from calculate_synsPerConn.py
'''

from cfg import cfg
import NetPyNE_BBP
import numpy as np
import collections
import numpy as np
import json

############################################################################################################################################
def merge(list1, list2):
    merged_list = [[(list1[i], list2[i])] for i in range(0, len(list1))]
    return merged_list

############################################################################################################################################

# target_pops = ['VPL_TC','Rt_RC']
target_pops = {'VPL_TC':['Rt_RC'], 'Rt_RC':['VPL_TC','Rt_RC']}

pathway_edges={}
for target_pop in target_pops.keys():
    pathway_edges.update({target_pop:{}})
    print(' --- Calculating syns per conn for population ', target_pop)
    pathway_edges[target_pop] = NetPyNE_BBP.LoadBBPCircuit.getNodeEdges(    cfg.sonataConfigFile,
                                                                {'population': 'thalamus_neurons', 'mtype': target_pop},
                                                                special_conditions={'inter_sources':target_pops[target_pop],'inter_targets':[],'select_microcircuit':None}
                                                                )

conn_pairs={} # convergence
conn_average=collections.defaultdict(int)

conn_pairs_={} # divergence
conn_average_=collections.defaultdict(int)

for target_pop in target_pops.keys():
    conn_pairs.update({target_pop:{}})
    conn_average.update({target_pop:{}})

    conn_pairs_.update({target_pop:{}})
    conn_average_.update({target_pop:{}})

    for edge_name in pathway_edges[target_pop].keys():
        ed = edge_name.split('__')[2]+'__'+edge_name.split('__')[1]
        print(ed)
        pathway_edges_sources = list(pathway_edges[target_pop][edge_name]['@source_node'])
        pathway_edges_targets = list(pathway_edges[target_pop][edge_name]['@target_node'])

        pop_conns=merge(pathway_edges_sources,pathway_edges_targets)
        count_conns = collections.defaultdict(int)

        # Using iteration
        for elem in pop_conns: count_conns[elem[0]] += 1 # Adds +1 for every time a pair of pre-post cell is identified
        # Averages the number of connections between pairs, to identify the average number of synaptic contacts (repeated connections between pre-post cells)
        conn_pairs[target_pop].update({ed:{'stats':(np.mean(list(count_conns.values())),np.std(list(count_conns.values()))),
                                        #    'vals':count_conns,
                                        #    'fibers':len(count_conns.values()), # debug this
                                           }})

        # Count the average number of connections received by a cell
        conn_repeats = collections.defaultdict(int)
        for (pre_cell,post_cell) in count_conns.keys():conn_repeats[post_cell]+=1

        conn_average[target_pop].update({ed:{'stats':(np.mean(list(conn_repeats.values())),np.std(list(conn_repeats.values()))),
                                            #  'vals':conn_repeats,
                                             }})
        
        # --- Divergence        
        # Count the average number of connections received by a cell
        conn_repeats_ = collections.defaultdict(int)
        for (pre_cell,post_cell) in count_conns.keys():conn_repeats_[pre_cell]+=1
        
        conn_average_[target_pop].update({ed:{'stats':(np.mean(list(conn_repeats_.values())),np.std(list(conn_repeats_.values()))),
                                            #  'vals':conn_repeats,
                                             }})


conn_stats=collections.defaultdict(int)
for target_pop in target_pops.keys():
    target_name='target__'+target_pop
    conn_stats.update({target_name:{}})
    
    for edge_name in pathway_edges[target_pop].keys():
        ed = edge_name.split('__')[2]+'__'+edge_name.split('__')[1]
        print(target_pop,'|',ed,'\tsyns per conn (mean, std): ',conn_pairs[target_pop][ed]['stats'],'\taverage conns per cell (mean, std): ',conn_average[target_pop][ed]['stats'])

        source_name = 'source__'+ed
        conn_stats[target_name].update({source_name:{'synsPerConn':0,'convergence':0}})
        conn_stats[target_name][source_name]['synsPerConn'] = conn_pairs[target_pop][ed]['stats']
        conn_stats[target_name][source_name]['convergence'] = conn_average[target_pop][ed]['stats']
        conn_stats[target_name][source_name]['divergence']  = conn_average_[target_pop][ed]['stats']
    
    print('\n')

update_conn_properties=False
if update_conn_properties:
    conn_stats_fileName='../conn/calculate_BBP_conn_properties/BBP_conn_propeties.json'
    json_object = json.dumps(conn_stats, indent=4)
    with open(conn_stats_fileName, "w") as outfile: outfile.write(json_object)