"""
batch.py 

Batch simulation for dend_thalamus project

Contributors:   joao.moreira@dowsntate.edu // salvadordura@gmail.com

OUTDATED SCRIPT

"""
from netpyne.batch import Batch
from netpyne import specs
import numpy as np


def func(simProperties):
    from cfg import cfg
    import NetPyNE_BBP
    
    thal_gid=40408

    # initial config
    initCfg = {}

    # initCfg['useTableVals']=True
    # initCfg['targetSoma']=True

    # initCfg['select_thal_gids']  = [[thal_gid]]

    edges_fileName      = cfg.NetPyNE_node_pathway_edges+'/'+str(thal_gid)+'_edges.json'
    try:
        node_pathway_edges  = NetPyNE_BBP.LoadBBPCircuit.loadNodeEdges(edges_fileName)
    except:
        print('Pathway Edges file missing - Please generate it running init.py or netParams.py')

    circuit_dict    = NetPyNE_BBP.LoadBBPCircuit.getDataFrames(cfg_file=cfg.sonataConfigFile, microcircuit_number=cfg.mc_number)
    edge_names      = NetPyNE_BBP.LoadBBPCircuit.getEdgeNames(cfg_file=cfg.sonataConfigFile)

    # Params config
    params = specs.ODict()

    params['edge_source']       = []

    # for i,edge_name in enumerate(edge_source_pops):
    #     edge_flag = edge_name.split('__')[2]
    #     for j in range(10):
    #         params['edge_source'].append((edge_name,str(edge_source_ids[i][j]),edge_flag))
    
    # for e_name in node_pathway_edges.keys():
    #     edge_sources = list(set(node_pathway_edges[e_name]['@source_node']))
    #     edge_sources.sort()
    #     print('edge_sources: ',edge_sources,'\n')
    # edge_source_pops = list(node_pathway_edges.keys())
    # edge_source_ids = [list(node_pathway_edges[k]['@source_node']) for k in list(node_pathway_edges.keys())]

    if simProperties['simType']=='CT':
        edge_name_ = 'external__chemical__CorticoThalamic_projections__'+str(thal_gid)
    elif simProperties['simType']=='ML':
        edge_name_ = 'external__chemical__MedialLemniscus_projections__'+str(thal_gid)
    elif simProperties['simType']=='RT':
        edge_name_ = 'internal__chemical__thalamus_neurons|Rt_RC__'+str(thal_gid)

    edge_source_ids = list(set(node_pathway_edges[edge_name_]['@source_node']))
    edge_source_ids.sort()

    print('edge_source_ids: ', edge_source_ids)


    edge_flag = edge_name_.split('__')[2]
    for j in range(1):
        params['edge_source'].append((edge_name_,str(edge_source_ids[j]),edge_flag))
    
    # params['edge_source']       = [
    #                                 ('MedialLemniscus_projections',1221,'ML'),
    #                                 ('thalamus_neurons|Rt_RC',33473,'RT'),
    #                                 ('CorticoThalamic_projections',83099,'CT'),
    #                                 ]
    
    groupedParams = [] 

    b = Batch(params=params, netParamsFile='netParams.py', cfgFile='cfg.py', initCfg=initCfg, groupedParams=groupedParams)

    b.batchLabel    = simProperties['simType']+'_'+simProperties['simCode']
    b.saveFolder    = simProperties['simFolder'] + '/' + simProperties['simLabel'] + '/' + simProperties['simDate'] + '/' +  b.batchLabel
    b.method        = 'grid'

    return b


# ----------------------------------------------------------------------------------------------
# Run configurations
# ----------------------------------------------------------------------------------------------
def setRunCfg(b, init_script = 'init.py'):
    b.runCfg = {'type': 'mpi_bulletin', 
                'script': init_script, 
                'skip': True} # skip sims that already exist

# ----------------------------------------------------------------------------------------------
# Main code
# ----------------------------------------------------------------------------------------------

if __name__ == '__main__':

    # simType = 'CT'
    # simType = 'ML'
    simType = 'RT'

    simProperties={
                    'simFolder':    '../data/synapses/batch_sims',
                    'simLabel':     'characterize_synapses',
                    'simDate':      '2023_03_23',
                    'simType':      simType,
                    'simCode':      'v00_tune04',
                    }
    
    b = func(simProperties)

    # ----------------------------------------------------------------------------------------------

    setRunCfg(b)
    b.run() # run batch

# ----------------------------------------------------------------------------------------------
# Command shortcuts
# ----------------------------------------------------------------------------------------------

'''
batch sim:
mpiexec -n 8 nrniv -python -mpi batch.py

init sim:
mpiexec -n 8 nrniv -python -mpi init.py

github repo commit:
sudo git commit -a -m " "; sudo git push; sudo git pull

commit_sim '<commit message>'
    using the alias setup in .bashrc file:
        commit_sim(){ git commit -a -m "$@"; sudo git push; sudo git pull}

'''