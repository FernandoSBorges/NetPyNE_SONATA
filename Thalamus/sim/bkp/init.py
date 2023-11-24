"""
init.py

Starting script to run NetPyNE-based thalamus model for thesis.

Usage:
    python init.py # Run simulation, optionally plot a raster

MPI usage:
    mpiexec -n 8 nrniv -python -mpi init.py

Contributors: joao.moreira@downstate.edu, salvadordura@gmail.com
"""


# snippet of code to import matplotlib and dynamically switch backend to "MacOSX" for plotting
from pydoc import source_synopsis
import sys
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
from netpyne import sim

# cfg, netParams = sim.readCmdLineArgs(simConfigDefault='cfg.py', netParamsDefault='pop_test.py')
cfg, netParams = sim.readCmdLineArgs(simConfigDefault='cfg.py', netParamsDefault='netParams.py')
# sim.create(netParams, cfg)
# sim.createSimulateAnalyze(netParams, cfg)

sim.initialize(
    simConfig = cfg, 	
    netParams = netParams)  				# create network object and set cfg and net params
sim.net.createPops()               			# instantiate network populations
sim.net.createCells()              			# instantiate network cells based on defined populations
sim.net.connectCells()            			# create connections between cells based on params
sim.net.addStims() 							# add network stimulation

try: sim.net.cells[0].secs['soma_0']['synMechs'][0]['hObj'].verboseLevel=1
except: print('skipping verbose in section soma_0')

# for i in range(10):
#     sec = 'dend_'+str(i)
#     try:
#         sim.net.cells[0].secs[sec]['synMechs'][0]['hObj'].verboseLevel=1
#         print('changing verbose for synapses in section ', sec)
#     except:
#         print('skipping section ', sec)


sim.setupRecording()              			# setup variables to record for each cell (spikes, V traces, etc)
sim.runSim()                      			# run parallel Neuron simulation  
sim.gatherData()                  			# gather spiking data and cell info from each node
sim.analyze()







# runPostAnalysis=True
# if runPostAnalysis:
#     import analyze_init
#     analyze_init.StoreSimInfo.dumpSimName(cfg)
#     if cfg.connType == 'original':
#         analysis = analyze_init.AnalyzeOriginal()
#         analysis.plotFig(filename=cfg.filename)
#     elif cfg.connType == 'testSyn':
#         analysis = analyze_init.AnalyzeTestSyn()
#         analysis.plotFig(filename=cfg.filename)



# def dumpSimName(cfg):
#     filePath=cfg.NetPyNE_data+'/_loadSimName.txt'
#     # Append-adds at last
#     file1 = open(filePath, "a")  # append mode
#     file1.write(cfg.filename+"\n")
#     file1.close()
# dumpSimName(cfg)


#############################################################################################################################
'''

runPostAnalysis=True
def getPathwayGids(cfg, netParams):
    skip_gids=0;skip_gids_=skip_gids
    pathway_gids={}
    for pop_ind, pop in enumerate(netParams.popParams.keys()):
        if cfg.th_targetMType in pop: skip_gids += netParams.popParams[pop]['numCells']
    pathway_gids.update({cfg.th_targetMType:{'1st':skip_gids_,'nth':skip_gids}})
    skip_gids_=skip_gids
        
    for pathway in cfg.th_select_pathways:
        pathway_gids.update({pathway:{'1st':skip_gids_,'nth':skip_gids}})
        for pop in netParams.popParams.keys():
            if pathway in pop: pathway_gids[pathway]['nth']+=netParams.popParams[pop]['numCells']
        skip_gids=pathway_gids[pathway]['nth']
        skip_gids_=skip_gids
    return pathway_gids

def groupSpikes(pathway_colors, pathway_gids):
    spike_id=[];spike_time=[];colors=[];markers=[]
    for pathway_ind, pathway in enumerate(pathway_gids.keys()):
        for spkid_ind,spkid in enumerate(sim.allSimData['spkid']):
            if (spkid>=pathway_gids[pathway]['1st']) and (spkid<pathway_gids[pathway]['nth']):
                spike_time.append(sim.allSimData['spkt'][spkid_ind])
                spike_id.append(sim.allSimData['spkid'][spkid_ind])
                colors.append(pathway_colors[pathway_ind])
    return spike_id,spike_time,colors

if runPostAnalysis:


    if cfg.connType == 'original':
        
        import analyze_init
        analyze_init.AnalyzeOriginal.plotFig()

        save_path = '../data/cell_inputs/init_sims/'
        simFlag   = '2023_05_10_ml_ct_rt_inputs'
        if len(cfg.select_thal_gids)==1:    
            saveData_name=save_path+simFlag+'/'+cfg.th_targetMType+'_'+str(cfg.select_thal_gids[0])+'__syn_eletrophys.json'
        else:                               
            saveData_name=save_path+simFlag+'/'+cfg.th_targetMType+'__syn_eletrophys.json'
        
        
        ####################################################################################################################################################################################
        
        # --- Convert BBP gids to NetPyNE gids
        pathway_gids = getPathwayGids(cfg, netParams)
        pathway_colors=['k','darkgreen','r','blue','orange','cyan']
        # --- Groups the spike time, id and color for plotting 
        spike_id, spike_time, colors = groupSpikes(pathway_colors, pathway_gids)
        # --- Handles for plot legend
        handles = []
        for x in pathway_colors: handles.append(plt.Line2D([], [], color=x, marker="o", linewidth=0))

        # --- Create figure
        figSize = (16,9)
        fig = plt.subplots(figsize=figSize)  # Open a new figure
        plt.suptitle('Response of cell '+ str(cfg.select_thal_gids[0]))
        
        # --- Subplot 1
        plt.subplot(2, 1, 1)
        plt.scatter(spike_time,spike_id,marker='.',s=2,color=colors)
        plt.ylabel('Cell Gids')
        plt.title('Raster plot of BBP model data')
        plt.legend(handles=handles,labels=list(pathway_gids.keys()),loc='upper right')
        plt.xlim([-100,cfg.duration])

        # --- Subplot 2
        plt.subplot(2, 1, 2)
        time_vec=sim.allSimData['t']
        cell_0 = sim.allSimData['V_soma']['cell_0']
        plt.plot(time_vec,cell_0,c='k')
        t_interval = [0,cfg.duration]
        spike_threshold=0 # mV
        cell_0_peaks=sig.find_peaks(cell_0,prominence=1)
        cell_0_peaks_t=cell_0_peaks[0]
        peaks_v=[];peaks_t=[]
        for ind_t,t in enumerate(cell_0_peaks_t):
            if cell_0[t]>spike_threshold:
                peaks_t.append(time_vec[t])
                peaks_v.append(cell_0[t])

        for t_ind,t in enumerate(peaks_t):
            plt.plot(t,peaks_v[t_ind],marker='*',c='r')

        import NetPyNE_BBP
        th_spk_dict = NetPyNE_BBP.LoadSimResults.getSpikesFromSim(filePath=cfg.th_spikes_file, cfg_file=cfg.sonataConfigFile, microcircuit_number=2, showFig=False)
        TC_cell_spkts_BBP=th_spk_dict['VPL_TC'][cfg.select_thal_gids[0]]
        index_TC_cell_spkts_BBP=[63 for i in range(len(TC_cell_spkts_BBP))]
        for t_ind,t in enumerate(TC_cell_spkts_BBP):
            plt.plot(t,index_TC_cell_spkts_BBP[t_ind],marker='*',c='b')

        plt.ylabel('Membrane voltage (mV)')
        plt.xlabel('Time (ms)')
        plt.xlim([-100,cfg.duration])
        plt.title('Voltage traces - NetPyNE Firing rate: '+str(round(len(peaks_v)*1000/(t_interval[1]-t_interval[0]),2))+' Hz ('+str(len(peaks_v))+' spikes)'+'\t|\t BBP Firing rate: '+str(round(len(index_TC_cell_spkts_BBP)*1000/(t_interval[1]-t_interval[0]),2))+' Hz ('+str(len(index_TC_cell_spkts_BBP))+' spikes)')
        plt.show()

    elif cfg.connType == 'testSyn':
        save_path = '../data/synapses/init_sims/'
        # simFlag = '08_Prob_S1_testUSE'
        simFlag = '2023_04_06_fixedSynsPerConn_200cells'
        # simFlag = '08_Prob_S1'
        # simFlag = '2023_04_03'
        Traces2 = sim.analysis.plotTraces(oneFigPer='trace', overlay=1, timeRange=[1400, sim.cfg.duration], saveData=save_path+simFlag+'/'+cfg.edge_source[2]+'_syn_eletrophys.json')
        
        ref_point=sim.cfg.spk_first-1400-1
        if ('CorticoThalamic_projections' in cfg.edge_source[2]) or ('MedialLemniscus_projections' in cfg.edge_source[2]):
            trace_color='r'
        else:
            trace_color='royalblue'
        figSize = (16,9)
        fig = plt.figure(figsize=figSize)  # Open a new figure
        # number=0
        for number in range(len(Traces2[1]['tracesData'])):
            trace_keys = list(Traces2[1]['tracesData'][number].keys())
            cell_name = [k for k in trace_keys if k.startswith('cell_')]
            print(cell_name)
            plt.plot(Traces2[1]['tracesData'][number]['t'][1:],Traces2[1]['tracesData'][number][cell_name[0]]-Traces2[1]['tracesData'][number][cell_name[0]][ref_point], label= cell_name[0], linewidth=1.0)
            # plt.plot(Traces2[1]['tracesData'][number]['t'][1:],Traces2[1]['tracesData'][number][cell_name[0]]-Traces2[1]['tracesData'][number][cell_name[0]][999]-1.0*number, label= cell_name[0], linewidth=1.0)
            # number+=1
        # plt.legend(loc='upper right', bbox_to_anchor=(1.0, 1.0))
        # plt.xlim(2400,4150)
        # plt.ylim(-4.25,1.0)
        # plt.xlabel('Time (ms)', fontsize=16)
        # plt.ylabel('EPSP (mV)', fontsize=16)
        # plt.xticks(range(2400,3900,200), range(0,1500,200), fontsize=14)
        # plt.yticks([0,1], fontsize=14);
        plt.savefig(save_path+simFlag+'/'+cfg.edge_source[2]+'_syn_eletrophys.png', facecolor = 'white' , dpi=300)
        # --- Mean trace
        import numpy as np
        alltraces = []
        gmax1 = []
        for number in range(len(Traces2[1]['tracesData'])):
            trace_keys = list(Traces2[1]['tracesData'][number].keys())
            cell_name = [k for k in trace_keys if k.startswith('cell_')]
            if np.max(Traces2[1]['tracesData'][number][cell_name[0]])-np.min(Traces2[1]['tracesData'][number][cell_name[0]]) > 0.01:
                alltraces.append(Traces2[1]['tracesData'][number][cell_name[0]])
                gmax1.append(np.max(Traces2[1]['tracesData'][number][cell_name[0]])-np.min(Traces2[1]['tracesData'][number][cell_name[0]]))
        figSize = (16,9)
        fig = plt.figure(figsize=figSize)  # Open a new figure
        for number in range(np.shape(alltraces)[0]):
            # plt.plot(alltraces[number]-alltraces[number][ref_point],color = trace_color,linewidth=0.25)
            plt.plot(alltraces[number]-alltraces[number][ref_point],color = 'lightgray',linewidth=0.25)
        plt.plot(np.mean(alltraces, axis=0)-np.mean(alltraces, axis=0)[ref_point],color = trace_color, linewidth=2.0)
        # plt.legend(loc='upper right', bbox_to_anchor=(0.95, 1.0))
        # plt.xlim(500,2500)
        # plt.ylim(-0.25,12.25)
        plt.xlabel('Time (ms)', fontsize=16)
        plt.ylabel('PSP (mV)', fontsize=16)
        # plt.xticks(range(500,3500,1000), range(0,300,100), fontsize=14);
        # plt.yticks([0,-0.5], fontsize=14); 
        plt.savefig(save_path+simFlag+'/'+cfg.edge_source[2]+'_mean_syn_eletrophys.png', facecolor = 'white' , dpi=300)
    
'''