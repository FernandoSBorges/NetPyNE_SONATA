import sys
import json
import matplotlib.patches as mpatches
from scipy import signal as sig
import numpy as np

from    matplotlib  import  pyplot  as plt
print("Matplotlib backend (default): %s" %plt.get_backend())
modules = []
for module in sys.modules:
    if module.startswith('matplotlib'):
        modules.append(module)
for module in modules:
    sys.modules.pop(module)
import matplotlib
matplotlib.use("MacOSX")
from    matplotlib  import  pyplot  as plt
print("Matplotlib backend (dynamic): %s" %plt.get_backend())

class Prompt():
    def headerMsg(msg):
        print(  '\n --------------------------------------------------------------------------------------------------------------------------------------\n',
                '\t\t',msg,'\t\t\n',
                '--------------------------------------------------------------------------------------------------------------------------------------')

class StoreSimInfo():
    def dumpSimName(cfg):
        filePath=cfg.saveFolder+'/_loadSimName.txt'
        # Append-adds at last
        file1 = open(filePath, "a")  # append mode
        file1.write(cfg.filename+"\n")
        file1.close()

class LoadSim():
    def loadSimObj(filename):
        with open(filename, 'r') as fileObj: simFile=fileObj.read()
        sim = json.loads(simFile)
        return sim

class AnalyzeOriginal():
    def plotFig(self,filename=None):
        # --- Load Sim object
        sim = LoadSim.loadSimObj(filename)
        # --- Adds a time shift in the reference spike times if skipTime exists (pre-run sim to stabilize variables)
        try:    skipTime = sim['simConfig']['skipTime']
        except: skipTime = 0
        # --- Convert BBP gids to NetPyNE gids
        pathway_gids = AnalyzeOriginal.getPathwayGids(sim)
        # pathway_colors = ['k','darkgreen','r','blue','orange','cyan']
        pathway_colors_dict = { 
                                'thalamus_neurons|VPL_TC':      'k',
                                'VPL_TC':                       'k',
                                'MedialLemniscus_projections':  'darkgreen',
                                'CorticoThalamic_projections':  'r',
                                'thalamus_neurons|Rt_RC':       'blue',
                                'Rt_RC':                        'blue',
                                'thalamus_neurons|VPL_IN':      'orange',
                                'VPL_IN':                       'orange',
                                'other':                        'cyan'
                                }
        
        # --- Groups the spike time, id and color for plotting 
        spike_id, spike_time, colors = AnalyzeOriginal.groupSpikes(pathway_colors_dict, pathway_gids, sim)
        # --- Handles for plot legend
        handles = []; labels = []
        for pathway in pathway_colors_dict.keys(): 
            if pathway in pathway_gids.keys():
                if 'thalamus_neurons|' in pathway: pathway.replace('thalamus_neurons|','')
                handles.append(plt.Line2D([], [], color=pathway_colors_dict[pathway], marker="o", linewidth=0))
                labels.append(pathway)

        # --- Create figure
        figSize = (9,9)
        fig = plt.subplots(figsize=figSize)  # Open a new figure
        plt.suptitle('Response of cell '+ str(sim['simConfig']['select_thal_gids'][0]))
        
        # --- Subplot 1
        plt.subplot(3, 1, 1)
        plt.scatter(spike_time,spike_id,marker='.',s=2,color=colors)
        plt.ylabel('Cell Gids')
        plt.title('Raster plot of BBP model data')
        plt.legend(handles=handles,labels=labels,loc='upper right')
        # plt.legend(handles=handles,labels=list(pathway_colors_dict.keys()),loc='upper right')
        plt.xlim([skipTime-100,sim['simConfig']['duration']])

        # --- Subplot 2
        plt.subplot(3, 1, 2)
        time_vec=sim['simData']['t']
        cell_0 = sim['simData']['V_soma']['cell_0']
        min_cell_0=min(cell_0)
        min_array= [min_cell_0 for i in range(len(cell_0))]
        
        # plt.plot(time_vec,cell_0,c='k')
        
        t_interval = [0,sim['simConfig']['duration']]
        spike_threshold=sim['net']['params']['defaultThreshold']
        # spike_threshold=0 # mV
        cell_0_peaks=sig.find_peaks(cell_0,prominence=1)
        cell_0_peaks_t=cell_0_peaks[0]
        peaks_v=[];peaks_t=[]
        for ind_t,t in enumerate(cell_0_peaks_t):
            if cell_0[t]>spike_threshold:
                peaks_t.append(time_vec[t])
                peaks_v.append(cell_0[t])

        # --- Array to hold spike times at a fixed y-axis location
        peaks_v_=[53 for i in range(len(peaks_v))]
        # --- Plots spike lines
        for t_ind,t in enumerate(peaks_t): plt.plot([t,t],[min_cell_0,peaks_v_[t_ind]],linestyle='--',linewidth=0.5,c='darkorange',zorder=1)
        # --- Plots spike markers
        # for t_ind,t in enumerate(peaks_t): plt.plot(t,peaks_v[t_ind],marker='*',markersize=2,c='darkorange',zorder=3)
        for t_ind,t in enumerate(peaks_t): plt.plot(t,peaks_v_[t_ind],marker='*',markersize=2,c='darkorange',zorder=3)

        import NetPyNE_BBP
        th_spk_dict = NetPyNE_BBP.LoadSimResults.getSpikesFromSim(filePath=sim['simConfig']['th_spikes_file'], cfg_file=sim['simConfig']['sonataConfigFile'], microcircuit_number=2, showFig=False)
        TC_cell_spkts_BBP=th_spk_dict[sim['simConfig']['target']][sim['simConfig']['select_thal_gids'][0]]
        
        TC_cell_spkts_BBP = [t+skipTime for t in TC_cell_spkts_BBP]
        
        peaks_v_BBP=[63 for i in range(len(TC_cell_spkts_BBP))]
        # --- Plots spike lines
        for t_ind,t in enumerate(TC_cell_spkts_BBP): plt.plot([t,t],[min_cell_0,peaks_v_BBP[t_ind]],linestyle='--',linewidth=0.5,c='slateblue',zorder=1)
        # --- Plots spike markers
        for t_ind,t in enumerate(TC_cell_spkts_BBP): plt.plot(t,peaks_v_BBP[t_ind],marker='*',markersize=2,c='slateblue',zorder=3)

        # --- Plot reference voltage values
        th_soma_voltage_dict = NetPyNE_BBP.LoadSimResults.loadVoltageReport(filePath=sim['simConfig']['Fig_4_1_traces'],cfg_file=sim['simConfig']['sonataConfigFile'],gids=sim['simConfig']['select_thal_gids'],dt = 0.1, timeWindow = [0,4000], microcircuit_number=sim['simConfig']['select_microcircuit'], showFig=False)
        # th_soma_voltage_dict = NetPyNE_BBP.LoadSimResults.loadVoltageReport(filePath=sim['simConfig']['Fig_4_1_traces']+'.h5',cfg_file=sim['simConfig']['sonataConfigFile'],gids=sim['simConfig']['select_thal_gids'],dt = 0.1, timeWindow = [0,4000], microcircuit_number=sim['simConfig']['select_microcircuit'], showFig=False)
        time_trace=th_soma_voltage_dict['t']
        time_trace=[t+skipTime for t in time_trace]
        voltage_trace=th_soma_voltage_dict['voltage'][sim['simConfig']['select_thal_gids'][0]]

        plt.fill_between(time_vec,min_array,cell_0,color='w',zorder=2)
        plt.plot(time_vec,cell_0,c='k',zorder=2)

        plt.plot(time_trace,voltage_trace,'k',linewidth=1,alpha=0.2)

        plt.ylabel('Membrane voltage (mV)')
        plt.xlabel('Time (ms)')
        plt.xlim([skipTime-100,sim['simConfig']['duration']])
        plt.title('Voltage traces - NetPyNE Firing rate: '+str(round(len(peaks_v)*1000/(t_interval[1]-t_interval[0]),2))+' Hz ('+str(len(peaks_v))+' spikes)'+'\t|\t BBP Firing rate: '+str(round(len(peaks_v_BBP)*1000/(t_interval[1]-t_interval[0]),2))+' Hz ('+str(len(peaks_v_BBP))+' spikes)')

        # --- Subplot 3
        plt.subplot(3, 1, 3)
        try:
            source_pops = set(list([pop.split('__')[0] for pop in sim['net']['params']['synMechParams'].keys()]))
            store_syn_drive = NetPyNE_BBP.Utils.calculateSynapticDrive(sim['net']['params'], source_pops,sim['simConfig']['duration'],plotFig=True)
            plt.xlim([skipTime-100,sim['simConfig']['duration']])
        except:print()
        ######################################################################################################################################################
        plotCurr=False
        for recordedTrace in sim['simData'].keys():
            if recordedTrace.startswith('i__'):
                plotCurr=True
                break
        if plotCurr:
            current_colors_dict = { 
                                    'SK_E2__ik'         : 'magenta',
                                    'TC_HH__ina'        : 'brown',
                                    'TC_HH__ik'         : 'teal',
                                    'TC_Nap_Et2__ina'   : 'yellow',
                                    'TC_iA__ik'         : 'green',
                                    'TC_iL__ica'        : 'r',
                                    'TC_iT_Des98__ica'  : 'darkblue',
                                    'TC_ih_Bud97__ih'   : 'orange',
                                    'other'             : 'cyan'
                                }
            # --- General fig config
            separateFigs=False
            figSize = (9,9)
            figAlpha = 0.5
            defaultLabel = 'other'
            defaultColor = 'cyan'

            current_handles=[]
            for current_color in current_colors_dict.keys(): current_handles.append(plt.Line2D([], [], color=current_colors_dict[current_color], marker="o", linewidth=0, alpha=figAlpha, markeredgecolor='k'))
                
            def plotSingleTraces(sim, time_vec, current_colors_dict):
                figSize = (9,9)
                for recordedTrace in sim['simData'].keys():
                    if recordedTrace.startswith('i__'):
                        c = defaultColor
                        for current_color in current_colors_dict.keys():
                            if current_color in recordedTrace: c=current_colors_dict[current_color]
                        plt.figure(figsize=figSize)
                        plt.plot(time_vec,sim['simData'][recordedTrace]['cell_0'],c=c)
                        plt.title(recordedTrace)

            def plotAllTraces(sim, time_vec, current_colors_dict, current_handles):
                figSize = (9,9)
                fig2 = plt.subplots(figsize=figSize)  # Open a new figure
                plt.suptitle('Currents of cell '+ str(sim['simConfig']['select_thal_gids'][0]))
                for recordedTrace in sim['simData'].keys():
                    if recordedTrace.startswith('i__'):
                        c = defaultColor
                        for current_color in current_colors_dict.keys():
                            if current_color in recordedTrace: c=current_colors_dict[current_color]
                        plt.plot(time_vec,sim['simData'][recordedTrace]['cell_0'],c=c)
                plt.legend(handles=current_handles,labels=list(current_colors_dict.keys()),loc='upper right')

            def storePlotInfo(sim, current_colors_dict):
                store_currents=[];store_labels=[];store_colors=[]
                for recordedTrace in sim['simData'].keys():
                    if recordedTrace.startswith('i__'):
                        l = defaultLabel; c = defaultColor
                        for current_color in current_colors_dict.keys():
                            if current_color in recordedTrace: 
                                l = current_color
                                c = current_colors_dict[current_color]
                        store_currents.append(sim['simData'][recordedTrace]['cell_0'])
                        store_labels.append(l)
                        store_colors.append(c)
                return store_currents, store_labels, store_colors
            
            def sortCurrents(store_currents, store_labels, store_colors):
                sum_abs_current=[]
                for current_i,current in enumerate(store_currents):
                    current_np  = np.array(current)
                    abs_current = np.abs(current_np)
                    sum_abs_current.append((sum(abs_current),current_i))
                sum_abs_current.sort(reverse=True)
                sorted_store_currents=[];sorted_store_labels=[];sorted_store_colors=[]
                for (sum_abs_curr,sorting_index) in sum_abs_current:
                    sorted_store_currents.append(store_currents[sorting_index])
                    sorted_store_labels.append(store_labels[sorting_index])
                    sorted_store_colors.append(store_colors[sorting_index])

                del store_currents, store_labels, store_colors
                    
                store_currents  = sorted_store_currents
                store_labels    = sorted_store_labels
                store_colors    = sorted_store_colors
                return store_currents, store_colors
            #########################################################################################################################

            # --- Plots each trace separately - using the mapped colors
            singleTraces=False
            if singleTraces:    plotSingleTraces(sim, time_vec, current_colors_dict)
            # --- Plots all traces together - using the mapped colors
            allTraces=False
            if allTraces:       plotAllTraces(sim, time_vec, current_colors_dict, current_handles)
            # --- Generates a list of traces and label/color info for plotting
            store_currents, store_labels, store_colors  = storePlotInfo(sim, current_colors_dict)
            store_currents, store_colors                = sortCurrents(store_currents, store_labels, store_colors)
            
            store_currents_np = np.array(store_currents)
            trace_0 = [0 for i in range(len(time_vec))]

            if separateFigs: plt.figure(figsize=figSize)
            else:
                fig = plt.subplots(figsize=figSize)  # Open a new figure
                plt.subplot(2, 1, 1)
                
            for ind,current in enumerate(store_currents_np): plt.fill_between(time_vec,trace_0,current,color=store_colors[ind],alpha=figAlpha)
            plt.legend(handles=current_handles,labels=list(current_colors_dict.keys()),loc='upper right')
            plt.title('Fill plot - current traces')

            abs_store_currents_np = np.abs(store_currents_np)
            ratio_currents = abs_store_currents_np / abs_store_currents_np.sum(axis=0)

            # --- Creates an array of currents that adds the values of the previous plotted current to it, to make the fill_between plot
            store_ratio_currents=[]
            for i, sublist in enumerate(ratio_currents): store_ratio_currents.append(sublist+np.sum(ratio_currents[0:i], axis=0))
            
            # --- Fill plot of currents ratio at each timepoint 
            if separateFigs: plt.figure(figsize=figSize)
            else:            plt.subplot(2, 1, 2)
            for ind,trace in enumerate(store_ratio_currents):
                if ind ==0: plt.fill_between(time_vec,trace_0,store_ratio_currents[ind],color=store_colors[ind],alpha=figAlpha)
                else:       plt.fill_between(time_vec,store_ratio_currents[ind-1],store_ratio_currents[ind],color=store_colors[ind],alpha=figAlpha)
            if separateFigs: plt.legend(handles=current_handles,labels=list(current_colors_dict.keys()),loc='upper right')
            plt.title('Fill plot - currents ratio')



    ################################################################################################################################
    def getPathwayGids(sim):
        skip_gids=0;skip_gids_=skip_gids
        pathway_gids={}
        # --- Counts the number of cells in the Target populations to allocate GIDs
        for pop_ind, pop in enumerate(sim['net']['params']['popParams'].keys()):
            if (sim['simConfig']['target'] in pop) and ('__pop' in pop): skip_gids += sim['net']['params']['popParams'][pop]['numCells']
        pathway_gids.update({sim['simConfig']['target']:{'1st':skip_gids_,'nth':skip_gids}})
        skip_gids_=skip_gids
            
        # --- Counts the number of cells in the Source populations to allocate GIDs
        for pathway in sim['simConfig']['th_select_pathways']:
            pathway_gids.update({pathway:{'1st':skip_gids_,'nth':skip_gids}})
            for pop in sim['net']['params']['popParams'].keys():
                if pathway in pop: pathway_gids[pathway]['nth']+=sim['net']['params']['popParams'][pop]['numCells']
            skip_gids=pathway_gids[pathway]['nth']
            skip_gids_=skip_gids
        return pathway_gids
    
    def groupSpikes(pathway_colors_dict, pathway_gids, sim):
        spike_id=[];spike_time=[];colors=[];markers=[]
        for pathway_ind, pathway in enumerate(pathway_gids.keys()):
            print(pathway)
            for spkid_ind,spkid in enumerate(sim['simData']['spkid']):
                if (spkid>=pathway_gids[pathway]['1st']) and (spkid<pathway_gids[pathway]['nth']):
                    spike_time.append(sim['simData']['spkt'][spkid_ind])
                    spike_id.append(sim['simData']['spkid'][spkid_ind])
                    colors.append(pathway_colors_dict[pathway])
        return spike_id,spike_time,colors

##############################################################################################################################################################################

class AnalyzeTestSyn():
    def plotFig(self,filename=None):
        import numpy as np
        # --- Load Sim object
        sim = LoadSim.loadSimObj(filename)
        
        time_vec=sim['simData']['t']
        cell_voltages={}
        for cell in sim['simData']['V_soma'].keys():
            cell_voltages.update({cell:sim['simData']['V_soma'][cell]})

        ref_point=sim['simConfig']['ts_spkFirst']-1400-1
        if ('CorticoThalamic_projections' in sim['simConfig']['ts_edgeSource'][2]) or ('MedialLemniscus_projections' in sim['simConfig']['ts_edgeSource'][2]):  
                trace_color='r'
        else:   trace_color='royalblue'

        figSize = (9,9)
        fig = plt.figure(figsize=figSize)  # Open a new figure
        plot_cell_voltages=[]
        for cell in cell_voltages.keys():
            cell_voltages[cell] = [v-cell_voltages[cell][ref_point] for v in cell_voltages[cell]]
            plt.plot(time_vec,cell_voltages[cell], label=cell, color='lightgray', linewidth=0.25)
            # --- Remove "zero" traces
            if np.max(cell_voltages[cell])-np.min(cell_voltages[cell]) > 0.01:
                plot_cell_voltages.append(cell_voltages[cell])

        plt.plot(np.mean(plot_cell_voltages, axis=0)-np.mean(plot_cell_voltages, axis=0)[ref_point],color = trace_color, linewidth=2.0)
        # plt.plot(np.mean(list(cell_voltages.values()), axis=0)-np.mean(list(cell_voltages.values()), axis=0)[ref_point],color = trace_color, linewidth=2.0)
        plt.xlabel('Time (ms)', fontsize=16)
        plt.ylabel('PSP (mV)', fontsize=16)
        plt.xlim([sim['simConfig']['ts_spkFirst']-500,sim['simConfig']['duration']])
        plt.savefig(sim['simConfig']['saveFolder']+'/'+sim['simConfig']['filename']+sim['simConfig']['ts_edgeSource'][2]+'_mean_syn_eletrophys.png', facecolor = 'white' , dpi=300)

##############################################################################################################################################################################

if __name__ == "__main__":
    
    dataFolder='../data/init_sims/'
    
    import tkinter as tk
    from tkinter import filedialog
    root = tk.Tk()
    root.withdraw()
    filenames = filedialog.askopenfilenames(initialdir=dataFolder,typevariable='.json',defaultextension='.json')
    root.update()

    for filename in filenames:
        if   '.json' not in filename: continue
        if   'single_cell_inputs'    in filename: foo = AnalyzeOriginal()
        elif 'single_synapse_inputs' in filename: foo = AnalyzeTestSyn()
        else: continue
        foo.plotFig(filename=filename)

    plt.show()
    root.update()