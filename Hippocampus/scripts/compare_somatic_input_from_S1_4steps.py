import h5py
import json
import numpy as np
import os
import sys
from matplotlib import pyplot as plt

from loadinfosfromBBP import *

rootFolder = '/home/fernando/CA1_netpyne/sim/'
os.chdir(rootFolder)

savedata = 1 # Save Netpyne and BBP soma_voltage

def runneuron(cellnumber):
        
    os.chdir(rootFolder)

    gid = cellnumber
    MorphoName = nodesinfo['morphology'][gid] + '.swc'
    hocName = nodesinfo['model_template'][gid][4:]  
    
    cellName = hocName + '_' + nodesinfo['morphology'][gid]

    from cellwrapper import loadCell
    cell = loadCell(hocName, MorphoName)

    soma = cell.soma[0]

    BBPTraces = []
    BBPTracesList = []

    i=0
    for x in range(4):
        i=i+1   

        stimulus = neuron.h.IClamp(0.5, sec=soma)

        stimulus.dur = durationstim # ms
        stimulus.delay = delaystim  # ms    
        
        stimulus.amp = ampstim[x]

        recordings = {}

        recordings['time'] = neuron.h.Vector()
        recordings['soma(0.5)'] = neuron.h.Vector()

        recordings['time'].record(neuron.h._ref_t, recordStep)
        recordings['soma(0.5)'].record(cell.soma[0](0.5)._ref_v, recordStep)

        neuron.h.dt = timeStep
        neuron.h.celsius = 34
        neuron.h.v_init = -65
        # neuron.h.cvode_active(0)
        neuron.h.tstop = timesimulation # ms
        neuron.h.run();

        time = np.array(recordings['time'])
        soma_voltage = np.array(recordings['soma(0.5)'])

        BBPTraces.append(soma_voltage)
        BBPTracesList.append(list(soma_voltage))
    
    return BBPTraces

def runnetpyne(cellnumber):
  
    os.chdir(rootFolder)

    from netpyne import sim
    from netpyne import specs
    import pickle

    cfg = specs.SimConfig()     
    
    cfg.duration = timesimulation ## Duration of the sim, in ms  
    cfg.dt = timeStep
    cfg.hParams = {'celsius': 34, 'v_init': -65}  
    cfg.verbose = False
    cfg.createNEURONObj = True
    cfg.createPyStruct = True
    # cfg.cvode_active = False
    
    cfg.includeParamsLabel = False
    cfg.printPopAvgRates = True
    cfg.checkErrors = False
    
    allpops = ['L1_1','L1_2','L1_3','L1_4']

    cfg.recordCells = allpops  # which cells to record from
    cfg.recordTraces = {'V_soma': {'sec':'soma', 'loc':0.5, 'var':'v'}}  ## Dict with traces to record
    cfg.recordStim = True
    cfg.recordTime = True
    cfg.recordStep = recordStep            

    cfg.simLabel = 'Ca1_comparation'
    cfg.saveFolder = '.'
    # cfg.filename =                	## Set file output name
    cfg.savePickle = False         	## Save pkl file
    cfg.saveJson = False           	## Save json file
    cfg.saveDataInclude = ['simConfig', 'netParams'] ## 'simData' , 'simConfig', 'netParams'
    cfg.backupCfgFile = None 		##  
    cfg.gatherOnlySimData = False	##  
    cfg.saveCellSecs = False			##  
    cfg.saveCellConns = False		##  

    #------------------------------------------------------------------------------
    # Analysis and plotting 
    #------------------------------------------------------------------------------
    # ~ cfg.analysis['plotTraces'] = {'include': [('L1_1',0)], 'saveFig': True, 'showFig': False, 'oneFigPer':'trace', 'overlay':False} 		
    #------------------------------------------------------------------------------
    # Cells
    #------------------------------------------------------------------------------
    cfg.cellmod =  {'L1_1': 'HH_full'}
    cfg.cellmod =  {'L1_2': 'HH_full'}
    cfg.cellmod =  {'L1_3': 'HH_full'}
    cfg.cellmod =  {'L1_4': 'HH_full'}

    #------------------------------------------------------------------------------
    # Current inputs 
    #------------------------------------------------------------------------------
    cfg.addIClamp = 1

    cfg.IClamp1 = {'pop': 'L1_1', 'sec': 'soma', 'loc': 0.5, 'start': delaystim, 'dur': durationstim, 'amp': ampstim[0]}
    cfg.IClamp3 = {'pop': 'L1_2', 'sec': 'soma', 'loc': 0.5, 'start': delaystim, 'dur': durationstim, 'amp': ampstim[1]}
    cfg.IClamp5 = {'pop': 'L1_3', 'sec': 'soma', 'loc': 0.5, 'start': delaystim, 'dur': durationstim, 'amp': ampstim[2]}
    cfg.IClamp7 = {'pop': 'L1_4', 'sec': 'soma', 'loc': 0.5, 'start': delaystim, 'dur': durationstim, 'amp': ampstim[3]}


    netParams = specs.NetParams()   # object of class NetParams to store the network parameters

    #------------------------------------------------------------------------------
    # Cell parameters
    #------------------------------------------------------------------------------
    gid = cellnumber
    MorphoName = nodesinfo['morphology'][gid] + '.swc'
    hocName = nodesinfo['model_template'][gid][4:]      
    cellName = hocName + '_' + nodesinfo['morphology'][gid]

    cellRule = netParams.importCellParams(label=cellName, somaAtOrigin=False,
        conds={'cellType': cellName, 'cellModel': 'HH_full'},
        fileName='cellwrapper.py',
        cellName='loadCell',
        cellInstance = True,
        cellArgs={'hocName': hocName, 'MorphoName': MorphoName})
    netParams.renameCellParamsSec(label=cellName, oldSec='soma_0', newSec='soma')
    #------------------------------------------------------------------------------
    # Population parameters
    #------------------------------------------------------------------------------

    netParams.popParams['L1_1'] = {'cellType': cellName, 'cellModel': 'HH_full', 'numCells': 1} 
    netParams.popParams['L1_2'] = {'cellType': cellName, 'cellModel': 'HH_full', 'numCells': 1} 
    netParams.popParams['L1_3'] = {'cellType': cellName, 'cellModel': 'HH_full', 'numCells': 1} 
    netParams.popParams['L1_4'] = {'cellType': cellName, 'cellModel': 'HH_full', 'numCells': 1} 

    #------------------------------------------------------------------------------
    # Current inputs (IClamp)
    #------------------------------------------------------------------------------
    if cfg.addIClamp:
         for key in [k for k in dir(cfg) if k.startswith('IClamp')]:
            params = getattr(cfg, key, None)
            [pop,sec,loc,start,dur,amp] = [params[s] for s in ['pop','sec','loc','start','dur','amp']]

            #cfg.analysis['plotTraces']['include'].append((pop,0))  # record that pop

            # add stim source
            netParams.stimSourceParams[key] = {'type': 'IClamp', 'delay': start, 'dur': dur, 'amp': amp}

            # connect stim source to target
            netParams.stimTargetParams[key+'_'+pop] =  {
                'source': key, 
                'conds': {'pop': pop},
                'sec': sec, 
                'loc': loc}

    sim.createSimulateAnalyze(netParams, cfg)
    
    netpyneTraces = []
    netpyneTracesList = []
    for c in range(4):
        netpyneTraces.append(np.array(sim.simData['V_soma']['cell_'+ str(c)]))
        netpyneTracesList.append(list(sim.simData['V_soma']['cell_'+ str(c)]))        
 
    return netpyneTraces

def compareTraces(cellnumber):
    BBPTraces = runneuron(cellnumber)
    netpyneTraces = runnetpyne(cellnumber)
    # plot both traces overlayed
    fontsiz=18
    timeRange = [0, timesimulation]
    # ~ ylim = [-100, 40]
    figSize = (15,10)
    fig = plt.figure(figsize=figSize)  # Open a new figure
     
    t = np.arange(timeRange[0], timeRange[1]+recordStep, recordStep) 
     
    for c in range(0,4):
        netpyneTrace = netpyneTraces[c]
        BBPTrace = BBPTraces[c]
        plt.subplot(4, 1, c+1)
        plt.ylabel('V (mV)', fontsize=fontsiz)
        plt.plot(t[:len(netpyneTrace)], netpyneTrace, linewidth=3.5, color='red', label='I = %.1f, NetPyNE' % (ampstim[c]))
        plt.plot(t[:len(BBPTrace)], BBPTrace, linewidth=2.0, color='blue', label='I = %.1f, BBPdet' % (ampstim[c]))  # linestyle=':'
        plt.xlabel('Time (ms)', fontsize=fontsiz)
        plt.xlim(0, timesimulation)
        # ~ plt.ylim(ylim)
        plt.grid(True)
        plt.legend(loc='upper right', bbox_to_anchor=(1.15, 1.0))
    plt.ion()
    plt.tight_layout()

    MEName = nodesinfo['mtype'][gid] + '_' + nodesinfo['etype'][gid]
    plt.savefig('../info/Figures-comparation-Vsoma_neuronXNetpyne/' + MEName + '/comparison_traces_soma_voltage_4steps_%s.png' % cellName)
    print ('Figure Saved in ../info/Figures-comparation-Vsoma_neuronXNetpyne/' + MEName + '/comparison_traces_soma_voltage_4steps_%s.png' % (cellName))    
    # ~ plt.show()

if __name__ == '__main__':
    if len(sys.argv) == 1:
        print ("Script need a cell number between 0 and 18197")
    elif int(sys.argv[1])>=0 and int(sys.argv[1])<=18197:
        print ("Comparing BBP and Netpyne Traces of:")
        cellnumber = int(sys.argv[1])

        gid_list_mc2 = load_gid_list_mc2_Column()
        mtype = 'SP_BS'
        gid_list = load_gid_list(mtype)		

        i=0
        gid_mtype = []
        size_gid_list = np.size(gid_list)
        while np.size(gid_mtype)<=cellnumber and i < size_gid_list:
            gid = gid_list[i]
            if gid in gid_list_mc2:
                gid_mtype.append(gid-1) # start from 1 in original list       
            i+=1    

        if np.size(gid_mtype)>cellnumber:
            gid = gid_mtype[cellnumber]

            MorphoName = nodesinfo['morphology'][gid] + '.swc'
            hocName = nodesinfo['model_template'][gid][4:]  

            cellName = hocName + '_' + nodesinfo['morphology'][gid]

            #combolist Morpholist_hoclist[0:3993]

            print ("gidNumber = %d" % gid)
            print ("CellName = %s" % nodesinfo['mtype'][gid] + '_' + nodesinfo['etype'][gid] + '_' + nodesinfo['region'][gid][:3])
            print ("comboName = %s" % cellName)

            # durationstim = 400.0
            # delaystim = 531.0
            # timesimulation = 1131.0
            # ampstim =  [-0.8, -0.6, -0.2, 0.2 ,0.4, 0.6 , 0.8, 1.0]

            timeStep = 0.025
            recordStep = 0.025
            durationstim = 200.0
            delaystim = 200.0
            timesimulation = 600.0
            ampstim =  [0.2, 0.4 ,0.6, 0.8]

            compareTraces(cellnumber)             
            
    else:
        raise Exception('Script need a cell number between 0 and 18197')
