from loadinfosfromBBP import *
from netpyne import specs, sim   

# Network parameters
netParams = specs.NetParams()  # object of class NetParams to store the network parameters

popLabel = {}
cellParamLabels = []
cellName = 'first_'
gid_list = [18097,18140,18163,18189,16950,17199,17411,17486,13513,17946,22,1]
gid_list = gid_list - np.ones_like(gid_list)

for gid in gid_list:
    MorphoName = nodesinfo['morphology'][gid] + '.swc'
    hocName = nodesinfo['model_template'][gid][4:]  
    mcName = nodesinfo['region'][gid][:3]  
    Mtype = nodesinfo['mtype'][gid]
    MEName = nodesinfo['mtype'][gid] + '_' + nodesinfo['etype'][gid]
    
    if cellName == nodesinfo['mtype'][gid] + '_' + nodesinfo['etype'][gid] + '_' + nodesinfo['region'][gid][:3]:
        cellName = nodesinfo['mtype'][gid] + '_' + nodesinfo['etype'][gid] + '_' + nodesinfo['region'][gid][:3] + '_2'
    else:  
        cellName = nodesinfo['mtype'][gid] + '_' + nodesinfo['etype'][gid] + '_' + nodesinfo['region'][gid][:3]

    cellParamLabels.append(cellName)
    popLabel[cellName] = Mtype

    print('\n\n%s \n gid = %d hoc = %s swc = %s\n' % (cellName,gid,hocName,MorphoName))

    cellRule = netParams.importCellParams(label=cellName, somaAtOrigin=False,
        conds={'cellType': cellName, 'cellModel': 'HH_full'},
        fileName='cellwrapper.py',
        cellName='loadCell',
        cellInstance = True,
        cellArgs={'hocName': hocName, 'MorphoName': MorphoName})
    netParams.renameCellParamsSec(label=cellName, oldSec='soma_0', newSec='soma')

    print(netParams.cellParams.keys())
    print(netParams.cellParams[cellName]['secs']['soma']['mechs'].keys())

    netParams.popParams[cellName] = {'cellType': cellName, 'numCells': 1, 'cellModel': 'HH_full'}

# Options
durationstim = 400.0
delaystim = 531.0
timesimulation = 1131.0
ampstim =  [-0.8, -0.6, -0.2, 0.2 ,0.4, 0.6 , 0.8, 1.0]
step_number = 4
netParams.stimSourceParams['Input'] = {'type': 'IClamp', 'del': delaystim, 'dur': durationstim, 'amp': ampstim[step_number]}

netParams.stimTargetParams['Input->all'] = {'source': 'Input', 'sec':'soma', 'loc': 0.5, 'conds': {'pop':cellParamLabels}}

## cfg  
cfg = specs.SimConfig()					            # object of class SimConfig to store simulation configuration
cfg.duration = timesimulation 						            # Duration of the simulation, in ms
cfg.dt = 0.01								                # Internal integration timestep to use
cfg.verbose = False							                # Show detailed messages 
cfg.recordTraces = {'V_soma':{'sec':'soma','loc':0.5,'var':'v'}}  # Dict with traces to record
cfg.recordStep = 0.01 
cfg.printRunTime = 0.1 # in sec			
cfg.filename = 'model_12cells'  			# Set file output name
cfg.saveJson = False
cfg.analysis['plotTraces'] = {'include': [0,1,2,3,4,5,6,7,8,9,10,11], 'timeRange': [400,1000], 'ylim': [-90,30], 'saveFig': True, 'showFig': True, 'figSize':(12,4)} # Plot recorded traces for this list of cells
cfg.analysis['plotShape'] = {'includePre': [0,1,2,3,4,5,6,7,8,9,10,11],'includePre': [0,1,2,3,4,5,6,7,8,9,10,11], 'saveFig': True, 'showFig': True, 'figSize':(12,12)}
cfg.hParams['celsius'] = 34.0

sim.createSimulateAnalyze(netParams = netParams, simConfig = cfg)


# sim.analysis.spikes.plotSpikeStats(include=cfg.allpops, saveData='../data/'+cfg.simLabel[0:9]+'/'+cfg.simLabel + '_rate.json', stats=['rate'], saveFig=False)
# sim.analysis.spikes.plotSpikeStats(include=cfg.allpops, saveData='../data/'+cfg.simLabel[0:9]+'/'+cfg.simLabel + '_CV.json', stats=['isicv'], saveFig=True)
# sim.analysis.spikes.plotSpikeStats(include=cfg.allpops, saveData='../data/'+cfg.simLabel[0:9]+'/'+cfg.simLabel + '_sync.json', stats=['sync'], saveFig=False)

# sim.analysis.plotRaster(include= ['OLM','BS','Basket','AA'], saveFig='../data/'+cfg.simLabel[0:9]+'/'+cfg.simLabel + '_rasterInh.png', labels = 'legend', popRates = True)
# sim.analysis.plot2Dnet(saveFig='../data/'+cfg.simLabel[0:9]+'/'+cfg.simLabel + '_xy.png', showConns=False)

# additional plots
#sim.analysis.plotTraces([('Pyramidal',0),('Basket',0)], saveFig=1, oneFigPer='trace', overlay=0)
# sim.analysis.plotConn(graphType='matrix', saveFig=1)
# sim.analysis.plotConn(graphType='bar', saveFig=1)
# sim.analysis.plotSpikeStats(stats = ['rate', 'isicv', 'sync', 'pairsync'], saveFig=1)
# sim.analysis.plotLFP(NFFT=256*10, noverlap=48*10, nperseg=64*10, saveFig=True)
# sim.analysis.granger(cells1=['EC'], cells2=['Pyramidal'], label1='EC', label2='Pyramidal')

# create and plot
#sim.saveData()
#sim.create()
#sim.analysis.plot2Dnet(include = ['AA', ('EC',[0,1,2]),('Pyramidal',[0,1,2]), ('CA3',[0,1,2])])
#sim.analysis.plotConn(include = ['allCells'], feature='strength', groupBy= 'pop', figSize=(9,9), showFig=True)
#sim.analysis.plotShape(includePost= ['Pyramidal','AA','Basket','BS','OLM'], showFig=True, includeAxon=True, showSyns=True)
#sim.analysis.plotRaster(saveFig=True, labels = 'legend', popRates = True)