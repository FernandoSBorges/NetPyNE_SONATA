import matplotlib; matplotlib.use('Agg')  # to avoid graphics error in servers
from netpyne import sim
import pyspike

cfg, netParams = sim.readCmdLineArgs()
sim.initialize(
    simConfig = cfg, 	
    netParams = netParams)  				# create network object and set cfg and net params
sim.net.createPops()               			# instantiate network populations
sim.net.createCells()              			# instantiate network cells based on defined populations
sim.net.connectCells()            			# create connections between cells based on params
sim.net.addStims() 							# add network stimulation
sim.setupRecording()              			# setup variables to record for each cell (spikes, V traces, etc)
sim.runSim()                      			# run parallel Neuron simulation  
sim.gatherData()                  			# gather spiking data and cell info from each node
sim.saveData()                    			# save params, cell info and sim output to file (pickle,mat,txt,etc)#
sim.analysis.plotData()         			# plot spike raster etc

# sim.analysis.spikes.plotSpikeStats(include=cfg.allpops, saveData='../data/'+cfg.simLabel[0:9]+'/'+cfg.simLabel + '_rate.json', stats=['rate'], saveFig=False)
# sim.analysis.spikes.plotSpikeStats(include=cfg.allpops, saveData='../data/'+cfg.simLabel[0:9]+'/'+cfg.simLabel + '_CV.json', stats=['isicv'], saveFig=True)
# sim.analysis.spikes.plotSpikeStats(include=cfg.allpops, saveData='../data/'+cfg.simLabel[0:9]+'/'+cfg.simLabel + '_sync.json', stats=['sync'], saveFig=False)

sim.analysis.plotRaster(include= ['OLM','BS','Basket','AA'], saveFig='../data/'+cfg.simLabel[0:9]+'/'+cfg.simLabel + '_rasterInh.png', labels = 'legend', popRates = True)
sim.analysis.plot2Dnet(saveFig='../data/'+cfg.simLabel[0:9]+'/'+cfg.simLabel + '_xy.png', showConns=False)

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