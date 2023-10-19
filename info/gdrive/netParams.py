from netpyne import sim, specs
# from neuron import gui
import matplotlib.pyplot as plt
import numpy as np

from loadinfosfromBBP import *

netParams = specs.NetParams()

try:
    from __main__ import cfg  # import SimConfig object with params from parent module
except:
    from cfg import cfg
#------------------------------------------------------------------------------
# POPULATION PARAMETERS	
#------------------------------------------------------------------------------
# layer_heights = (0, 170, 230, 510, 660)
# layers = ('SO', 'SP', 'SR', 'SLM')
layer = {'SO':[0.0, 0.258], 'SP':[0.258, 0.348], 'SR':[0.348, 0.773], 'SLM':[0.773, 1.0], 'CA3':[2.0,3.0], 'EC':[3.0, 4.0]}  # normalized layer boundaries

# full import
scalepopNum = 0.002

# mc2 import
netParams.shape = 'cylinder' # cylindrical (column-like) volume
area = (3*np.sqrt(3)*240**2)/2 # six triag l = 240, h =  240*np.sqrt(3)/2
diam = 2*np.sqrt(area/np.pi) #equiv circle
netParams.sizeX=diam
netParams.sizeZ=diam
netParams.sizeY=660

#------------------------------------------------------------------------------
Ca1_cellNumber = 0
for mtype in Mtypelist:    
    popNumber[mtype] = int(np.ceil(popNumber[mtype]*scalepopNum))
    Ca1_cellNumber = Ca1_cellNumber + popNumber[mtype]

print('scalepopNum =',scalepopNum,'  Ca1_cellNumber =',Ca1_cellNumber)
print (popNumber)

#------------------------------------------------------------------------------
# to debug
popNumber['SP_PC'] = int(1)
Ca1_cellNumber = 0
for mtype in Mtypelist:    
    Ca1_cellNumber = Ca1_cellNumber + popNumber[mtype]

print('scalepopNum =',scalepopNum,'  Ca1_cellNumber =',Ca1_cellNumber)
print (popNumber)

#------------------------------------------------------------------------------
# IMPORT CELL PARAMETERS
#------------------------------------------------------------------------------
popLabel = {}
cellParamLabels = []
cellNumber = {}
# gid_list = [18097, 18109, 18140, 18149, 18163, 18177,18189, 18191, 16950, 16963, 17199, 17202, 
# 						17411, 17416, 17486, 17497, 13513, 14311, 17946, 17958, 22, 25, 1, 3]
# gid_list = gid_list - np.ones_like(gid_list) # start from 1 in original list

gid_list_mc2 = load_gid_list_mc2_Column()

for mtype in Mtypelist:

	gid_list = load_gid_list(mtype)		

	gid_mtype = []
	for gid in gid_list:
		if gid in gid_list_mc2 and np.size(gid_mtype)<popNumber[mtype]:
			gid_mtype.append(gid-1) # start from 1 in original list

	for gid in gid_mtype:
		MorphoName = nodesinfo['morphology'][gid] + '.swc'
		hocName = nodesinfo['model_template'][gid][4:]  
		cellName = nodesinfo['mtype'][gid] + '_' + nodesinfo['etype'][gid] + '_' + nodesinfo['region'][gid][:3] + '_0'

		if cellName in cellParamLabels:
			cellNumber[cellName[:-2]] = cellNumber[cellName[:-2]] + 1
			cellName = cellName[:-2] + '_' + str(cellNumber[cellName[:-2]]-1) # [cell_0,cell_1] cellNumber = 2
		else:
			cellNumber[cellName[:-2]] = 1

		cellParamLabels.append(cellName)
		popLabel[cellName] = nodesinfo['mtype'][gid]

		print('%s  gid = %d hoc = %s swc = %s' % (cellName,gid,hocName,MorphoName))

		cellRule = netParams.importCellParams(label=cellName, somaAtOrigin=False,
			conds={'cellType': cellName, 'cellModel': 'HH_full'},
			fileName='cellwrapper.py',
			cellName='loadCell',
			cellInstance = True,
			cellArgs={'hocName': hocName, 'MorphoName': MorphoName})
		netParams.renameCellParamsSec(label=cellName, oldSec='soma_0', newSec='soma')

# 'SLM_PPA': # No SLM in mc2, import from mc3 (mc5 or mc6 are avaliable too)
gid = 0  
MorphoName = nodesinfo['morphology'][gid] + '.swc'
hocName = nodesinfo['model_template'][gid][4:]  
cellName = nodesinfo['mtype'][gid] + '_' + nodesinfo['etype'][gid] + '_' + nodesinfo['region'][gid][:3] + '_0'
cellNumber[cellName[:-2]] = 1

cellParamLabels.append(cellName)
popLabel[cellName] = nodesinfo['mtype'][gid]

print('%s  gid = %d hoc = %s swc = %s' % (cellName,gid,hocName,MorphoName))

cellRule = netParams.importCellParams(label=cellName, somaAtOrigin=False,
	conds={'cellType': cellName, 'cellModel': 'HH_full'},
	fileName='cellwrapper.py',
	cellName='loadCell',
	cellInstance = True,
	cellArgs={'hocName': hocName, 'MorphoName': MorphoName})
netParams.renameCellParamsSec(label=cellName, oldSec='soma_0', newSec='soma')

# 'SO_BP': # No SO_BP in mc2, import from mc3 (mcX are avaliable too)
gid = 18193  
MorphoName = nodesinfo['morphology'][gid] + '.swc'
hocName = nodesinfo['model_template'][gid][4:]  
cellName = nodesinfo['mtype'][gid] + '_' + nodesinfo['etype'][gid] + '_' + nodesinfo['region'][gid][:3] + '_0'
cellNumber[cellName[:-2]] = 1

cellParamLabels.append(cellName)
popLabel[cellName] = nodesinfo['mtype'][gid]

print('%s  gid = %d hoc = %s swc = %s' % (cellName,gid,hocName,MorphoName))

cellRule = netParams.importCellParams(label=cellName, somaAtOrigin=False,
	conds={'cellType': cellName, 'cellModel': 'HH_full'},
	fileName='cellwrapper.py',
	cellName='loadCell',
	cellInstance = True,
	cellArgs={'hocName': hocName, 'MorphoName': MorphoName})
netParams.renameCellParamsSec(label=cellName, oldSec='soma_0', newSec='soma')

print(netParams.cellParams[cellName]['secs']['soma']['mechs'].keys())

#------------------------------------------------------------------------------
#Amount of cells in the network
#------------------------------------------------------------------------------
for cellName in cellParamLabels:
	celllayer = cellName.split('_')[0]
	netParams.popParams[cellName] = {'cellType': cellName, 'numCells': 1, 'ynormRange': layer[celllayer], 'cellModel': 'HH_full'}

#------------------------------------------------------------------------------
#	NETWORK CONNECTIONS	
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# DESCRIPTION OF SYNAPTIC MECHANISMS	
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Stimulus
#------------------------------------------------------------------------------
durationstim = 400.0
delaystim = 531.0
timesimulation = 1131.0
ampstim =  [-0.8, -0.6, -0.2, 0.2 ,0.4, 0.6 , 0.8, 1.0]
step_number = 4
netParams.stimSourceParams['Input'] = {'type': 'IClamp', 'del': delaystim, 'dur': durationstim, 'amp': ampstim[step_number]}
netParams.stimTargetParams['Input->all'] = {'source': 'Input', 'sec':'soma', 'loc': 0.5, 'conds': {'pop':cellParamLabels}}
