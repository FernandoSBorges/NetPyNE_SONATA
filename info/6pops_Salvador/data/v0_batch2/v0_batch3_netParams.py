from netpyne import sim, specs
from neuron import gui
import matplotlib.pyplot as plt
import numpy as np

netParams = specs.NetParams()


#############################################
####		POPULATION PARAMETERS		#####
#############################################
# layer_heights = (0, 170, 230, 510, 660)
# layers = ('SO', 'SP', 'SR', 'SLM')
layer = {'SO':[0.0, 0.258], 'SP':[0.258, 0.348], 'SR':[0.348, 0.773], 'SLM':[0.773, 1.0], 'CA3':[2.0,3.0], 'EC':[3.0, 4.0]}  # normalized layer boundaries

scalepopNum = 0.1
netParams.shape = 'cylinder' # cylindrical (column-like) volume

netParams.sizeX=440*scalepopNum
netParams.sizeZ=440*scalepopNum
netParams.sizeY=660

popNumber = {}
popNumber['SO_BS'] = np.ceil(3.22*scalepopNum)
popNumber['SO_BP'] = np.ceil(1.38*scalepopNum)
popNumber['SO_Tri'] = np.ceil(5.16*scalepopNum)
popNumber['SO_OLM'] = np.ceil(11.60*scalepopNum)
popNumber['SP_CCKBC'] = np.ceil(25.73*scalepopNum)
popNumber['SP_PVBC'] = np.ceil(38.70*scalepopNum)
popNumber['SP_BS'] = np.ceil(11.92*scalepopNum)
popNumber['SP_PC'] = np.ceil(2147.68*scalepopNum)
popNumber['SP_AA'] = np.ceil(10.56*scalepopNum)
popNumber['SP_Ivy'] = np.ceil(61.46*scalepopNum)
popNumber['SR_SCA'] = np.ceil(4.47*scalepopNum)
popNumber['SLM_PPA'] = np.ceil(0.44*scalepopNum)

print (popNumber)

layerVolume = {}
layerVolume['SO'] = 25849.714*scalepopNum
layerVolume['SP'] = 9123.429*scalepopNum
layerVolume['SR'] = 42576.000*scalepopNum
layerVolume['SLM'] = 22808.571*scalepopNum

#Amount of cells in the network
nPyramidal = int(popNumber['SP_PC'])
nOLM = int(popNumber['SO_OLM'])
nBS = int(popNumber['SP_CCKBC'] + popNumber['SP_PVBC'])
nB = int(popNumber['SO_BS'] + popNumber['SP_BS'])
nAA = int(popNumber['SP_AA'])

nCA3 = int(0.2*nPyramidal)
nEC = int(0.1*nPyramidal)
GAMMA = 10.0	 #  Hz
ECCA3DEL = 9.0	# msecs

# Population parameters
netParams.popParams['OLM'] = {'cellType': 
	'OLMcell', 'numCells': nOLM, 
	'cellModel': 'OLM_model', 
	'yRange':layer['SO']}
netParams.popParams['Pyramidal'] = {'cellType': 'Pyramidalcell', 
	'numCells': nPyramidal, 
	'cellModel': 'Pyramidal_model', 
	'yRange': layer['SP']}
netParams.popParams['BS'] = {'cellType': 
	'BScell', 'numCells': nBS, 
	'cellModel': 'BS_model',
	'yRange':layer['SP']}
netParams.popParams['Basket'] = {'cellType':
 	'Basketcell', 'numCells': nB, 
 	'cellModel': 'B_model', 
 	'yRange':layer['SP']}
netParams.popParams['AA'] = {'cellType': 
	'AAcell', 'numCells': nAA, 
	'cellModel': 'AA_model',
	'yRange':layer['SP']}

#############################################
netParams.popParams['EC']={'cellModel': 
	'NetStim', 'numCells': nEC, 'rate': GAMMA,
	'start': 0, 'noise': 0.2, 
	'yRange':layer['EC']}
netParams.popParams['CA3']={'cellModel': 
	'NetStim', 'numCells': nCA3, 'rate': GAMMA,
	'start': 0, 'noise': 0.2,
	'yRange':layer['CA3']}


#############################################
####		IMPORT CELL PARAMETERS  #####
#############################################
netParams.importCellParams(label='Pyramidalcell', conds={'cellType': 'Pyramidalcell', 'cellModel': 'Pyramidal_model'}, \
fileName='cells/pyramidal_cell_14Vb.hoc', cellName='PyramidalCell', importSynMechs=False)

netParams.importCellParams(label='OLMcell', conds={'cellType': 'OLMcell', 'cellModel': 'OLM_model'}, \
fileName='cells/olm_cell2.hoc', cellName='OLMCell', importSynMechs=False)

netParams.importCellParams(label='BScell', conds={'cellType': 'BScell', 'cellModel': 'BS_model'}, \
fileName='cells/bistratified_cell13S.hoc', cellName='BistratifiedCell', importSynMechs=False)

netParams.importCellParams(label='Basketcell', conds={'cellType': 'Basketcell', 'cellModel': 'B_model'}, \
fileName='cells/basket_cell17S.hoc', cellName='BasketCell', importSynMechs=False)

netParams.importCellParams(label='AAcell', conds={'cellType': 'AAcell', 'cellModel': 'AA_model'}, \
fileName='cells/axoaxonic_cell17S.hoc', cellName='AACell', importSynMechs=False)


##Setting thresholds
cells=['Pyramidalcell','OLMcell','BScell','Basketcell','AAcell']
for i in cells:
	for sec in netParams.cellParams[i].secs:
 		netParams.cellParams[i].secs[sec].threshold = -10.0

#############################################
####		NETWORK CONNECTIONS	#####
#############################################

weights={'Pyramidalcell2Pyramidalcell': 0.001, 'Pyramidalcell2AAcell':0.0005, 'Pyramidalcell2Basketcell':0.0005, 'Pyramidalcell2BScell':0.0005,'Pyramidalcell2OLMcell': 0.00005, \
'AAcell2Pyramidalcell': 0.04,\
'Basketcell2Pyramidalcell': 0.02, 'Basketcell2Basketcell': 0.001, 'Basketcell2BScell': 0.02,\
'BScell2Pyramidalcell': 0.002, 'BScell2Pyramidal_GABABasketcell': 0.0004, 'BScell2Basketcell': 0.01, \
'OLMcell2Pyramidalcell': 0.04, 'OLMcell2Pyramidal_GABABasketcell': 0.0004,'OLMcell2Basketcell': 0.01, }

delays={'Pyramidalcell2Pyramidalcell': 1., 'Pyramidalcell2AAcell':1., 'Pyramidalcell2Basketcell':1., 'Pyramidalcell2BScell':1.,'Pyramidalcell2OLMcell': 1., \
'AAcell2Pyramidalcell': 1., \
'Basketcell2Pyramidalcell': 1., 'Basketcell2Basketcell': 1., 'Basketcell2BScell': 1., \
'BScell2Pyramidalcell': 1., 'BScell2Pyramidal_GABABasketcell': 1., 'BScell2Basketcell': 1., \
'OLMcell2Pyramidalcell': 1., 'OLMcell2Pyramidal_GABABasketcell': 1.,'OLMcell2Basketcell': 1. }

# Cue (CA3) excitation
CHWGT = 0.0015	#// cue weight
CLWGT = 0.0005	#// unlearnt weight (usually 0)
CNWGT = 0.0005	#// excitatory weights (NMDA)
CDEL = 1.	#// cue delay

#EC excitation
ECWGT = 0.0	# EC weight to PCs
#ECWGT = 0.001	# EC weight to PCs
ECDEL = 1.	# EC delay
EIWGT = 0.00015	# excitatory weights to INs
EIDEL = 1.	# delay (msecs)


#############################################
#### DESCRIPTION OF SYNAPTIC MECHANISMS	#####
#############################################

netParams.synMechParams['GABAA']={'mod':'MyExp2Syn', 'tau1':1.0, 'tau2':8.0, 'e':-75.0}
netParams.synMechParams['GABAB']={'mod':'MyExp2Syn', 'tau1':35.0, 'tau2':100.0, 'e':-75.0}
netParams.synMechParams['AMPA']={'mod':'MyExp2Syn', 'tau1':0.5, 'tau2':3.0, 'e':0.0}
netParams.synMechParams['NMDA']={'mod':'NMDA', 'tcon': 2.3, 'tcoff': 100.0, 'enmda': 0.0, 'gNMDAmax': 1.0, 'tauD': 800.0, 'tauF': 800.0, 'util': 0.3}

netParams.synMechParams['OLM_GABAA']={'mod':'Exp2Syn', 'tau1':1.0, 'tau2':8.0, 'e':-75.0}
netParams.synMechParams['OLM_GABAB']={'mod':'Exp2Syn', 'tau1':35.0, 'tau2':100.0, 'e':-75.0}
netParams.synMechParams['OLM_AMPA']={'mod':'Exp2Syn', 'tau1':0.5, 'tau2':3.0, 'e':0.0}

#######################
##presyn = Pyramidal CHECKED
#######################

postsynList=['Pyramidal','AA','Basket','BS','OLM']
postsynDict={'Pyramidal':['radTprox'], 'AA': ['oriT1','oriT2'], 'Basket':['oriT1','oriT2'], 'BS':['oriT1','oriT2'], 'OLM':['dend1','dend2']}

for i in range(len(postsynList)):
	k='Pyramidalcell2'+postsynList[i]+'cell'
	netParams.connParams['Pyramidal->'+postsynList[i]] = {
		'preConds': {'pop': 'Pyramidal'},
		'postConds': {'pop': postsynList[i]},
		'sec': postsynDict[postsynList[i]],
		'synsPerConn':len(postsynDict[postsynList[i]]),
		'synMech': 'AMPA',
		'weight': weights[k],
		'delay': delays[k]
		#'threshold': -10.0
		}
	if postsynList[i]=='Pyramidal':
		netParams.connParams['Pyramidal->Pyramidal']['convergence'] = 100. # PC_PC = 1 percent in full scale
	if postsynList[i]=='OLM':
		netParams.connParams['Pyramidal->OLM']['synMech'] = 'OLM_AMPA'

#######################
##presyn == AA CHECKED
#######################

netParams.connParams['AA->Pyramidal'] = {
        'preConds': {'pop': 'AA'},
        'postConds': {'pop': 'Pyramidal'},
        'sec': 'axon',
		'loc': 0.1,
        'synMech': 'GABAA',
        'weight': weights['AAcell2Pyramidalcell'],
        'delay': delays['AAcell2Pyramidalcell']
		}

#######################
##presyn == B CHECKED
#######################

postsynList=['Pyramidal','Basket','BS']		##B->AA not connected

for i in range(len(postsynList)):
	k='Basketcell2'+postsynList[i]+'cell'
	netParams.connParams['B->'+postsynList[i]] = {
			'preConds': {'pop': 'Basket'},
			'postConds': {'pop': postsynList[i]},
			'sec': 'soma',
			'synMech': 'GABAA',   #GABA-A
			'weight': weights[k],
			'delay': delays[k]
			}
	if postsynList[i]=='BS': netParams.connParams['B->BS']['loc'] = 0.6

#######################
##presyn == BS  CHECKED
#######################

##BS->AA & BS->BS not connected

netParams.connParams['BS->B'] = {
	'preConds': {'pop': 'BS'},
	'postConds': {'pop': 'Basket'},
	'sec': 'soma',
	'synMech': 'GABAA',
	'loc':0.6,
	'weight': weights['BScell2Basketcell'],
	'delay': delays['BScell2Basketcell']
	}

netParams.connParams['BS->Pyramidal'] = {
		'preConds': {'pop': 'BS'},
		'postConds': {'pop': 'Pyramidal'},
		'sec': 'radTmed',
		'synsPerConn':7,
		'loc':[[0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2],[0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2]],
		'synMech': ['GABAA','GABAB'],
		'weight': [weights['BScell2Pyramidalcell'], weights['BScell2Pyramidal_GABABasketcell']],
		'delay': [delays['BScell2Pyramidalcell'],delays['BScell2Pyramidal_GABABasketcell']]
		}

#######################
##presyn == OLM  CHECKED
#######################

netParams.connParams['OLM->Pyramidal'] = {
		'preConds': {'pop': 'OLM'},
		'postConds': {'pop': 'Pyramidal'},
		'sec': ['lm_thick1','lm_thick2'],
		'synMech': ['GABAA','GABAB'],  #GABA-A,GABA-B
		'weight': [weights['OLMcell2Pyramidalcell'], weights['OLMcell2Pyramidal_GABABasketcell']],
		'delay': [delays['OLMcell2Pyramidalcell'],delays['OLMcell2Pyramidal_GABABasketcell']],
		'synsPerConn':21
		}


#############################################
####		STIMULATION - INPUTS		#####
#############################################
#####################################
#####EC 
#####################################

netParams.connParams['EC->Pyramidal'] = {
		'preConds': {'pop': 'EC'},
		'postConds': {'pop': 'Pyramidal'},
		'convergence': 100.0,
		'sec': ['lm_thick1','lm_thick2'],
		'synMech': 'AMPA',
		'loc':0.5,
		'weight': ECWGT,
		'delay': ECDEL,
		'synsPerConn':2
		}
netParams.connParams['EC->IN'] = {
		'preConds': {'pop': 'EC'},
		'postConds': {'pop': ['Basket','AA']},
		'sec': ['lmM1','lmM2'],
		'convergence': 100.0,
		'synMech': 'AMPA',
		'weight': EIWGT,
		'delay': EIDEL,
		'synsPerConn':2
		}

#############################
####CA3
############################
postsynList=['AA','Basket','BS']
for i in range(len(postsynList)):
	netParams.connParams['CA3->'+postsynList[i]] = {
		'preConds': {'pop': 'CA3'},
		'postConds': {'pop': postsynList[i]},
		'convergence': 100.0,
		'sec': postsynDict[postsynList[i]],
		'synsPerConn':len(postsynDict[postsynList[i]]),
		'synMech': 'AMPA',
		'weight': EIWGT,
		'delay': EIDEL,
		'loc':0.5
		}

netParams.connParams['CA3_highW->Pyramidal'] = {
		'preConds': {'pop': 'CA3'},
		'postConds': {'pop': 'Pyramidal'},
		'convergence': 100.0,
		'sec': 'radTmed',
		'synMech': 'AMPA',
		'loc':0.5,
		'weight': CHWGT,
		'delay': CDEL
		}
netParams.connParams['CA3_NMDA->Pyramidal'] = {
		'preConds': {'pop': 'CA3'},
		'postConds': {'pop': 'Pyramidal'},
		'sec': 'radTmed',
		'convergence': 100.0,
		# 'probability':0.01,
		'synMech': 'NMDA',
		'loc':0.5,
		'weight': CNWGT,
		'delay': CDEL
		}