import os, sys

import h5py
import numpy as np
import pandas as pd
import json

import matplotlib.pyplot as plt

def compareTraces(cellnumber):

	f = h5py.File('/home/fernando/Documentos/ca1data20191017/sonata/edges/edges.h5', 'r') # 2.4 Gb not in github
	f1 = h5py.File('../info/data-bbp/20191017/sonata/nodes/nodes.h5', 'r')

	mtypes = ['SLM_PPA','SO_BP','SO_BS','SO_OLM','SO_Tri','SP_AA','SP_BS','SP_CCKBC','SP_Ivy','SP_PC','SP_PVBC','SR_SCA']
	syn_type_id = np.array(f['edges']['hippocampus_neurons__hippocampus_neurons__chemical']['0']['syn_type_id'])

	histog = plt.hist(syn_type_id, bins=500)

	listofsyntypes = []
	synNumber = {}

	number = 0
	for i,y in enumerate(histog[0]):
	    if y>0:
	    	number += y
	    	# print('%d [%.1f %.1f] %d %d' % (int(histog[1][1+i]),histog[1][i],histog[1][1+i],y,number))
	    	listofsyntypes.append(int(histog[1][1+i]))
	    	synNumber[int(histog[1][1+i])] = int(y)

	conductance = []
	conductance.append(np.array(f['edges']['hippocampus_neurons__hippocampus_neurons__chemical']['0']['conductance'])[np.where(syn_type_id==4)])
	histg = plt.hist(conductance, bins=200)

	syns_mtypes = {}
	for j in listofsyntypes:
	    syns_mtypes[j] = []

	for j in listofsyntypes[cellnumber:cellnumber+1]:

	    for i in list([np.where(syn_type_id==j)])[0][0]:
	    	source_id = f['edges']['hippocampus_neurons__hippocampus_neurons__chemical']['source_node_id'][i]
	    	target_id = f['edges']['hippocampus_neurons__hippocampus_neurons__chemical']['target_node_id'][i]

	    	pre_mtype = mtypes[f1['nodes']['hippocampus_neurons']['0']['mtype'][source_id]]
	    	post_mtype = mtypes[f1['nodes']['hippocampus_neurons']['0']['mtype'][target_id]]
	    	if pre_mtype + ':' + post_mtype not in syns_mtypes[j]:
	    		syns_mtypes[j].append(pre_mtype + ':' + post_mtype)
	print(syns_mtypes)		    
		    

if __name__ == '__main__':
    if len(sys.argv) == 1:
        print ("Script need a cell number between 0 and 18197")
    elif int(sys.argv[1])>=0 and int(sys.argv[1])<=18197:
        # print ("Comparing BBP and Netpyne Traces of:")
        cellnumber = int(sys.argv[1])

        print (cellnumber)

        compareTraces(cellnumber)             
            
    else:
        raise Exception('Script need a cell number between 0 and 18197')		    
