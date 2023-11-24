"""
init.py

Starting script to run NetPyNE-based thalamus model for thesis.

Usage:
    python barreloid_init.py # Run simulation, optionally plot a raster

MPI usage:
    mpiexec -n 8 nrniv -python -mpi barreloid_init.py

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

############################################################################################################
# --- Running simulation
cfg, netParams = sim.readCmdLineArgs(simConfigDefault='barreloid_cfg.py', netParamsDefault='barreloid_netParams.py')
createOnly=False
if createOnly:
    print('\n\n\n ---- Running CREATE ONLY mode to inspect for bugs during network creation --- \n\n\n')
    sim.create(netParams = netParams, simConfig = cfg)
else:
    # sim.create(netParams = netParams, simConfig = cfg)
    # sim.createSimulate(netParams = netParams, simConfig = cfg)
    # sim.createSimulateAnalyze(netParams = netParams, simConfig = cfg)
    # sim.create(netParams = netParams, simConfig = cfg)

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
    sim.analyze()

    extra_analysis = False
    if extra_analysis:
        rasterData = sim.analysis.prepareRaster()
        spk_color=[]
        mle_gids = []
        for pop in sim.net.pops.keys():
            if 'MLe' in pop: mle_gids+=sim.net.pops[pop].cellGids
        for spk_gid in rasterData['spkGids']:
            if   spk_gid in sim.net.pops['VPM__pop'].cellGids:  spk_color.append('b')
            elif spk_gid in sim.net.pops['TRN__pop'].cellGids:  spk_color.append('r')
            elif spk_gid in sim.net.pops['TRN__pop'].cellGids:  spk_color.append('g')
            elif spk_gid in mle_gids:                           spk_color.append('m')
            else:                                               spk_color.append('cyan')

        gids_y=[]
        for gid in rasterData['spkGids']:
            gids_y.append(sim.net.cells[int(gid)].tags['y'])

        # - ordered by gid
        plt.figure(figsize=(25,20))
        plt.scatter(rasterData['spkTimes'],rasterData['spkGids'],c = spk_color,marker='.')
        plt.savefig('figs_analysis/000001_raster_byGID.png')

        # - ordered by y
        plt.figure(figsize=(25,20))
        plt.scatter(rasterData['spkTimes'],gids_y,c = spk_color,marker='.')
        plt.gca().invert_yaxis()
        plt.savefig('figs_analysis/000001_raster_byY.png')

        sim.analysis.plotRaster(rasterData=rasterData, dpi=1000, orderBy='y', orderInverse=True, labels = False, figSize=(25, 20), saveFig= 'figs_analysis/000001_raster.png')

        '''
        sim.analysis.plotConn(includePre=['VPM__pop', 'TRN__pop', 'L6A__pop'],includePost=['VPM__pop', 'TRN__pop', 'L6A__pop'],figSize=(10, 10),saveFig='figs_analysis/fullConn_topological_2023_09_21___connMatrix_uninverted2.png')

        sim.analysis.plot2Dnet(figSize=(15, 20), saveFig='figs_analysis/MLe_to_VPM_2dnet_uninverted2.png')
        '''
        # --- Selects the source of spike IDs, because they have different names in the sim object and in rasterData
        skipPlot=False
        from collections import Counter
        try:
            spkids = rasterData['spkGids']
        except:
            try: spkids = sim.simData['spkid']
            except:
                skipPlot=True
                print('No spkids')

        if not skipPlot:
            # --- Creates bins based on the y range of the cells
            import numpy as np
            n_bins = 8
            bins = np.arange(min(sim.net.pops['VPM__pop'].tags['yRange']),max(sim.net.pops['VPM__pop'].tags['yRange'])+1,(max(sim.net.pops['VPM__pop'].tags['yRange'])-min(sim.net.pops['VPM__pop'].tags['yRange']))/n_bins)

            # --- Creates a dictionary to store cell gids based on position
            store_cell_pos = {bins[i]:0 for i in range(len(bins)-1)}

            # --- Sorts and converts spike ids to int
            spkids.sort()
            spkids_int = [int(spkid) for spkid in spkids]
            # --- Counts the number of occurences of each spike id
            spkcounts = Counter(spkids_int)

            # --- Counts the number of spikes for each cell based on its position
            for cell in sim.net.cells:
                if 'VPM' in cell.tags['pop']:
                    cell_y = cell.tags['y']
                    for i in range(len(bins)-1):
                        if cell_y > bins[i] and cell_y<=bins[i+1]:
                            store_cell_pos[bins[i]]+=spkcounts[cell.gid]

            plt.figure()
            plt.bar([i*(360/n_bins) for i in range(len(store_cell_pos.keys()))],store_cell_pos.values())
            plt.savefig(sim.cfg.saveFigPath+'/'+sim.cfg.filename+'_'+sim.cfg.sim_tag+'_spkCount'+'.png',dpi=500)

    inspect_VPM_conns=False
    if inspect_VPM_conns:
        inspect_conns_dict={}
        for cell_ind, cell in enumerate(sim.net.cells):
            conns=0
            for conn in cell.conns:
                if 'MLe' in conn['synMech']:
                    conns+=1
            if 'VPM' in cell.tags['pop']:
                print(cell, conns)
                inspect_conns_dict.update({cell_ind:conns})

        import numpy as np
        np.mean(list(inspect_conns_dict.values()))

    # sim.cfg.saveFigPath


    # if not skipPlot:
    #     spkids.sort()
    #     spkcounts = Counter(spkids)

    #     VPM_spks_ = [[],[],[],[],[],[],[],[]]
    #     for i in range(318):    
    #         div = int(318/8)+1
    #         print(i,'\t',div,'\t',int(i/div))
    #         VPM_spks_[int(i/div)].append(spkcounts[i])
    #     VPM_spks = [sum(spks) for spks in VPM_spks_]
    #     plt.figure()
    #     plt.bar([i*90 for i in range(len(VPM_spks))],VPM_spks)
    #     plt.savefig('figs_analysis/VPM_spks_uninverted2.png',dpi=500)



    '''
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
    sim.analyze()
    '''

############################################################################################################
runPostAnalysis=False
if runPostAnalysis:
    # --- Verify conns
    verifyConns=False
    if verifyConns:
        for cell_ind in range(len(sim.net.cells)):
            print('cell ', cell_ind)
            secList_proximal    = sim.net.cells[cell_ind].secLists['inputs__proximal']
            secList_intermediate= sim.net.cells[cell_ind].secLists['inputs__intermediate']
            secList_distal      = sim.net.cells[cell_ind].secLists['inputs__distal']
            if len(sim.net.cells[cell_ind].conns)==0:
                print('No conns')
                continue
            for conn_ind in range(len(sim.net.cells[cell_ind].conns)):
                conn_sec = sim.net.cells[cell_ind].conns[conn_ind]['sec']
                if   conn_sec in secList_proximal:      print('0 - proximal conn')
                elif conn_sec in secList_intermediate:  print('1 - intermediate conn')
                elif conn_sec in secList_distal:        print('2 - distal conn')

    ############################################################################################################
    # === Evaluate connectivity
    print('\t>>\tRunning network evaluation - Divergence and Convergence values')
    from EvaluateModel import Evaluate as Eval

    pop_gids  = Eval.getPopGids(sim)
    pop_names = list(pop_gids.keys())

    pop_convergence  = Eval.getPopConvergence(sim)
    pop_divergence   = Eval.getPopDivergence(sim)

    cell_convergence = Eval.getCellConvergence(sim)
    cell_divergence  = Eval.getCellDivergence(sim)

    # pop_gids,pop_names,pop_convergence,pop_divergence,cell_convergence,cell_divergence = Eval.runAll(sim)

    # --- Pathway ratios for cell convergence/divergence
    convergence_ratio = Eval.getConvergenceRatio(sim)
    divergence_ratio  = Eval.getDivergenceRatio(sim)

    print('\n------------------------------------------------------------------------------------------------------------------------------\n')

    print('\t>>\t--- Convergence Data ---')
    print('\t>>\tConvergence Absolute [(number of conns * syns per conn) for each cell]:  ')
    [print('\t\t',key,cell_convergence[key]) for key in cell_convergence.keys()]
    print('\t>>\tRatio - Convergence Absolute [(number of conns * syns per conn) for each cell]:  ')
    [print('\t\t',key,convergence_ratio[key]) for key in convergence_ratio.keys()]

    print('\n\t\t - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n')
    print('\t>>\tPathway ratios:  ')
    try: print('\t>>\t (L6A->VPM)/(L6A->MLe) ratio: ', cell_convergence['VPM']['L6A']/cell_convergence['VPM']['MLe'], '\t-- target - BBP: ', 581.4493258129902/122.97733554449974)
    except: 0
    try: print('\t>>\t (TRN->VPM)/(MLe->VPM) ratio: ', cell_convergence['VPM']['TRN']/cell_convergence['VPM']['MLe'], '\t-- target - BBP: ', 454.1682034044653/122.97733554449974)
    except: 0

    try: print('\t>>\t (L6A->VPM)/(L6A->TRN) ratio: ', cell_convergence['VPM']['L6A']/cell_convergence['TRN']['L6A'], '\t-- target - BBP: ', 581.4493258129902/255.23635294982)
    except: 0
    try: print('\t>>\t (L6A->TRN)/(VPM->TRN) ratio: ', cell_convergence['TRN']['L6A']/cell_convergence['TRN']['VPM'], '\t-- target - BBP: ', 255.23635294982/161.43545983636517)
    except: 0

    print('\n\t\t - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n')

    print('\t>>\t--- Divergence Data ---')
    print('\t>>\tDivergence Absolute [(number of conns * syns per conn) for each cell]:  ')
    [print('\t\t',key,cell_divergence[key]) for key in cell_divergence.keys()]
    print('\t>>\tRatio - Divergence Absolute [(number of conns * syns per conn) for each cell]:  ')
    [print('\t\t',key,divergence_ratio[key]) for key in divergence_ratio.keys()]

    print('\n------------------------------------------------------------------------------------------------------------------------------\n')

    # --- BBP model stats
    print('\t>>\t--- BBP Convergence Data ---')
    print('\t>>\tBBP Convergence Absolue:  ')
    [print('\t\t',key2,'-to->',key,': ',cfg.conn_data[key]['chem'][key2]['synsPerConn'][0]*cfg.conn_data[key]['chem'][key2]['convergence'][0]) for key in cfg.conn_data.keys() for key2 in cfg.conn_data[key]['chem'].keys()]

    print('\n\t\t - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n')
    print('\t>>\tPathway ratios BBP:  ')
    print('\t>>\t (L6A->VPM)/(L6A->MLe) ratio: ', 581.4493258129902/122.97733554449974)
    print('\t>>\t (TRN->VPM)/(MLe->VPM) ratio: ', 454.1682034044653/122.97733554449974)
    print('\t>>\t (L6A->VPM)/(L6A->TRN) ratio: ', 581.4493258129902/255.23635294982)
    print('\t>>\t (L6A->TRN)/(VPM->TRN) ratio: ', 255.23635294982/161.43545983636517)

    print('\n------------------------------------------------------------------------------------------------------------------------------\n')

    # --- Ratio of synapses into VB thalamus from (Van Horn and Sherman, 2007)
    driver_RL   = 28
    modul_RS    = 146
    trn_F       = 48
    sum_all = driver_RL+modul_RS+trn_F
    print('\t>>\t--- Convergence Data (Van Horn and Sherman, 2007) ---')
    print('\t>>\tConvergence Proportional:  ')
    print('\t\t MLe-to->VPM: ', driver_RL)
    print('\t\t L6A-to->VPM: ', modul_RS)
    print('\t\t TRN-to->VPM: ', trn_F)

    print('\t>>\tConvergence Ratio:  ')
    print('\t\t MLe-to->VPM: ', driver_RL/sum_all)
    print('\t\t L6A-to->VPM: ', modul_RS/sum_all)
    print('\t\t TRN-to->VPM: ', trn_F/sum_all)

    print('\n\t\t - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n')
    print('\t>>\tPathway ratios:  ')
    print('\t>>\t (L6A->VPM)/(MLe->VPM) ratio: ', modul_RS/driver_RL)
    print('\t>>\t (TRN->VPM)/(MLe->VPM) ratio: ', trn_F/driver_RL)
    # print('\t>>\tBBP Convergence ratio:  ')
    # [print('\t\t',key,conn_data[key]) for key in conn_data.keys()]

    print('\n------------------------------------------------------------------------------------------------------------------------------\n')
    # print('\n\t\t - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n')

    '''
    (çavdar, 2011) - VB terminals and synapses - Fig 8
                    driver - RL             modulator - RS          trn - F
    terminals       4.91                    86.81818181818181       6.092307692307692
    synapses        12.248803827751196      79.6969696969697        7.569230769230771       (USE THIS VALUE)
    '''
    driver_RL   = 12.248803827751196
    modul_RS    = 79.6969696969697
    trn_F       = 7.569230769230771
    sum_all = driver_RL+modul_RS+trn_F

    print('\t>>\t--- Convergence Data (çavdar, 2007) ---')
    print('\t>>\tConvergence Proportional:  ')
    print('\t\t MLe-to->VPM: ', driver_RL)
    print('\t\t L6A-to->VPM: ', modul_RS)
    print('\t\t TRN-to->VPM: ', trn_F)

    print('\t>>\tConvergence Ratio:  ')
    print('\t\t MLe-to->VPM: ', driver_RL/sum_all)
    print('\t\t L6A-to->VPM: ', modul_RS/sum_all)
    print('\t\t TRN-to->VPM: ', trn_F/sum_all)

    print('\n\t\t - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n')
    print('\t>>\tPathway ratios:  ')
    print('\t>>\t (L6A->VPM)/(MLe->VPM) ratio: ', modul_RS/driver_RL)
    print('\t>>\t (TRN->VPM)/(MLe->VPM) ratio: ', trn_F/driver_RL)

    print('\n------------------------------------------------------------------------------------------------------------------------------\n')


    # --- Bar Plot with conn values
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 36})
    plt.figure( figsize=(20,20))

    pathway_ratios={
                    '(toVPM)_L6A/MLe':  {'BBP':581.4493258129902/122.97733554449974,    'model':cell_convergence['VPM']['L6A']/cell_convergence['VPM']['MLe']},
                    '(toVPM)_TRN/MLe':  {'BBP':454.1682034044653/122.97733554449974,    'model':cell_convergence['VPM']['TRN']/cell_convergence['VPM']['MLe']},
                    '(fromL6A)_VPM/TRN':{'BBP':581.4493258129902/255.23635294982,       'model':cell_convergence['VPM']['L6A']/cell_convergence['TRN']['L6A']},
                    '(toTRN)_L6A/VPM':  {'BBP':255.23635294982/161.43545983636517,      'model':cell_convergence['TRN']['L6A']/cell_convergence['TRN']['VPM']},
                    }

    for pathway_ratio_ind,pathway_ratio in enumerate(pathway_ratios.keys()):
        for source in pathway_ratios[pathway_ratio].keys():
            width = 0.25
            if source == 'model':
                c   = 'green'
                loc =  pathway_ratio_ind+(width/2)
            elif source == 'BBP':
                c   = 'purple'
                loc =  pathway_ratio_ind-(width/2)

            plt.bar(loc,pathway_ratios[pathway_ratio][source],width=width,color=c)

    plt.xticks([0,1,2,3],['L6A/MLe ->VPM','TRN/MLe ->VPM','L6A-> VPM/TRN','L6A/VPM ->TRN'],rotation=10)
    plt.xlabel('Pathway')
    plt.ylabel('Conn Ratio')
    plt.legend(['Thal_NB model + literature','NetPyNE model'])

    plt.savefig('test_bar_plot.png',dpi=500)

    '''
    ############################################################################################################

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

    plotConns=False
    selectGids=[419,619,800,850,1071]
    if plotConns:
        for cell_gid in selectGids: 
            fig = plt.figure()
            ax = fig.add_subplot(projection='3d')
            # --- Plot Soma center - SWC files are loaded with cells centered at (0,0,0)
            ax.scatter(0,0,0, marker='o',color='cyan',linewidths=0, alpha=1)
            
            # --- Plot 3d locations of cell morphology sections
            x3d=[];     y3d=[];     z3d=[]
            for sec in sim.net.cells[cell_gid].secs.keys():
                if 'pt3d' not in sim.net.cells[cell_gid].secs[sec]['geom'].keys():continue
                if len(sim.net.cells[cell_gid].secs[sec]['geom']['pt3d'])>0 and type(sim.net.cells[cell_gid].secs[sec]['geom']['pt3d'])==list:
                    for pt3d in sim.net.cells[cell_gid].secs[sec]['geom']['pt3d']:
                        x3d.append(pt3d[0]);    y3d.append(pt3d[1]);    z3d.append(pt3d[2])
            # # --- sec/pt3d points
            ax.scatter(x3d, y3d, z3d, marker='o',color='k',s=0.25,linewidths=0, alpha=0.1)
            
            
            # syn_keys = list(sim.net.params.synMechParams.keys())
            # for syn_key in sim.net.params.synMechParams.keys():
            
            conn_dict = {}
            colors = ['b','r','g','lightblue','orange','magenta','limegreen','grey']
            for ind,conn_key in enumerate(['conn|TRN|TRN|chem','conn|TRN|VPM|chem','conn|VPM|TRN|chem','conn|MLe|VPM|chem','conn|L6A|TRN|chem','conn|L6A|VPM|chem','conn|VPM|L6A|chem']):
                # try:    conn_key = 'conn|'+syn_key.split('|')[1]+'|'+syn_key.split('|')[2]
                # except: continue
                conn_dict.update({conn_key:{'c':colors[ind],'syn':[]}})

            for conn_ind in range(len(sim.net.cells[cell_gid].conns)):
                
                postSec = sim.net.cells[cell_gid].conns[conn_ind]['sec']
                postLoc = sim.net.cells[cell_gid].conns[conn_ind]['loc']
                preSyn  = sim.net.cells[cell_gid].conns[conn_ind]['label']
                
                
                p0=sim.net.cells[cell_gid].secs[postSec]['geom']['pt3d'][0]
                p1=sim.net.cells[cell_gid].secs[postSec]['geom']['pt3d'][-1]
                p_x = (p1[0]-p0[0])*postLoc
                p_y = (p1[1]-p0[1])*postLoc
                p_z = (p1[2]-p0[2])*postLoc
                
                conn_dict[preSyn]['syn'].append((p_x,p_y,p_z))

            for syn_ind,syn_key in enumerate(conn_dict.keys()):
                try:
                    x3d = [x[0] for x in conn_dict[syn_key]['syn']]
                    y3d = [x[1] for x in conn_dict[syn_key]['syn']]
                    z3d = [x[2] for x in conn_dict[syn_key]['syn']]
                except:continue

                ax.scatter(x3d, y3d, z3d, marker='o',color=conn_dict[syn_key]['c'],s=1,linewidths=0, alpha=1)

            plt.savefig('./figs/cell_projections_'+str(cell_gid)+'.png',dpi=1000)

    plt.show()
    '''
    