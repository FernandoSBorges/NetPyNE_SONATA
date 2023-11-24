'''
Class to load and convert the BBP thalamus model into NetPyNE

Contributors: joao.moreira@downstate.edu, salvadordura@gmail.com
'''
##################################################################################################################################################################################################################

# --- Import packages
import numpy as np
import pandas as pd
import h5py
import json

##################################################################################################################################################################################################################

class Global():
    figsFolder='../figs/'

##################################################################################################################################################################################################################

class Prompt():
    def headerMsg(msg):
        print(  '\n -------------------------------------------------------------------------------------------------\n',
                '\t\t',msg,'\t\t\n',
                '-------------------------------------------------------------------------------------------------')

##################################################################################################################################################################################################################

class Conversion():
    def convertCellMorphology(inputFolder_h5,outputFolder_swc,inputFolder_asc=None):
        import os
        morphologies_h5=list(os.listdir(inputFolder_h5))
        morphologies_swc=list(os.listdir(outputFolder_swc))
        from morph_tool import converter as cv
        success=0;fail=0;already_existing_files=0
        missing_files=[];still_missing_files=[]
        for ind,morphology_h5 in enumerate(morphologies_h5):
            # --- Check which files already exist in SWC format
            if already_existing_files%1000==0:print('Already existing files: ',already_existing_files)
            if morphology_h5[:-3]+'.swc' in morphologies_swc:
                already_existing_files+=1
                continue
            else:
                print('Missing: ', morphology_h5)
                missing_files.append(morphology_h5)

            # --- Converts non-existing files into SWC and save
            if ind%1000==0:print(len(morphologies_h5)-ind,' remaining')
            morphology_swc = morphology_h5[:-3]+'.swc'
            try:
                original_file=inputFolder_h5+'/'+morphology_h5
                converted_file = outputFolder_swc+'/'+morphology_swc
                converted_morphology_swc = cv.convert(original_file,converted_file)
                success+=1
            except:
                # --- Attempts to convert from ASC if SWC conversion fails
                if inputFolder_asc is not None:
                    try:
                        print('Trying ASC conversion')
                        retry_original_file_asc = inputFolder_asc+'/'+morphology_h5[:-3]+'.asc'
                        retry_converted_file = outputFolder_swc+'/'+morphology_swc
                        retry_converted_morphology_swc = cv.convert(retry_original_file_asc,retry_converted_file)
                        success+=1
                        print('Converted: ',retry_converted_file)
                    # --- Conversion failed - Cell morphology will be missing
                    except:
                        fail+=1
                        still_missing_files.append(retry_converted_file)
                        print('Failed to convert ASC: ',retry_converted_file)
                        pass
        print('successes: ', success,'\t|\tfails: ',fail)
        print('\nOriginally missing files:\n',missing_files)
        print('\nStill missing files:\n',still_missing_files)

    def getAfferentSectionIds(morphologyFullPath):
        # --- Method to map the sections in the full cell morphology to individual IDs, used to connect the circuit edges to a specific section
        '''
        # Read more at: https://github.com/AllenInstitute/sonata/blob/master/docs/SONATA_DEVELOPER_GUIDE.md
        For the purpose of mapping information to morphologies, each section has a designated section ID, which is a integer in [0, #sections - 1]. The soma is always section 0. The rest of the sections are first grouped by section type in this order: 
        1 = axon, 2 = basal and 3 = apical. 
        All sections in a group are given IDs incrementally, starting at the last ID from the previous section plus 1. The order in which sections are assigned an ID is the order in which their first segment appears in the file. 
        '''
        # --- Import NetPyNE library that loads SWC files
        from netpyne.support import morphology
        # --- Loads the full cell morphology directly from the SWC file to recreate the indexing done using 'afferent_section_pos'
        raw_morphology = morphology.load(morphologyFullPath)
        raw_all_secs_ = raw_morphology.all
        raw_all_secs=[]
        for rsec in raw_all_secs_:
            # --- Extracts and formats the NEURON section name (sec[NUM]) to NetPyNE fomat ('sec_NUM')
            rsec_name_0 = rsec.hname(); rsec_name_1=rsec_name_0.split('['); rsec_name_2=rsec_name_1[1].split(']');  rsec_name = rsec_name_1[0]+'_'+rsec_name_2[0]
            raw_all_secs.append(rsec_name)
        # --- Dictionary that counts and store the number of sections from each type
        rsec_name_dict={'soma':0,'axon':0,'dend':0} # Note: the order of these keys should remain the same, to preserve the indexing SOMA->AXON->DEND(BASAL->APICAL)
        # --- Finds the number of IDs (final_id_) for each section type (key)
        final_id_=0;final_ids=[]
        for key in rsec_name_dict.keys():
            final_id_ = sum(key in rsec_name for rsec_name in raw_all_secs)+final_id_; final_ids.append((key,final_id_))
        # --- Assigns an ID to each section type in the full cell morphology and indexes with a NetPyNE section
        start_id=0;full_ids_list=[]
        for key,final_id in final_ids:
            for sec_ind,id in enumerate(range(start_id,final_id)):full_ids_list.append((key+'_'+str(sec_ind),id))
            start_id = final_id
        # --- Return a list of sections and IDs for the full cell morphology
        return full_ids_list

##################################################################################################################################################################################################################
class Utils():
    def pathDistance(cell_dict,root_sec='soma_0',ignoreAxon=True):
        # --- Creates a dictionary where the key:value are (parent Section, section):<path distance between center points (parent section center to connecting extremity + current section center to connecting extremity)>
        map_secs={}
        for sec in cell_dict['secs'].keys():
            if ignoreAxon:
                if 'axon' in sec:continue

            if 'parentSec' in cell_dict['secs'][sec]['topol'].keys():
                parSec=cell_dict['secs'][sec]['topol']['parentSec']
                map_secs.update({(parSec,sec):  cell_dict['secs'][parSec]['geom']['L']/2+
                                                cell_dict['secs'][sec]['geom']['L']/2})

        # --- Creates a dictionary where the key:value (section:<path distance to soma>), calculated using the previous sections
        soma_pathDist={root_sec:0}
        failed=[]
        for (parSec,sec) in map_secs.keys():
            if parSec == root_sec:soma_pathDist.update({sec:map_secs[(parSec,sec)]})
            else:
                try:    soma_pathDist.update({sec:map_secs[(parSec,sec)]+soma_pathDist[parSec]})
                except: failed.append(sec)
        if len(failed)>0:print('Failed sections: ',failed)
        return soma_pathDist
    
    def findBasalDends(cell_dict,root_sec='soma_0',ignoreAxon=True):
        basalDends=[]
        for sec in cell_dict['secs'].keys():
            if ignoreAxon:
                if 'axon' in sec:continue
            if 'parentSec' in cell_dict['secs'][sec]['topol'].keys():
                if root_sec == cell_dict['secs'][sec]['topol']['parentSec']:basalDends.append(sec)
        return basalDends
    
    def secListFromPathDistance(soma_pathDist,distThresh=[10,75],basalDendrites=[]):
        secLists_dict={
                        'inputs__proximal':[],
                        'inputs__intermediate':[],
                        'inputs__distal':[],
                        }
        for sec in soma_pathDist.keys():
            if 'soma' in sec:                       secLists_dict['inputs__proximal'].append(sec)
            elif sec in basalDendrites:             secLists_dict['inputs__proximal'].append(sec)
            elif soma_pathDist[sec]< distThresh[0]: secLists_dict['inputs__proximal'].append(sec)
            elif soma_pathDist[sec]>=distThresh[0] \
              and soma_pathDist[sec]<distThresh[1]: secLists_dict['inputs__intermediate'].append(sec)
            elif soma_pathDist[sec]>=distThresh[1]: secLists_dict['inputs__distal'].append(sec)
            else: continue
        return secLists_dict
    
    # --- Recursive function to replace the keys of a dictionary at any nested level 
    # (source: https://stackoverflow.com/questions/38491318/replace-keys-in-a-nested-dictionary)
    # old_dict:     dictionary to be modified
    # key_dict:     dictionary with <current keys> as <keys>, and <new keys> as <values> (e.g.: {old_key1:new_key1, old_key2:new_key2})
    def replace_keys(old_dict, key_dict):
        new_dict = { }
        for key in old_dict.keys():
            new_key = key_dict.get(key, key)
            if isinstance(old_dict[key], dict): new_dict[new_key] = Utils.replace_keys(old_dict[key], key_dict)
            else:                               new_dict[new_key] = old_dict[key]
        return new_dict

    #define moving average function
    def moving_avg(x, n):
        cumsum = np.cumsum(np.insert(x, 0, 0)) 
        return (cumsum[n:] - cumsum[:-n]) / float(n)

    def calculateSynapticDrive(netParams, source_pops, duration, movingAvgWindow=25, plotFig=False):
        if type(netParams)==dict:
            popParams       = netParams['popParams']
            synMechParams   = netParams['synMechParams']
            connParams      = netParams['connParams']
        else:    
            popParams       = netParams.popParams
            synMechParams   = netParams.synMechParams
            connParams      = netParams.connParams
            
        store_syn_drive = {}    # --- Store synaptic drive (spike * weight)
        for source_pop in source_pops:
            print('Storing synaptic drive for ', source_pop)
            
            edge_indexes=[]
            for mechName in synMechParams.keys():
                if source_pop in mechName:
                    edge_indexes.append(mechName.split('__')[1])
            
            # --- Array of zeros with time step == cfg.dt
            store_syn_drive.update({source_pop:[0 for a in np.arange(0,duration)]})  # --- Store synaptic drive (spike * weight)
            for prePop in popParams.keys():
                if prePop.split('__')[0] in source_pop:
                    spike_times     = popParams[prePop]['spkTimes']
                    spike_times_int = [int(spkt) for spkt in spike_times if spkt<duration]

                    for edge_index in edge_indexes:
                        connName    = prePop.split('__')[0]+'__'+edge_index+'__vecstim_conn'
                        synMechName = prePop.split('__')[0]+'__'+edge_index+'__vecstim_mech'
                        if 'GABA' in synMechParams[synMechName]['mod']: weightValence = -1 # inhibitory mech
                        else:                                           weightValence =  1 # excitatory mech
                        try:    synConductance = synMechParams[synMechName]['gmax']
                        except: synConductance = 1
                        connWeight   = connParams[connName]['weight']

                        for spkt in spike_times_int:
                            store_syn_drive[source_pop][spkt]+=(weightValence*connWeight*synConductance)
        if plotFig:
            legend_list=[]
            color_list=[]
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

            import matplotlib.pyplot as plt
            # plt.figure(figsize=(20,10))
            for source_pop in source_pops:
                print('Plotting synaptic drive for ', source_pop)
                color = pathway_colors_dict[source_pop]
                plt.plot(store_syn_drive[source_pop],color,alpha=0.3)
                plt.plot(Utils.moving_avg(store_syn_drive[source_pop], 25),color)

                legend_list.append(source_pop)
                color_list.append(color)
            # plt.legend(legend_list)
            
            # --- Handles for plot legend
            handles = []; labels = []
            for pathway in pathway_colors_dict.keys(): 
                if pathway in source_pops:
                    if 'thalamus_neurons|' in pathway: pathway.replace('thalamus_neurons|','')
                    handles.append(plt.Line2D([], [], color=pathway_colors_dict[pathway], marker="o", linewidth=0))
                    labels.append(pathway)
            plt.legend(handles=handles,labels=labels,loc='lower right')
            # plt.savefig('syn_drive_test.png')
        return store_syn_drive

##################################################################################################################################################################################################################

class ConvertSynapses():
    def getTargetSections(cell,node_pathway_edges,full_ids_list):
        SONATA_dict = {'SONATA_sec_id':{},'SONATA_id_sec':{}}
        # cell['secLists'].update({'SONATA_sec_id':{},'SONATA_id_sec':{}})
        for sec,id in full_ids_list:
            # --- Only creates IDs for the secs that are present in the model (IDs for removed axonal sections are ignored)
            if sec in cell['secs'].keys():
                SONATA_dict['SONATA_sec_id'].update({sec:id})
                SONATA_dict['SONATA_id_sec'].update({id:sec})
            
        import pandas as pd
        for edge_name in node_pathway_edges.keys():
            edge_indexes            = node_pathway_edges[edge_name].index
            afferent_section_ids    = node_pathway_edges[edge_name]['afferent_section_id']
            secList = [SONATA_dict['SONATA_id_sec'][aff_sec_id] for aff_sec_id in afferent_section_ids]
            df_sec = pd.DataFrame(secList, columns=['sec'],index=edge_indexes)    
            for column in list(df_sec.columns):
                if column not in node_pathway_edges[edge_name].columns:
                    value=df_sec[column].values
                    node_pathway_edges[edge_name].insert(len(node_pathway_edges[edge_name].columns),column=column,value=value)   
        return node_pathway_edges
    
    def checkSecs(cell,cell_properties,node_pathway_edges,full_ids_list):
        # --- Checks if it needs to  create the list of sections before attempting to plot
        remap_secs=False
        for edge_name in node_pathway_edges.keys():
            if 'sec' not in list(node_pathway_edges[edge_name].columns):
                remap_secs=True
                break
        if remap_secs: print('missing data - auto-generating secs'); node_pathway_edges = ConvertSynapses.getTargetSections(cell,node_pathway_edges,full_ids_list)
        return node_pathway_edges

    def getInverseOrientation(cell_properties):
        from scipy.spatial.transform import Rotation as R
        # --- Cell orientation is based on Quaternions
        orientation = R.from_quat([cell_properties['orientation_x'],cell_properties['orientation_y'],cell_properties['orientation_z'],cell_properties['orientation_w']])
        # --- Obtaining the inverse orientation to reposition synapses on top of the cell morphology
        orientation_inv = orientation.inv()
        return orientation_inv

    def convert3DLocation(cell,cell_properties,node_pathway_edges):
        print('Converting Synapse Locations')
        # --- Package imports
        import math
        import pandas as pd
        # from scipy.spatial.transform import Rotation as R

        # --- Reference position within cell_properties to match the synapses to
        cell_3d_position = [cell_properties['x'],cell_properties['y'],cell_properties['z']] 
        # # --- Cell orientation is based on Quaternions
        # orientation = R.from_quat([cell_properties['orientation_x'],cell_properties['orientation_y'],cell_properties['orientation_z'],cell_properties['orientation_w']])
        # # --- Obtaining the inverse orientation to reposition synapses on top of the cell morphology
        # orientation_inv = orientation.inv()

        orientation_inv = ConvertSynapses.getInverseOrientation(cell_properties)

        # --- Iterate over different edge sources
        for edge_name in node_pathway_edges.keys():
            edge_indexes                = node_pathway_edges[edge_name].index
            afferent_positions          = node_pathway_edges[edge_name][['afferent_center_x','afferent_center_y','afferent_center_z']]
            afferent_positions_         = afferent_positions.sub(cell_3d_position, axis='columns')
            afferent_positions_inverted = orientation_inv.apply(afferent_positions_)
            df_afferent_positions_inverted = pd.DataFrame(afferent_positions_inverted, columns=['x_remap','y_remap','z_remap'],index=edge_indexes)    
            for column in list(df_afferent_positions_inverted.columns):
                if column not in node_pathway_edges[edge_name].columns:
                    value=df_afferent_positions_inverted[column].values
                    node_pathway_edges[edge_name].insert(len(node_pathway_edges[edge_name].columns),column=column,value=value)
        # --- Returns updated version of <node_pathway_edges> including remapped synaptic locations
        print('Synapse Locations converted')
        return node_pathway_edges

    def checkSyns(cell,cell_properties,node_pathway_edges):
        # --- Checks if it needs to remap the synapses before attempting to plot
        remap_syns=False
        for edge_name in node_pathway_edges.keys():
            if ('x_remap' and 'y_remap' and 'z_remap') not in list(node_pathway_edges[edge_name].columns):
                remap_syns=True
                break
        if remap_syns: print('missing data - auto-remapping syns'); node_pathway_edges = ConvertSynapses.convert3DLocation(cell,cell_properties,node_pathway_edges)
        return node_pathway_edges
    def plotCellShape(cell,extra_flag='_cellShape'):
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches

        try:cell_label = str(cell['conds']['cellType'])
        except:cell_label = ''

        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        # --- Plot Soma center - SWC files are loaded with cells centered at (0,0,0)
        # ax.scatter(0,0,0, marker='o',color='r',linewidths=0, alpha=1)

        # --- Plot 3d locations of cell morphology sections
        x3d=[];     y3d=[];     z3d=[];     d3d=[]
        for sec in cell['secs'].keys():
            if 'soma' in sec: 
                pt3d_diams = [pt3d[3] for ind_pt3d, pt3d in enumerate(cell['secs'][sec]['geom']['pt3d'])]
                ax.scatter(0,0,0, s= max(pt3d_diams), marker='o',color='k', alpha=1)
                continue
            if 'pt3d' not in cell['secs'][sec]['geom'].keys():continue
            if len(cell['secs'][sec]['geom']['pt3d'])>0 and type(cell['secs'][sec]['geom']['pt3d'])==list:
                for ind_pt3d, pt3d in enumerate(cell['secs'][sec]['geom']['pt3d']):
                    if ind_pt3d==0:
                        parent_sec = cell['secs'][sec]['topol']['parentSec']
                        previous_pt3d = cell['secs'][parent_sec]['geom']['pt3d'][-1]
                    else:
                        previous_pt3d = cell['secs'][sec]['geom']['pt3d'][ind_pt3d-1]
                    ax.plot([previous_pt3d[0],pt3d[0]],[previous_pt3d[1],pt3d[1]],[previous_pt3d[2],pt3d[2]],linewidth=pt3d[3],color='k', alpha=0.7)


        ax.set_xlabel('X Label');ax.set_ylabel('Y Label');ax.set_zlabel('Z Label')
        ax.view_init(90, 0)
        ax.set_xlim(-150,150)
        ax.set_ylim(-200,200)
        ax.set_axis_off() 
        plt.savefig(Global.figsFolder+'Conversion_cellSynapses_'+cell_label+extra_flag+'.png',dpi=1000)

    def plotSynapseLocation(cell,cell_properties,node_pathway_edges,full_ids_list, extra_flag='',showEstimadedSyns=True,plotSimplified=True):
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches
        # --- Checks if it needs to  create the list of sections before attempting to plot
        node_pathway_edges = ConvertSynapses.checkSecs(cell,cell_properties,node_pathway_edges,full_ids_list)
        # --- Checks if it needs to remap the synapses before attempting to plot
        node_pathway_edges = ConvertSynapses.checkSyns(cell, cell_properties, node_pathway_edges)
        
        try:cell_label = str(cell['conds']['cellType'])
        except:cell_label = ''

        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        # --- Plot Soma center - SWC files are loaded with cells centered at (0,0,0)
        ax.scatter(0,0,0, marker='o',color='cyan',linewidths=0, alpha=1)
        
        # --- Plot 3d locations of cell morphology sections
        x3d=[];     y3d=[];     z3d=[]
        for sec in cell['secs'].keys():
            if 'pt3d' not in cell['secs'][sec]['geom'].keys():continue
            if len(cell['secs'][sec]['geom']['pt3d'])>0 and type(cell['secs'][sec]['geom']['pt3d'])==list:
                for pt3d in cell['secs'][sec]['geom']['pt3d']:
                    x3d.append(pt3d[0]);    y3d.append(pt3d[1]);    z3d.append(pt3d[2])
        # --- sec/pt3d points
        ax.scatter(x3d, y3d, z3d, marker='o',color='k',s=0.25,linewidths=0, alpha=0.1)
        # --- remapped syn locations
        for edge_name in node_pathway_edges.keys():
            if   'MedialLemniscus_projections'  in edge_name:   c_ = 'darkorange'
            elif 'thalamus_neurons|VPL_TC'      in edge_name:   c_ = 'r'
            elif 'thalamus_neurons|VPL_IN'      in edge_name:   c_ = 'hotpink'
            elif 'thalamus_neurons|Rt_RC'       in edge_name:   c_ = 'slateblue'
            elif 'CorticoThalamic_projections'  in edge_name:   c_ = 'limegreen'
            else:                                               c_ = 'deepskyblue'
            ax.scatter( node_pathway_edges[edge_name]['x_remap'], node_pathway_edges[edge_name]['y_remap'],node_pathway_edges[edge_name]['z_remap'], 
                        marker='o',color=c_,s=1.0)

        ax.set_xlabel('X Label');ax.set_ylabel('Y Label');ax.set_zlabel('Z Label')
        ax.view_init(90, 0)
        # k_patch = mpatches.Patch(color='k',     label='sec pt3d points')
        # b_patch = mpatches.Patch(color='b',     label='Remapped synapses')
        # r_patch = mpatches.Patch(color='r',     label='Estimated section')
        # plt.legend(handles=[k_patch,b_patch,r_patch])
        plt.savefig(Global.figsFolder+'Conversion_cellSynapses_'+cell_label+extra_flag+'.png',dpi=1000)

        # --- Simplified plot - with no pt3d details
        if plotSimplified:
            fig = plt.figure()
            ax = fig.add_subplot(projection='3d')
            
            validate_x3d_lines=[];  validate_y3d_lines=[];  validate_z3d_lines=[]
            for edge_name in node_pathway_edges.keys():
                validate_x3d_=[];       validate_y3d_=[];       validate_z3d_=[]
                sec_list                    = list(node_pathway_edges[edge_name]['sec'])                    # estimated section
                afferent_section_pos_list   = list(node_pathway_edges[edge_name]['afferent_section_pos'])   # syn location from edge information
                for ind,sec in enumerate(sec_list):
                    afferent_section_pos =  afferent_section_pos_list[ind]   # syn location from edge information

                    pt3d_i=cell['secs'][sec]['geom']['pt3d'][0];    pt3d_f=cell['secs'][sec]['geom']['pt3d'][-1]
                    pt3d_d = [] # 1st and last pt3d distance
                    for i, f in zip(pt3d_i,pt3d_f):pt3d_d.append(f - i)
                    pt3d_g=[]   # loc equivalent distance
                    for i in range(3):pt3d_g.append((afferent_section_pos*(pt3d_f[i]-pt3d_i[i]))+pt3d_i[i])
                    validate_x3d_.append(pt3d_g[0]);    validate_y3d_.append(pt3d_g[1]);    validate_z3d_.append(pt3d_g[2])

                    if [pt3d_i[0],pt3d_f[0]] not in validate_x3d_lines: # avoid plotting the same line multiple times 
                        validate_x3d_lines.append([pt3d_i[0],pt3d_f[0]]);validate_y3d_lines.append([pt3d_i[1],pt3d_f[1]]);validate_z3d_lines.append([pt3d_i[2],pt3d_f[2]])
                    
                if   'MedialLemniscus_projections'  in edge_name:   c_ = 'darkorange'
                elif 'thalamus_neurons|VPL_TC'      in edge_name:   c_ = 'r'
                elif 'thalamus_neurons|VPL_IN'      in edge_name:   c_ = 'hotpink'
                elif 'thalamus_neurons|Rt_RC'       in edge_name:   c_ = 'slateblue'
                elif 'CorticoThalamic_projections'  in edge_name:   c_ = 'limegreen'
                else:                                               c_ = 'deepskyblue'
                ax.scatter(validate_x3d_, validate_y3d_, validate_z3d_, marker='o',color=c_,s=0.5)

            # ax.scatter(validate_x3d_, validate_y3d_, validate_z3d_, marker='o',color='g',s=0.5)
            # --- Plot morphology without curving based on 3d points
            for i in range(len(validate_x3d_lines)):
                ax.plot(validate_x3d_lines[i],validate_y3d_lines[i],validate_z3d_lines[i],color='k',alpha=0.1)

            ax.set_xlabel('X Label');ax.set_ylabel('Y Label');ax.set_zlabel('Z Label')
            ax.view_init(90, 0)
            # g_patch = mpatches.Patch(color='g',     label='Validated section')
            # y_patch = mpatches.Patch(color='grey',  label='Validated section')
            # plt.legend(handles=[g_patch,y_patch])
            plt.savefig(Global.figsFolder+'Conversion_simplified_'+cell_label+extra_flag+'.png',dpi=1000)
        
        plt.show()


    def reorderSynapseLocation(cell,cell_properties,node_pathway_edges,full_ids_list,plotFigure=False):
        import math
        # --- Checks if it needs to  create the list of sections before attempting to plot
        node_pathway_edges = ConvertSynapses.checkSecs(cell,cell_properties,node_pathway_edges,full_ids_list)
        # --- Checks if it needs to remap the synapses before attempting to plot
        node_pathway_edges = ConvertSynapses.checkSyns(cell, cell_properties, node_pathway_edges)
        
        # --- Reorder indexes based on distance
        store_syn_to_soma_distances=[]
        for edge_name in node_pathway_edges.keys():
            edge_indexes = node_pathway_edges[edge_name].index
            for edge_index in edge_indexes:
                edge = node_pathway_edges[edge_name].loc[edge_index]
                syn_to_soma_distance = math.sqrt((edge['x_remap'])**2+(edge['y_remap'])**2+(edge['z_remap'])**2)                # obs: soma position = (0,0,0)
                store_syn_to_soma_distances.append((syn_to_soma_distance,edge_index,edge))
        # --- Sorts the array based on soma distance
        store_syn_to_soma_distances.sort()

        # --- Ranking edge_names (1st: ML, 2nd: IN, 3rd: RT, 4th: CT)
        rank_edge_names=[[],[],[],[],[],[]]
        for edge_name in node_pathway_edges.keys():
            if   'MedialLemniscus_projections'  in edge_name:   rank=0
            elif 'thalamus_neurons|VPL_TC'      in edge_name:   rank=1
            elif 'thalamus_neurons|VPL_IN'      in edge_name:   rank=2
            elif 'thalamus_neurons|Rt_RC'       in edge_name:   rank=3
            elif 'CorticoThalamic_projections'  in edge_name:   rank=4
            else:                                               rank=5
            rank_edge_names[rank].append(edge_name)

        # --- Create a list of ranked edge names to use as dict keys
        ranked_edge_names=[]
        node_pathway_edges_r={}
        for rank in rank_edge_names:
            for edge_name_r in rank:
                ranked_edge_names.append(edge_name_r)
                node_pathway_edges_r.update({edge_name_r:node_pathway_edges[edge_name_r].copy()}) # it has to be a COPY of the dataframe, so that the original is preserved
        # print(node_pathway_edges_r)
              
        change_parameters   = [ 'sec','x_remap','y_remap','z_remap','afferent_section_pos',
                                # 'afferent_center_x','afferent_center_y','afferent_center_z',
                                # 'afferent_section_id','afferent_section_type',
                                # 'afferent_segment_id','afferent_segment_offset',
                                # 'afferent_surface_x','afferent_surface_y','afferent_surface_z'
                                ]
        # all_parameters      = node_pathway_edges[edge_name].columns
        # keep_parameters     = [par for par in all_parameters if par not in change_parameters];keep_parameters.sort()
        # edge_change         = node_pathway_edges[edge_name].loc[edge_index][change_parameters]
        # edge_keep           = node_pathway_edges[edge_name].loc[edge_index][keep_parameters]

        resample_stored_syn=0
        for edge_name_r in node_pathway_edges_r.keys():
            if (len(store_syn_to_soma_distances)-resample_stored_syn)%100==0:print('\t | \t syns left: ', len(store_syn_to_soma_distances)-resample_stored_syn)
            # print(edge_name_r)
            edge_indexes_r = node_pathway_edges_r[edge_name_r].index
            for edge_index_r in edge_indexes_r:
                edge_all_r  = node_pathway_edges_r[edge_name_r].loc[edge_index_r]
                new_syn_dist = store_syn_to_soma_distances[resample_stored_syn][0]
                new_syn_edge = store_syn_to_soma_distances[resample_stored_syn][1]
                new_syn_position = store_syn_to_soma_distances[resample_stored_syn][2][change_parameters]
                # --- Changing values one by one
                for col in change_parameters:
                    node_pathway_edges_r[edge_name_r].at[edge_index_r,col] = new_syn_position[col]
                resample_stored_syn+=1
        
        if plotFigure: 
            # ConvertSynapses.plotSynapseLocation(cell, cell_properties, node_pathway_edges,   full_ids_list, extra_flag='_original',     showEstimadedSyns=False, plotSimplified=False)
            ConvertSynapses.plotSynapseLocation(cell, cell_properties, node_pathway_edges_r, full_ids_list, extra_flag='_repositioned', showEstimadedSyns=False, plotSimplified=False)

        return node_pathway_edges_r

    def replaceSynParams(cell_properties,node_pathway_edges):
        TableS1 = StoreParameters.getTableS1()
        print('Converting edge values to TableS1')
        for edge_name in node_pathway_edges.keys():
            # --- Getting edge name properties
            edge_type,edge_mech,edge_source,edge_target = edge_name.split('__')
            pre_pop=edge_source
            post_pop = cell_properties['mtype']
            standard_edge_dict={}
            # --- Loops through the properties of the TableS1
            # print(node_pathway_edges[edge_name].columns)
            for key in TableS1[pre_pop][post_pop]['paper_reference_values'].keys():
                # print(key)
                # --- Checks if it matches the properties in the edges columns
                if key in list(node_pathway_edges[edge_name].columns): standard_edge_dict.update({key:TableS1[pre_pop][post_pop]['paper_reference_values'][key]})
                # else: print('miss')
            # --- Gets the indexes of the edges
            edge_indexes = node_pathway_edges[edge_name].index
            for edge_index in edge_indexes:
                # --- Reasign new values
                for col in standard_edge_dict.keys():
                    node_pathway_edges[edge_name].at[edge_index,col] = standard_edge_dict[col]
        return node_pathway_edges

##################################################################################################################################################################################################################

class LoadBBPCircuit():
    def loadModel(cfg_file):
        # --- BBP package - load the circuit and store the node population
        import bluepysnap
        from bluepysnap import morph
        circuit = bluepysnap.Circuit(cfg_file)
        return circuit

    # --- Get information from a single node
    def getNode(circuit,node_name):
        node_population = circuit.nodes[node_name]
        node_property_names_list = list(circuit.nodes[node_name].property_names)
        node_dataFrame=circuit.nodes[node_name].get(properties=node_property_names_list)
        node_dict={node_population.name:node_dataFrame}
        return node_dict

    # --- Get information from all nodes
    def getNodes(circuit):
        nodes_dict={}
        for node_name in circuit.nodes:
            node_dict=LoadBBPCircuit.getNode(circuit,node_name)
            nodes_dict.update(node_dict)
        return nodes_dict

    # def getNodes(circuit):
    #     nodes_dict={}
    #     for node_name in circuit.nodes:
    #         node_population = circuit.nodes[node_name]
    #         node_property_names_list = list(circuit.nodes[node_name].property_names)
    #         node_dataFrame=circuit.nodes[node_name].get(properties=node_property_names_list)
    #         nodes_dict.update({node_population.name:node_dataFrame})
    #     return nodes_dict

    def getNodeNames(circuit):
        node_names=[]
        for node_name in circuit.nodes:
            node_population = circuit.nodes[node_name]
            node_names.append(node_population.name)
        return node_names
    
    def getEdgeNames(cfg_file):
        circuit     = LoadBBPCircuit.loadModel(cfg_file=cfg_file)
        edge_names=[]
        for edge_name in circuit.edges:
            edge_population = circuit.edges[edge_name]
            edge_names.append(edge_population.name)
        return edge_names

    def getDataFrames(cfg_file, microcircuit_number=2):
        # --- Get dictionary of nodes
        circuit     = LoadBBPCircuit.loadModel(cfg_file=cfg_file)
        nodes_dict  = LoadBBPCircuit.getNodes(circuit)

        # --- Select microcircuit
        circuit_dict={}
        # circuit_dict.update({'thalamus_neurons': nodes_dict['thalamus_neurons']})
        population_names = circuit.nodes.population_names
        print(population_names)
        for pop_name in population_names: circuit_dict.update({pop_name: nodes_dict[pop_name]})
        
        microcircuits_list = list(set(circuit_dict['thalamus_neurons']['region'])); microcircuits_list.sort()
        circuit_dict.update({'microcircuits_list': microcircuits_list})

        mc_VPL  = 'mc'+str(microcircuit_number)+';VPL'
        mc_RT   = 'mc'+str(microcircuit_number)+';Rt'
        mc_num  = 'mc'+str(microcircuit_number)

        # --- Extracting subpopulation information
        circuit_dict.update({mc_num:{}})
        circuit_dict[mc_num].update({'VPL_TC': circuit_dict['thalamus_neurons'].loc[(circuit_dict['thalamus_neurons']['mtype'] == 'VPL_TC') & (circuit_dict['thalamus_neurons']['region'] == mc_VPL)]})
        circuit_dict[mc_num].update({'VPL_IN': circuit_dict['thalamus_neurons'].loc[(circuit_dict['thalamus_neurons']['mtype'] == 'VPL_IN') & (circuit_dict['thalamus_neurons']['region'] == mc_VPL)]})
        circuit_dict[mc_num].update({'Rt_RC':  circuit_dict['thalamus_neurons'].loc[(circuit_dict['thalamus_neurons']['mtype'] == 'Rt_RC') &  (circuit_dict['thalamus_neurons']['region'] == mc_RT)]})

        # --- Extracting subtype information
        circuit_dict[mc_num].update({'VPL_TC_dAD_ltb':     circuit_dict['thalamus_neurons'].loc[(circuit_dict['thalamus_neurons']['mtype'] == 'VPL_TC') & (circuit_dict['thalamus_neurons']['region'] == mc_VPL) & (circuit_dict['thalamus_neurons']['etype'] == 'dAD_ltb')]})
        circuit_dict[mc_num].update({'VPL_TC_dNAD_ltb':    circuit_dict['thalamus_neurons'].loc[(circuit_dict['thalamus_neurons']['mtype'] == 'VPL_TC') & (circuit_dict['thalamus_neurons']['region'] == mc_VPL) & (circuit_dict['thalamus_neurons']['etype'] == 'dNAD_ltb')]})
        circuit_dict[mc_num].update({'Rt_RC_cAD_noscltb':  circuit_dict['thalamus_neurons'].loc[(circuit_dict['thalamus_neurons']['mtype'] == 'Rt_RC')  & (circuit_dict['thalamus_neurons']['region'] == mc_RT)  & (circuit_dict['thalamus_neurons']['etype'] == 'cAD_noscltb')]})
        circuit_dict[mc_num].update({'Rt_RC_cNAD_noscltb': circuit_dict['thalamus_neurons'].loc[(circuit_dict['thalamus_neurons']['mtype'] == 'Rt_RC')  & (circuit_dict['thalamus_neurons']['region'] == mc_RT)  & (circuit_dict['thalamus_neurons']['etype'] == 'cNAD_noscltb')]})
        
        return circuit_dict

    def getDataFrames_new(cfg_file):
        # --- Get dictionary of nodes
        circuit     = LoadBBPCircuit.loadModel(cfg_file=cfg_file)
        nodes_dict  = LoadBBPCircuit.getNodes(circuit)

        # --- Select microcircuit
        circuit_dict={}
        # circuit_dict.update({'thalamus_neurons': nodes_dict['thalamus_neurons']})
        population_names = circuit.nodes.population_names
        for pop_name in population_names: circuit_dict.update({pop_name: nodes_dict[pop_name]})
        
        microcircuits_list = list(set(circuit_dict['thalamus_neurons']['region'])); microcircuits_list.sort()
        circuit_dict.update({'microcircuits_list': microcircuits_list})

        mc_nums_=[]
        for mc in microcircuits_list:
            mc_ = mc.split(';')
            mc_nums_.append(mc_[0])
        mc_nums=list(set(mc_nums_));mc_nums.sort()

        # --- Iterates over all microcircuits
        for mc_num in mc_nums:
            mc_VPL=mc_num+';VPL'
            mc_RT=mc_num+';Rt'
            # --- Extracting subpopulation information
            circuit_dict.update({mc_num:{}})
            circuit_dict[mc_num].update({'VPL_TC': circuit_dict['thalamus_neurons'].loc[(circuit_dict['thalamus_neurons']['mtype'] == 'VPL_TC') & (circuit_dict['thalamus_neurons']['region'] == mc_VPL)]})
            circuit_dict[mc_num].update({'VPL_IN': circuit_dict['thalamus_neurons'].loc[(circuit_dict['thalamus_neurons']['mtype'] == 'VPL_IN') & (circuit_dict['thalamus_neurons']['region'] == mc_VPL)]})
            circuit_dict[mc_num].update({'Rt_RC':  circuit_dict['thalamus_neurons'].loc[(circuit_dict['thalamus_neurons']['mtype'] == 'Rt_RC') &  (circuit_dict['thalamus_neurons']['region'] == mc_RT)]})

            # --- Extracting subtype information
            circuit_dict[mc_num].update({'VPL_TC_dAD_ltb':     circuit_dict['thalamus_neurons'].loc[(circuit_dict['thalamus_neurons']['mtype'] == 'VPL_TC') & (circuit_dict['thalamus_neurons']['region'] == mc_VPL) & (circuit_dict['thalamus_neurons']['etype'] == 'dAD_ltb')]})
            circuit_dict[mc_num].update({'VPL_TC_dNAD_ltb':    circuit_dict['thalamus_neurons'].loc[(circuit_dict['thalamus_neurons']['mtype'] == 'VPL_TC') & (circuit_dict['thalamus_neurons']['region'] == mc_VPL) & (circuit_dict['thalamus_neurons']['etype'] == 'dNAD_ltb')]})
            circuit_dict[mc_num].update({'Rt_RC_cAD_noscltb':  circuit_dict['thalamus_neurons'].loc[(circuit_dict['thalamus_neurons']['mtype'] == 'Rt_RC')  & (circuit_dict['thalamus_neurons']['region'] == mc_RT)  & (circuit_dict['thalamus_neurons']['etype'] == 'cAD_noscltb')]})
            circuit_dict[mc_num].update({'Rt_RC_cNAD_noscltb': circuit_dict['thalamus_neurons'].loc[(circuit_dict['thalamus_neurons']['mtype'] == 'Rt_RC')  & (circuit_dict['thalamus_neurons']['region'] == mc_RT)  & (circuit_dict['thalamus_neurons']['etype'] == 'cNAD_noscltb')]})
        
        return circuit_dict

    def getNodeProperties(circuit_dict,node_id):
        # # --- Get dictionary of nodes
        # circuit     = LoadBBPCircuit.loadModel(cfg_file=cfg_file)
        # nodes_dict  = LoadBBPCircuit.getNodes(circuit)
        # # --- Select microcircuit
        # circuit_dict={}
        # circuit_dict.update({'thalamus_neurons': nodes_dict['thalamus_neurons']})
        cell_properties = circuit_dict['thalamus_neurons'].loc[node_id]
        return cell_properties

    # --- Method to run stats on the pathway_edges properties to compare with reported values
    def evaluateEdgeParams(pathway_edges):
        print('\n# of edges: ', len(pathway_edges))
        edge_params = list(pathway_edges.columns)
        edge_params.sort()
        for edge_param in edge_params:
            edge_param_mean = np.mean(pathway_edges[edge_param])
            edge_param_std  = np.std(pathway_edges[edge_param])
            print(edge_param,':',edge_param_mean,' +/-',edge_param_std)

    def getNodeEdges(cfg_file,node_id,special_conditions={
        'inter_sources':['Rt_RC'],
        'inter_targets':['Rt_RC'],
        # 'inter_conns':[ ('Rt_RC','VPL_TC'),
        #                 ('VPL_TC','Rt_RC'),
        #                 ('Rt_RC','Rt_RC')
        #                 ],
        # 'remove_recurrent_conns':True,
        'select_microcircuit':None
        }):

        # --- Name of the target (supports GID <int> and dict <{'population': edge_target_pop, 'mtype': t_mtype}>)
        if type(node_id)==dict:
            try:    node_id_name = node_id['population']+'|'+node_id['mtype']
            except: node_id_name = ''
        else:
            try:        node_id_name = str(node_id)
            except:     node_id_name = ''

        print('node_id_name: ',node_id_name)
        # --- Get circuit edges names
        circuit     = LoadBBPCircuit.loadModel(cfg_file=cfg_file)
        edge_names  = LoadBBPCircuit.getEdgeNames(cfg_file=cfg_file)

        # # --- Dictionary structure to store the edges
        node_pathway_edges={}

        for edge_name_ind,edge_name in enumerate(edge_names):
            print('\n________________________________________________________________________________________________')
            print('Processing Edge: \t', edge_name)
            edge_population = circuit.edges[edge_name]
            edge_property_names_list = list(circuit.edges[edge_name].property_names)
            
            edge_source_pop = edge_population.source.name
            edge_target_pop = edge_population.target.name
            s,t,edge_mech=edge_name.split('__')

            # --- Process internal connectivity
            if edge_source_pop==edge_target_pop:
                edge_type = 'internal'
            
                # --- Node id as target
                source_node_dict = LoadBBPCircuit.getNode(circuit,edge_source_pop)
                source_mtypes = list(set(source_node_dict[edge_source_pop]['mtype']))
                for s_mtype in source_mtypes:
                    if 'Rt' in s_mtype:mtype_flag='Rt'
                    elif 'VPL' in s_mtype:mtype_flag='VPL'
                    if s_mtype in special_conditions['inter_sources']:
                        conn_source = {'population': edge_source_pop, 'mtype': s_mtype}
                        conn_target = node_id
                    
                        # --- Selects only connections within the same microcircuit (Defalt: None - selects projections from the whole circuit)
                        if special_conditions['select_microcircuit'] is not None:
                            if type(special_conditions['select_microcircuit'])==int:
                                mc_region='mc'+str(special_conditions['select_microcircuit'])+';'+mtype_flag
                                conn_source.update({'region':mc_region})
                        
                        print('conn_source: ', conn_source, '\t conn_target: ',conn_target)
                        pathway_edges_name  = edge_type+'__'+edge_mech+'__'+edge_source_pop+'|'+s_mtype+'__'+node_id_name
                        pathway_edges       = edge_population.pathway_edges(source=conn_source, target=conn_target, properties=edge_property_names_list)
                        if pathway_edges.empty:
                            print('Empty dataframe')
                        else:
                            node_pathway_edges.update({pathway_edges_name:pathway_edges})
                            LoadBBPCircuit.evaluateEdgeParams(pathway_edges)
                
                # --- Node id as source
                target_node_dict = LoadBBPCircuit.getNode(circuit,edge_target_pop)
                target_mtypes = list(set(target_node_dict[edge_target_pop]['mtype']))
                for t_mtype in target_mtypes:
                    if 'Rt' in t_mtype:mtype_flag='Rt'
                    elif 'VPL' in t_mtype:mtype_flag='VPL'
                    if t_mtype in special_conditions['inter_targets']:
                        conn_source = node_id
                        conn_target = {'population': edge_target_pop, 'mtype': t_mtype}
                    
                        # --- Selects only connections within the same microcircuit
                        if special_conditions['select_microcircuit'] is not None:
                            if type(special_conditions['select_microcircuit'])==int:
                                mc_region='mc'+str(special_conditions['select_microcircuit'])+';'+mtype_flag
                                conn_target.update({'region':mc_region})
                        
                        print('conn_source: ', conn_source, '\t conn_target: ',conn_target)

                        pathway_edges_name  = edge_type+'__'+edge_mech+'__'+node_id_name+'__'+edge_target_pop+'|'+t_mtype
                        pathway_edges       = edge_population.pathway_edges(source=conn_source, target=conn_target, properties=edge_property_names_list)
                        if pathway_edges.empty:
                            print('Empty dataframe')
                        else:
                            node_pathway_edges.update({pathway_edges_name:pathway_edges})
                            # node_pathway_edges[edge_type][edge_mech].update({pathway_edges_name:pathway_edges})
                            LoadBBPCircuit.evaluateEdgeParams(pathway_edges)
            # --- Process external connectivity
            else:
                edge_type = 'external'
                conn_source = {'population': edge_source_pop}
                conn_target = node_id

                print('\nConn_source: ', conn_source, '\t conn_target: ',conn_target)
                pathway_edges_name  = edge_type+'__'+edge_mech+'__'+edge_source_pop+'__'+node_id_name
                pathway_edges       = edge_population.pathway_edges(source=conn_source, target=conn_target, properties=edge_property_names_list)
                if pathway_edges.empty:
                    print('Empty dataframe')
                else:
                    node_pathway_edges.update({pathway_edges_name:pathway_edges})
                    # node_pathway_edges[edge_type][edge_mech].update({pathway_edges_name:pathway_edges})
                    LoadBBPCircuit.evaluateEdgeParams(pathway_edges)
        
        return node_pathway_edges

    def saveNodeEdges(node_pathway_edges,edges_fileName):
        import json
        node_pathway_edges_dict={}
        for edge_name in node_pathway_edges.keys():
            node_pathway_edges_dict.update({edge_name:{}})
            node_pathway_edges_dict[edge_name]=node_pathway_edges[edge_name].to_dict(orient='index')
        json_object = json.dumps(node_pathway_edges_dict, indent=4)
        with open(edges_fileName, "w") as outfile: outfile.write(json_object)
        return node_pathway_edges_dict
    
    def loadNodeEdges(edges_fileName):
        print('\t\t- \tLoading edge information')
        import json
        import pandas as pd
        with open(edges_fileName, 'r') as node_pathway_edges_dictObj: node_pathway_edges_dict = json.load(node_pathway_edges_dictObj)        
        node_pathway_edges={}
        for edge_name in node_pathway_edges_dict.keys():df = pd.DataFrame.from_dict(node_pathway_edges_dict[edge_name],orient='index');node_pathway_edges.update({edge_name:df})
        return node_pathway_edges

    def selectPathways(select_pathways,edges):
        edges_={}
        for edge_name in edges.keys():
            edge_type,edge_mech,edge_source_pop,edge_target = edge_name.split('__')
            if edge_source_pop in select_pathways: edges_.update({edge_name:edges[edge_name]})
        edges=edges_
        return edges

##################################################################################################################################################################################################################

class ConfigCfg():
    def sampleGids(pop,cfg_file,gids=None):
        import sys
        import random
        circuit_dict = LoadBBPCircuit.getDataFrames_new(cfg_file=cfg_file)
        
        if  gids is None:
            sample_gids = random.sample(list(circuit_dict['mc2'][pop].index),1)
            print('hi')
        elif type(gids)==int:
            sample_gid = random.sample(list(circuit_dict['mc2'][pop].index),gids)
            # sample_gids = [sample_gid]
            sample_gids = sample_gid

        elif gids == 'test_single':
            if   'VPL_TC' in pop:   sample_gids = [38906]
            elif 'Rt_RC'  in pop:   sample_gids = [30550]
            elif 'VPL_IN' in pop:   sample_gids = [42503]
        elif gids == 'test_reduced':
            if   'VPL_TC' in pop:   sample_gids = [35103, 36243, 37082, 33702, 37625, 41429]
            elif 'Rt_RC'  in pop:   sample_gids = [30550, 30541, 29582, 29174, 29908, 33109]
            elif 'VPL_IN' in pop:   sample_gids = [42503, 42486, 42483, 42470, 42502, 42506]
        elif gids == 'test':
            if   'VPL_TC' in pop:   sample_gids = [35103, 36243, 37082, 33702, 37625, 41429, 35879, 41240, 41615, 34361, 37543, 37177, 41876, 34569, 36963, 41912, 39985, 37055, 36484, 35847, 33798, 34368, 36219, 39232, 34389, 34242, 38102, 35703, 38487, 41067, 37463, 38468, 36711, 34932, 38346, 34503, 36248, 41454, 36721, 33741, 40602, 34274, 41534, 33640, 36881, 34859, 36169, 38276, 37409, 34707, 38440, 41237, 38052, 36302, 33602, 41247, 38036, 39429, 38474, 35824, 38651, 37968, 40213, 42177, 40168, 40215, 41723, 36655, 38134, 41695, 42422, 42460, 36521, 38775, 35220, 35162, 34349, 36440, 35739, 34954, 37256, 41168, 39751, 38748, 33967, 35343, 40876, 39755, 36185, 41399, 39299, 38971, 37093, 37917, 37599, 34471, 39745, 39477, 42073, 36043, 41388, 38169, 34773, 34401, 41379, 37475, 38090, 40659, 37782, 38709, 42405, 41353, 41307, 40641, 37685, 39390, 39239, 35684, 34363, 37548, 36748, 36059, 35158, 40735, 35483, 42198, 34433, 41390, 39229, 40044, 37740, 40122, 36364, 35113, 38793, 40560, 36857, 37553, 41271, 39981, 41439, 38171, 39183, 41890, 37925, 37824, 38002, 35649, 41579, 38806, 37520, 40430, 33822, 39202, 37863, 41253, 33571, 35332, 35748, 39340, 33774, 41571, 42273, 41996, 38098, 36368, 41395, 37033, 39864, 39123, 36611, 40153, 39451, 35662, 42357, 40624, 40363, 36612, 36499, 33806]
            elif 'Rt_RC'  in pop:   sample_gids = [30550, 30541, 29582, 29174, 29908, 33109, 31521, 32579, 32893, 32954, 32696, 32933, 33334, 31927, 30299, 29934, 30694, 31191, 31989, 32369, 30242, 30823, 29379, 31241, 31793, 31492, 32974, 30653, 29993, 30022, 29770, 32501, 29195, 29892, 30730, 30655, 32740, 32640, 28671, 28831, 28660, 29828, 31704, 28988, 29183, 29690, 31254, 30838, 31637, 30922, 30182, 33200, 28663, 31412, 31625, 31778, 29791, 31120, 30543, 29184, 28612, 30652, 32453, 32047, 29522, 32049, 29342, 31907, 30072, 32729, 29735, 32221, 30986, 33224, 31309, 30551, 31296, 29803, 29007, 30947, 28805, 30849, 33463, 29657, 30946, 32631, 31840, 30892, 31646, 31738, 31315, 29086, 29040, 28852, 29608, 30025, 31528, 32662, 32781, 31170, 32479, 33190, 31420, 28785, 30084, 31972, 30225, 30872, 30506, 32036, 33089, 33362, 32299, 32620, 29371, 32292, 32978, 32313, 32267, 30174, 33014, 30007, 31239, 28733, 32470, 31044, 28694, 29087, 29476, 29687, 30990, 29126, 31800, 28834, 31881, 28925, 30252, 29621, 29094, 29304, 31400, 29526, 31674, 32147, 31113, 29861, 32413, 29052, 30152, 29731, 29205, 31864, 31393, 33031, 30772, 28731, 30090, 33325, 30891, 29863, 30403, 31638, 32406, 33043, 30905, 32926, 30014, 30813, 30854, 29679, 29049, 31751, 31816, 29689, 32540, 29846, 28833, 32411, 32730, 29805, 32846, 29328, 30216, 32641, 29663, 30936, 32371, 29722, 31923, 30609, 32591, 30670, 31012, 31181, 33204, 31924, 32040, 28873, 33230, 31602]
            elif 'VPL_IN' in pop:   sample_gids = [42503, 42486, 42483, 42470, 42502, 42506, 42488, 42484, 42495, 42465, 42477, 42485, 42489, 42468, 42480, 42476, 42469, 42464, 42509, 42487, 42507, 42471, 42473, 42474, 42494, 42479, 42505, 42496, 42482, 42490, 42492, 42501, 42510, 42475, 42500, 42499, 42497, 42493, 42472, 42508, 42491, 42467, 42498, 42504, 42466]
        elif    type(gids)==list:
            if gids[0]==int:  sample_gids = gids
            else:
                try:    sample_gids = [int(gid) for gid in gids]
                except: sys.exit('Gids must be integers')
        else:   return

        print(sample_gids)

        for gid in sample_gids:
            cell_properties = LoadBBPCircuit.getNodeProperties(circuit_dict,gid)
            if cell_properties['mtype']!=pop:sys.exit('Gids from many pops')

        return sample_gids

##################################################################################################################################################################################################################      

class CreateNetPyNE():
    def groupEdges(edges):
        print('Grouping edges into a summarized dictionary')
        edges_dict={}
        for edge_name in edges.keys():
            edges_dict.update({edge_name:{}})
            for edge_index in edges[edge_name].index:
                edges_dict[edge_name].update({edge_index:{}})
                edge_data = edges[edge_name].loc[edge_index]
                keys = edge_data.index
                vals = edge_data.values
                for i,k in enumerate(keys):
                    edges_dict[edge_name][edge_index].update({k:vals[i]})
    
        return edges_dict

                # print(edge_index,edge_data)

    def modTemplate(edge_data,mod_file):
        ##################################################################################################################################################################################################################
        # --- Original Prob mod from BBP model
        if mod_file =='ProbAMPANMDA_EMS_original':
            edge_dict={ 'mod':          'ProbAMPANMDA_EMS_original',
                        'Use':          edge_data['u_syn'],
                        'Dep':          edge_data['depression_time'],
                        'Fac':          edge_data['facilitation_time'],
                        'gmax':         edge_data['conductance']*.001,# IMPORTANT: .001 IS A CONVERSION FACTOR FROM nS to uS AND MUST BE ADDED
                        'tau_d_AMPA':   edge_data['decay_time'], 
                        # 'Nrrp':         edge_data['n_rrp_vesicles'], 
                        'NMDA_ratio':   edge_data['NMDA_ratio'], 
                        }
        elif mod_file =='ProbGABAAB_EMS_original':
            edge_dict={ 'mod':          'ProbGABAAB_EMS_original',
                        'Use':          edge_data['u_syn'],
                        'Dep':          edge_data['depression_time'],
                        'Fac':          edge_data['facilitation_time'],
                        'gmax':         edge_data['conductance']*.001,# IMPORTANT: .001 IS A CONVERSION FACTOR FROM nS to uS AND MUST BE ADDED
                        'tau_d_GABAA':  edge_data['decay_time'], 
                        # 'Nrrp':         edge_data['n_rrp_vesicles'], 
                        }
        ##################################################################################################################################################################################################################
        # --- Prob mod
        elif mod_file =='ProbAMPANMDA_EMS':
            edge_dict={ 'mod':          'ProbAMPANMDA_EMS',
                        'Use':          edge_data['u_syn'],
                        'Dep':          edge_data['depression_time'],
                        'Fac':          edge_data['facilitation_time'],
                        'gmax':         edge_data['conductance']*.001,
                        'tau_d_AMPA':   edge_data['decay_time'], 
                        # 'Nrrp':         edge_data['n_rrp_vesicles'], 
                        'NMDA_ratio':   edge_data['NMDA_ratio'], 
                        }
        elif mod_file =='ProbGABAAB_EMS':
            edge_dict={ 'mod':          'ProbGABAAB_EMS',
                        'Use':          edge_data['u_syn'],
                        'Dep':          edge_data['depression_time'],
                        'Fac':          edge_data['facilitation_time'],
                        'gmax':         edge_data['conductance']*.001,
                        'tau_d_GABAA':  edge_data['decay_time'], 
                        # 'Nrrp':         edge_data['n_rrp_vesicles'], 
                        }
        ##################################################################################################################################################################################################################
        # --- Prob S1 mod
        elif mod_file =='ProbAMPANMDA_EMS_S1':
            edge_dict={ 'mod':          'ProbAMPANMDA_EMS_S1',
                        'Use':          edge_data['u_syn'],
                        'Dep':          edge_data['depression_time'],
                        'Fac':          edge_data['facilitation_time'],
                        'gmax':         edge_data['conductance']*.001, # IMPORTANT: .001 IS A CONVERSION FACTOR FROM nS to uS AND MUST BE ADDED
                        'tau_d_AMPA':   edge_data['decay_time'], 
                        'NMDA_ratio':   edge_data['NMDA_ratio'],
                        }
        elif mod_file =='ProbGABAAB_EMS_S1':
            edge_dict={ 'mod':          'ProbGABAAB_EMS_S1',
                        'Use':          edge_data['u_syn'],
                        'Dep':          edge_data['depression_time'],
                        'Fac':          edge_data['facilitation_time'],
                        'gmax':         edge_data['conductance']*.001, # IMPORTANT: .001 IS A CONVERSION FACTOR FROM nS to uS AND MUST BE ADDED
                        'tau_d_GABAA':  edge_data['decay_time'], 
                        }
        ##################################################################################################################################################################################################################
        # --- Prob mod CORENEURON
        elif mod_file =='ProbAMPANMDA_EMS_CORENEURON':
            edge_dict={ 'mod':          'ProbAMPANMDA_EMS_CORENEURON',
                        'Use':          edge_data['u_syn'],
                        'Dep':          edge_data['depression_time'],
                        'Fac':          edge_data['facilitation_time'],
                        'gmax':         edge_data['conductance']*.001,
                        'tau_d_AMPA':   edge_data['decay_time'], 
                        # 'Nrrp':         edge_data['n_rrp_vesicles'], 
                        'NMDA_ratio':   edge_data['NMDA_ratio'], 
                        }
        elif mod_file =='ProbGABAAB_EMS_CORENEURON':
            edge_dict={ 'mod':          'ProbGABAAB_EMS_CORENEURON',
                        'Use':          edge_data['u_syn'],
                        'Dep':          edge_data['depression_time'],
                        'Fac':          edge_data['facilitation_time'],
                        'gmax':         edge_data['conductance']*.001,
                        'tau_d_GABAA':  edge_data['decay_time'], 
                        # 'Nrrp':         edge_data['n_rrp_vesicles'], 
                        }

        ##################################################################################################################################################################################################################
        # --- Prob S1 mod CORENEURON
        elif mod_file =='ProbAMPANMDA_EMS_S1_CORENEURON':
            edge_dict={ 'mod':          'ProbAMPANMDA_EMS_S1_CORENEURON',
                        'Use':          edge_data['u_syn'],
                        'Dep':          edge_data['depression_time'],
                        'Fac':          edge_data['facilitation_time'],
                        # 'g':         edge_data['conductance'],
                        # 'gmax':         edge_data['conductance']*.001, # IMPORTANT: .001 IS A CONVERSION FACTOR FROM nS to uS AND MUST BE ADDED
                        'tau_d_AMPA':   edge_data['decay_time'], 
                        'NMDA_ratio':   edge_data['NMDA_ratio'],
                        }
        elif mod_file =='ProbGABAAB_EMS_S1_CORENEURON':
            edge_dict={ 'mod':          'ProbGABAAB_EMS_S1_CORENEURON',
                        'Use':          edge_data['u_syn'],
                        'Dep':          edge_data['depression_time'],
                        'Fac':          edge_data['facilitation_time'],
                        # 'g':         edge_data['conductance'],
                        # 'gmax':         edge_data['conductance']*.001, # IMPORTANT: .001 IS A CONVERSION FACTOR FROM nS to uS AND MUST BE ADDED
                        'tau_d_GABAA':  edge_data['decay_time'], 
                        }

        ##################################################################################################################################################################################################################
        # --- Det mod
        elif mod_file =='DetAMPANMDA':
            edge_dict={ 'mod':          'DetAMPANMDA',
                        'Use':          edge_data['u_syn'],
                        'Dep':          edge_data['depression_time'],
                        'Fac':          edge_data['facilitation_time'],
                        # 'g':            edge_data['conductance'],
                        # 'gmax':         edge_data['conductance']*.001, # IMPORTANT: .001 IS A CONVERSION FACTOR FROM nS to uS AND MUST BE ADDED
                        'tau_d_AMPA':   edge_data['decay_time'], 
                        'NMDA_ratio':   edge_data['NMDA_ratio'],
                        }
        elif mod_file =='DetGABAAB':
            edge_dict={ 'mod':          'DetGABAAB',
                        'Use':          edge_data['u_syn'],
                        'Dep':          edge_data['depression_time'],
                        'Fac':          edge_data['facilitation_time'],
                        # 'g':            edge_data['conductance'],
                        # 'gmax':         edge_data['conductance']*.001, # IMPORTANT: .001 IS A CONVERSION FACTOR FROM nS to uS AND MUST BE ADDED
                        'tau_d_GABAA':  edge_data['decay_time'], 
                        }
        # --- Do not use
        elif mod_file == 'Gap':
            edge_dict={ 'mod':          'Gap',
                        'g':            edge_data['conductance'], 
                        }

        # --- Use this gap junction instead
        elif mod_file == 'ElectSyn':
            edge_dict=  {   'mod': 'ElectSyn',
                            'g':            edge_data['conductance'], 
                            'pointerParams': {
                                'target_var':  'vpeer',
                                'source_var': 'v', # already there by default:
                                'bidirectional': True # already there by default:
                                }}    
        
        else:
            raise Exception("Mod file not recognized") 
        
        return edge_dict
    
    def modifySyn(circuit_dict,edge_name):
        # Spreadsheet with the analysis of paramters for each pathway
        # https://docs.google.com/spreadsheets/d/1SsnAR_bq_BPaWQPI9gDqif51KUVrtCssAtf8ErUx4Sk/edit?usp=sharing 
        modif_dict={}
        edge_type,edge_mech,edge_source_pop,edge_target = edge_name.split('__')
        edge_target_pop = circuit_dict['thalamus_neurons']['mtype'].loc[int(edge_target)]
        if edge_mech == 'chemical': 
            if ('CorticoThalamic' in edge_source_pop):
                modif_dict = {'mg':1, 'NMDA_ratio': 1.91}
            elif ('MedialLemniscus' in edge_source_pop):
                modif_dict = {'mg':1, 'NMDA_ratio': 0.41}
            else:
                if 'thalamus_neurons' in edge_source_pop:
                    # --- Identifies if the synapse SOURCE is VPL_TC
                    if   edge_source_pop.split('|')[1]=='VPL_TC':   
                        modif_dict = {'mg':1, 'NMDA_ratio': 0.57}
                    elif (edge_source_pop.split('|')[1]=='VPL_IN'):
                        if edge_target_pop == 'VPL_TC':             modif_dict = {'e_GABAA':-94,'e_GABAB':-97, 'GABAB_ratio': 0, 'tau_d_GABAB':77} # 2023-09-05: testing if midification in 'tau_d_GABAB':77 should apply to all synapses, or just to the ones that are not Rt_RC and VPL_IN (if any exist)
                    elif (edge_source_pop.split('|')[1]=='Rt_RC'):
                        if edge_target_pop == 'VPL_TC':             modif_dict = {'e_GABAA':-94,'e_GABAB':-97, 'GABAB_ratio': 0, 'tau_d_GABAB':77}
                        else:                                       modif_dict = {'e_GABAA':-82,'e_GABAB':-97, 'GABAB_ratio': 0, 'tau_d_GABAB':77}
                    else:
                        modif_dict={}
                        print('no matching syn')
                
            # if ('CorticoThalamic' in edge_source_pop) or ('MedialLemniscus' in edge_source_pop):
            #     modif_dict = {'mg':1}
            #     # modif_dict = {}
            # else:
            #     if 'thalamus_neurons' in edge_source_pop:
            #         # --- Identifies if the synapse SOURCE is VPL_TC
            #         if   edge_source_pop.split('|')[1]=='VPL_TC':   
            #             # modif_dict = {}
            #             modif_dict = {'mg':1}
            #         elif (edge_source_pop.split('|')[1]=='VPL_IN'):
            #             # if edge_target_pop == 'VPL_TC':             modif_dict = {'e_GABAA':-94,'e_GABAB':-97}
            #             if edge_target_pop == 'VPL_TC':             modif_dict = {'e_GABAA':-94,'e_GABAB':-97, 'tau_d_GABAB':77} # 2023-09-05: testing if midification in 'tau_d_GABAB':77 should apply to all synapses, or just to the ones that are not Rt_RC and VPL_IN (if any exist)
            #         elif (edge_source_pop.split('|')[1]=='Rt_RC'):
            #             # if edge_target_pop == 'VPL_TC':             modif_dict = {'e_GABAA':-94,'e_GABAB':-97}
            #             if edge_target_pop == 'VPL_TC':             modif_dict = {'e_GABAA':-94,'e_GABAB':-97, 'tau_d_GABAB':77}
            #             else:                                       modif_dict = {'e_GABAA':-82,'e_GABAB':-97, 'tau_d_GABAB':77}
            #         else:
            #             modif_dict={}
            #             print('no matching syn')
                
        return modif_dict
    
    def createConns(node_pathway_edges,thal_gid,cell_properties,stimType='VecStim',numCells=1,mod_type='Prob'):
        mechs_dict      = {}
        preCells_dict   = {}
        prePops_dict    = {}
        conns_dict      = {}

        if   mod_type == 'Det':                 modAMPANMDA = 'DetAMPANMDA';                    modGABA     = 'DetGABAAB'
        elif mod_type == 'Prob_S1':             modAMPANMDA = 'ProbAMPANMDA_EMS_S1';            modGABA     = 'ProbGABAAB_EMS_S1'
        elif mod_type == 'Prob_original':       modAMPANMDA = 'ProbAMPANMDA_EMS_original';         modGABA     = 'ProbGABAAB_EMS_original' # original MOD from BBP model
        elif mod_type == 'Prob_CORENEURON':     modAMPANMDA = 'ProbAMPANMDA_EMS_CORENEURON';    modGABA     = 'ProbGABAAB_EMS_CORENEURON'
        elif mod_type == 'Prob_S1_CORENEURON':  modAMPANMDA = 'ProbAMPANMDA_EMS_S1_CORENEURON'; modGABA     = 'ProbGABAAB_EMS_S1_CORENEURON'
        else:                                   modAMPANMDA = 'ProbAMPANMDA_EMS';               modGABA     = 'ProbGABAAB_EMS'

        for edge_name in node_pathway_edges.keys():
            mechs_dict.update(      {edge_name:{}})
            preCells_dict.update(   {edge_name:{}})
            prePops_dict.update(    {edge_name:{}})
            conns_dict.update(      {edge_name:{}})
            print('Processing connections: ', edge_name)
            edge_type,edge_mech,edge_source,edge_target = edge_name.split('__')
            
            if edge_mech == 'chemical':
                if      'Rt_RC'  in edge_source:    mod_file=modGABA;       mod_label = 'Inh'
                elif    'VPL_IN' in edge_source:    mod_file=modGABA;       mod_label = 'Inh'
                else:                               mod_file=modAMPANMDA;   mod_label = 'Exc'
            elif edge_mech == 'electrical_synapse': mod_file='ElectSyn';    mod_label = 'Gap'
            # elif edge_mech == 'electrical_synapse': mod_file='Gap';         mod_label = 'Gap'
            
            # --- Processing edges to create conns
            for edge_index in node_pathway_edges[edge_name].index:
                edge_data = node_pathway_edges[edge_name].loc[edge_index]
                TableS1 = StoreParameters.getTableS1()
                syn_values = TableS1[edge_source][cell_properties['mtype']]
                # --- Updates the Edge dataset to include updated NMDA_ratio values, which the authors set manually afterwards during the model run
                if syn_values['paper_reference_values']['NMDA_ratio'] is not None:
                    nmda_ratio_df = pd.Series([syn_values['paper_reference_values']['NMDA_ratio']], index=['NMDA_ratio'])
                    edge_data = edge_data.append(nmda_ratio_df)
                
                edge_dict=CreateNetPyNE.modTemplate(edge_data,mod_file)


                mechs_dict[edge_name].update({mod_label+'__'+str(edge_index):edge_dict})
                # --- VecStim artificial spike generator
                preCells_dict[edge_name].update({stimType+'_cell__'+str(edge_index):{   'cellModel': 'VecStim', 'numCells': numCells}})
                prePops_dict[edge_name].update({stimType+'_pop__'+str(edge_index): {    'cellType': stimType+'_cell__'+str(edge_index), 
                                                                                        'numCells': 1, 
                                                                                        'interval': 50, 
                                                                                        'noise': 0, 
                                                                                        'start': 8000,
                                                                                        # 'pulses': [ {'start': 500,  'end': 501,  'rate': 10, 'noise': 0}, 
                                                                                        #             {'start': 2500, 'end': 2901, 'rate': 40, 'noise': 0},
                                                                                        #             {'start': 3400, 'end': 3401, 'rate': 10, 'noise': 0},
                                                                                        #         ]
                                                                                        }})
                conns_dict[edge_name].update({stimType+'_conn__'+str(edge_index):{      'preConds':     {'pop': stimType+'_pop__'+str(edge_index)}, 
                                                                                        'postConds':    {'pop': 'pop__'+str(thal_gid)},
                                                                                        'probability':  1,
                                                                                        'weight':       1,
                                                                                        'synsPerConn':  1,
                                                                                        'synMech':      mod_label+'__'+str(edge_index),
                                                                                        'delay':        edge_data['delay'],
                                                                                        'sec':          node_pathway_edges[edge_name].loc[edge_index]['sec'], # from <NetPyNE_BBP.ConvertSynapses.getTargetSections>
                                                                                        'loc':          node_pathway_edges[edge_name].loc[edge_index]['afferent_section_pos'], # from circuit edges,
                                                                                        }})
    
        return mechs_dict, preCells_dict, prePops_dict, conns_dict

##################################################################################################################################################################################################################      

class ConvertMorphology():
    # def convert_netpyne():
    #     from netpyne import specs, sim
    #     import netpyne
    #     # Network parameters
    #     netParams = specs.NetParams()  # object of class NetParams to store the network parameters
    #     netParams.importCellParams(
    #         label=df_postCell['morphology'], 
    #         fileName=df_postCell['etype']+'.hoc', 
    #         cellName=df_postCell['etype'], 
    #         conds={}, 
    #         cellArgs=[gid,morphology_sourceFolder,df_postCell['morphology']+'.swc'], 
    #         importSynMechs=True, 
    #         somaAtOrigin=False, 
    #         cellInstance=False,
    #     )


    def convert_cell(netParams,gid,df_postCell,morphology_sourceFolder,morphology_outputFolder=None):
        netParams.importCellParams( 
            label=df_postCell['morphology'], 
            fileName=df_postCell['etype']+'.hoc', 
            cellName=df_postCell['etype'], 
            conds={}, 
            cellArgs=[gid,morphology_sourceFolder,df_postCell['morphology']+'.swc'], 
            importSynMechs=True, 
            somaAtOrigin=False, 
            cellInstance=False,
        )

##################################################################################################################################################################################################################

class StoreParameters():
    def getTableS1(suppress_output=True):
        if not suppress_output:print('Loading data from TableS1')
        TableS1={  'thalamus_neurons|Rt_RC':{
                        'Rt_RC':    {'syn_type':'I2',   'mod_file':'ProbGABAAB_EMS',     'paper_reference_values':{'conductance':0.9, 'conductance_std':0.23,   'decay_time':8.3, 'decay_time_std':2.2,   'NMDA_ratio':None, 'u_syn':0.41, 'u_syn_std':0.14, 'depression_time':464, 'depression_time_std':339, 'facilitation_time':54,  'facilitation_time_std':71, }},
                        'VPL_TC':   {'syn_type':'I2',   'mod_file':'ProbGABAAB_EMS',     'paper_reference_values':{'conductance':1.1, 'conductance_std':0.4,    'decay_time':8.3, 'decay_time_std':2.2,   'NMDA_ratio':None, 'u_syn':0.32, 'u_syn_std':0.18, 'depression_time':352, 'depression_time_std':46,  'facilitation_time':2,   'facilitation_time_std':209, }},
                        'VPL_IN':   {'syn_type':'I2',   'mod_file':'ProbGABAAB_EMS',     'paper_reference_values':{'conductance':0.9, 'conductance_std':0.23,   'decay_time':8.3, 'decay_time_std':2.2,   'NMDA_ratio':None, 'u_syn':0.41, 'u_syn_std':0.14, 'depression_time':464, 'depression_time_std':339, 'facilitation_time':54,  'facilitation_time_std':71, }},
                        },
                    'thalamus_neurons|VPL_TC':{
                        'Rt_RC':    {'syn_type':'E2',   'mod_file':'ProbAMPANMDA_EMS',   'paper_reference_values':{'conductance':2.8, 'conductance_std':0.1,    'decay_time':1.58, 'decay_time_std':0.26, 'NMDA_ratio':0.57, 'u_syn':0.86, 'u_syn_std':0.09, 'depression_time':671, 'depression_time_std':17,  'facilitation_time':17,  'facilitation_time_std':5, }},
                        },
                    'thalamus_neurons|VPL_IN':{
                        'VPL_TC':   {'syn_type':'I2',   'mod_file':'ProbGABAAB_EMS',     'paper_reference_values':{'conductance':0.4, 'conductance_std':0.4,    'decay_time':8.3, 'decay_time_std':2.2,   'NMDA_ratio':None, 'u_syn':0.47, 'u_syn_std':0.18, 'depression_time':137, 'depression_time_std':46,  'facilitation_time':239, 'facilitation_time_std':209, }},
                        'VPL_IN':   {'syn_type':'I2',   'mod_file':'ProbGABAAB_EMS',     'paper_reference_values':{'conductance':2.7, 'conductance_std':0.4,    'decay_time':8.3, 'decay_time_std':2.2,   'NMDA_ratio':None, 'u_syn':0.41, 'u_syn_std':0.14, 'depression_time':464, 'depression_time_std':339, 'facilitation_time':54,  'facilitation_time_std':71, }},
                        },
                    'MedialLemniscus_projections':{
                        'VPL_TC':   {'syn_type':'E2',   'mod_file':'ProbAMPANMDA_EMS',   'paper_reference_values':{'conductance':1.15,'conductance_std':0.2,    'decay_time':1.74, 'decay_time_std':0.18, 'NMDA_ratio':0.41, 'u_syn':0.3,  'u_syn_std':0.21, 'depression_time':2350,'depression_time_std':315, 'facilitation_time':1,   'facilitation_time_std':2, }},
                        'VPL_IN':   {'syn_type':'E2',   'mod_file':'ProbAMPANMDA_EMS',   'paper_reference_values':{'conductance':1.15,'conductance_std':0.2,    'decay_time':1.74, 'decay_time_std':0.18, 'NMDA_ratio':0.41, 'u_syn':0.48, 'u_syn_std':0.21, 'depression_time':690, 'depression_time_std':315, 'facilitation_time':57,  'facilitation_time_std':53, }},
                        },
                    'CorticoThalamic_projections':{
                        'Rt_RC':    {'syn_type':'E1',   'mod_file':'ProbAMPANMDA_EMS',   'paper_reference_values':{'conductance':0.16, 'conductance_std':0.016, 'decay_time':2.74, 'decay_time_std':0.25, 'NMDA_ratio':0.99, 'u_syn':0.09, 'u_syn_std':0.12, 'depression_time':138, 'depression_time_std':211, 'facilitation_time':670, 'facilitation_time_std':830, }},
                        'VPL_TC':   {'syn_type':'E1',   'mod_file':'ProbAMPANMDA_EMS',   'paper_reference_values':{'conductance':0.16, 'conductance_std':0.016, 'decay_time':1.74, 'decay_time_std':0.18, 'NMDA_ratio':1.91, 'u_syn':0.09, 'u_syn_std':0.12, 'depression_time':138, 'depression_time_std':211, 'facilitation_time':670, 'facilitation_time_std':830, }},
                        'VPL_IN':   {'syn_type':'E1',   'mod_file':'ProbAMPANMDA_EMS',   'paper_reference_values':{'conductance':0.16, 'conductance_std':0.016, 'decay_time':1.74, 'decay_time_std':0.18, 'NMDA_ratio':0.99, 'u_syn':0.09, 'u_syn_std':0.12, 'depression_time':138, 'depression_time_std':211, 'facilitation_time':670, 'facilitation_time_std':830, }},
                        },
                    }
        return TableS1
    
    def getTableS1_barreloidThalamus(suppress_output=True):
        if not suppress_output:print('Loading data from TableS1 - Formatted to the barreloid thalamus pops')
        TableS1={  'TRN':{
                        'TRN':      {'syn_type':'I2',   'mod_file':'ProbGABAAB_EMS',     'paper_reference_values':{'conductance':0.9, 'conductance_std':0.23,   'decay_time':8.3, 'decay_time_std':2.2,   'NMDA_ratio':None, 'u_syn':0.41, 'u_syn_std':0.14, 'depression_time':464, 'depression_time_std':339, 'facilitation_time':54,  'facilitation_time_std':71, }},
                        'VPM':      {'syn_type':'I2',   'mod_file':'ProbGABAAB_EMS',     'paper_reference_values':{'conductance':1.1, 'conductance_std':0.4,    'decay_time':8.3, 'decay_time_std':2.2,   'NMDA_ratio':None, 'u_syn':0.32, 'u_syn_std':0.18, 'depression_time':352, 'depression_time_std':46,  'facilitation_time':2,   'facilitation_time_std':209, }},
                        },
                    'VPM':{
                        'TRN':      {'syn_type':'E2',   'mod_file':'ProbAMPANMDA_EMS',   'paper_reference_values':{'conductance':2.8, 'conductance_std':0.1,    'decay_time':1.58, 'decay_time_std':0.26, 'NMDA_ratio':0.57, 'u_syn':0.86, 'u_syn_std':0.09, 'depression_time':671, 'depression_time_std':17,  'facilitation_time':17,  'facilitation_time_std':5, }},
                        'L6A':      {'syn_type':'E2',   'mod_file':'ProbAMPANMDA_EMS',   'paper_reference_values':{'conductance':0.72,'conductance_std':0.5,    'decay_time':1.74, 'decay_time_std':0.18, 'NMDA_ratio':0.71, 'u_syn':0.72, 'u_syn_std':0.065, 'depression_time':227, 'depression_time_std':38,  'facilitation_time':14,  'facilitation_time_std':12, }},
                        },
                    'MLe':{
                        'VPM':      {'syn_type':'E2',   'mod_file':'ProbAMPANMDA_EMS',   'paper_reference_values':{'conductance':1.15,'conductance_std':0.2,    'decay_time':1.74, 'decay_time_std':0.18, 'NMDA_ratio':0.41, 'u_syn':0.3,  'u_syn_std':0.21, 'depression_time':2350,'depression_time_std':315, 'facilitation_time':1,   'facilitation_time_std':2, }},
                        },
                    'L6A':{
                        'TRN':      {'syn_type':'E1',   'mod_file':'ProbAMPANMDA_EMS',   'paper_reference_values':{'conductance':0.16, 'conductance_std':0.016, 'decay_time':2.74, 'decay_time_std':0.25, 'NMDA_ratio':0.99, 'u_syn':0.09, 'u_syn_std':0.12, 'depression_time':138, 'depression_time_std':211, 'facilitation_time':670, 'facilitation_time_std':830, }},
                        'VPM':      {'syn_type':'E1',   'mod_file':'ProbAMPANMDA_EMS',   'paper_reference_values':{'conductance':0.16, 'conductance_std':0.016, 'decay_time':1.74, 'decay_time_std':0.18, 'NMDA_ratio':1.91, 'u_syn':0.09, 'u_syn_std':0.12, 'depression_time':138, 'depression_time_std':211, 'facilitation_time':670, 'facilitation_time_std':830, }},
                        },
                    }
        return TableS1

##################################################################################################################################################################################################################

class LoadSimResults():
    def getSpikesFromConnectedCells(th_spk_dict, thal_gids, edge_sources_dict, showFig=False):
        
        # --- Allocates gids into a single list
        edge_sources_dict_gids=[]
        for pop in edge_sources_dict.keys(): 
            if 'thalamus_neurons|' in pop:
                edge_sources_dict_gids+=edge_sources_dict[pop]

        # --- Adds presyn cell gids
        edge_sources_dict_gids+=thal_gids
        valid_gids = list(set(edge_sources_dict_gids))
        valid_gids.sort()

        # --- Adds gid and spike times into the new dictionary
        th_connected_cells_dict={}
        for pop in th_spk_dict.keys():
            print(pop,len(th_spk_dict[pop].keys()))
            th_connected_cells_dict.update({pop:{}})
            for gid in th_spk_dict[pop]:
                if gid in valid_gids:
                    th_connected_cells_dict[pop].update({gid:th_spk_dict[pop][gid]})

        # --- Comparing the unfiltered vs filtered dictionary
        for pop in th_spk_dict.keys():
            spiking_cells=[];non_spiking_cells=[]
            print(pop, 'all  cells: ', len(th_spk_dict[pop]))
            print(pop, 'conn cells: ', len(th_connected_cells_dict[pop]))
            for cell_gid in th_connected_cells_dict[pop]:
                if len(th_connected_cells_dict[pop][cell_gid])>0:   spiking_cells.append(cell_gid)
                else:                                               non_spiking_cells.append(cell_gid)

            print(pop, '\tspiking cells: ', len(spiking_cells), '\tnon-spiking cells: ', len(non_spiking_cells))

        return th_connected_cells_dict

        
        # f       = h5py.File(filePath, 'r')
        # dset    = f['spikes']
        # spkids  = list(dset['All']['node_ids'])
        # spkts   = list(dset['All']['timestamps'])
        # spkids_set = list(set(spkids)); spkids_set.sort()
        # print('# of individual spiking cells: ',len(spkids_set))
        # spk_dict={}                                                         # --- Creates a dictionary with spike times by cell GID 
        # for spk in spkids_set:spk_dict.update({spk:[]})                     # --- Allocates the GIDs as the keys 
        # for ind,spkt in enumerate(spkts):spk_dict[spkids[ind]].append(spkt) # --- Assigns the spike times to each matching GID
        
        # connected_cells={}
        # spiking_cells=[];non_spiking_cells=[]
        # for pop in edge_sources_dict.keys():
        #     if 'thalamus_neurons|' in pop:  pop_name = pop.split('|')[1]
        #     else:                           pop_name = pop
        #     connected_cells.update({pop_name:{}})
        #     for cell_gid in edge_sources_dict[pop]:

        #         if cell_gid in spk_dict.keys():
        #             cell_spks = spk_dict[cell_gid]
        #             spiking_cells.append(cell_gid)
        #             # print(spk_dict[cell_gid])
        #         else:
        #             cell_spks =[]
        #             non_spiking_cells.append(cell_gid)
        #             # print('cell # ', cell_gid, ' not spiking')

        #         connected_cells[pop_name].update({cell_gid:cell_spks})

        # print('\nspiking_cells: ', spiking_cells,'\n\nnon_spiking_cells: ',non_spiking_cells)

        # return connected_cells

    def getSpikesFromSim(filePath, cfg_file, microcircuit_number=None, showFig=False):
        f       = h5py.File(filePath, 'r')
        dset    = f['spikes']
        spkids  = list(dset['All']['node_ids'])
        spkts   = list(dset['All']['timestamps'])
        spkids_set = list(set(spkids)); spkids_set.sort()
                
        spk_dict={}                                                         # --- Creates a dictionary with spike times by cell GID 
        for spk in spkids_set:spk_dict.update({spk:[]})                     # --- Allocates the GIDs as the keys 
        for ind,spkt in enumerate(spkts):spk_dict[spkids[ind]].append(spkt) # --- Assigns the spike times to each matching GID
        
        #### FIX THIS CODE TO TAKE THE CONNECTED NODES INSTEAD OF ALL NODES
        circuit_dict    = LoadBBPCircuit.getDataFrames_new(cfg_file=cfg_file)

        th_spk_dict={}
        for th_pop in list(set(circuit_dict['thalamus_neurons']['mtype'])):
            if microcircuit_number is None: gids = list(circuit_dict['thalamus_neurons'].loc[(circuit_dict['thalamus_neurons']['mtype'] == th_pop)].index)
            else:                           gids = list(circuit_dict['mc'+str(microcircuit_number)][th_pop].index)
            gids.sort()
            store_gid_spkts={}
            # missing_gids=[]
            for gid in gids:
                store_gid_spkts.update({gid:[]})
                # if gid in spkids_set:
                #     store_gid_spkts.update({gid:[]})
                # else:
                #     missing_gids.append(gid)
            th_spk_dict.update({th_pop:store_gid_spkts})
            # print(th_pop,'missing gids \n',missing_gids)

        for th_pop in list(set(circuit_dict['thalamus_neurons']['mtype'])):
            for cell_gid in spk_dict.keys():
                if cell_gid in th_spk_dict[th_pop].keys(): th_spk_dict[th_pop][cell_gid]=spk_dict[cell_gid]

        #### 
        if showFig:
            import matplotlib.pyplot as plt
            plt.figure(figsize=(20,15))

            th_pop_colors=['royalblue','r','g']
            spkcs=[]
            for spkid_ind,spkid in enumerate(spkids):
                for th_pop_ind,th_pop in enumerate(list(set(circuit_dict['thalamus_neurons']['mtype']))):
                    if spkid in th_spk_dict[th_pop].keys(): spkcs.append(th_pop_colors[th_pop_ind])
            plt.scatter(spkts,spkids,marker='.',s=2,color=spkcs)
            plt.ylabel('Cell Gids')
            plt.xlabel('Time (ms)')
            plt.title('Raster plot of BBP model data')
            plt.show()
            plt.savefig('../data/bbp_model_raster.png',dpi=500)

        return th_spk_dict
    
    def getVirtualInputSpikes(filePath, filterGids = None, skipFirstLine=True,setIndexToZero=True,showFig=False,saveName=None,loadFromSourceFile=True):
        if not loadFromSourceFile:
            try:
                if '.json' not in saveName: saveName+='.json'
                with open(saveName, 'r') as openfile: noise_input_dict = json.load(openfile)                
                filterGids_str=[str(gid_) for gid_ in filterGids]
                # --- formatting the lists before comparing
                filterGids_str.sort()
                noise_input_dict_keys=list(noise_input_dict.keys()); noise_input_dict_keys.sort()
                if filterGids_str!=noise_input_dict_keys: loadFromSourceFile=True # Loaded GIDs dont match. Gerating new GIDs dictionary
                else: print('\t\t- \tLoading noise from saved spike file: ', saveName)
            except: loadFromSourceFile=True # Converted noise file doesnt exist

        if loadFromSourceFile:
            print('\t\t- \tLoading noise from source file')
            f = open(filePath, "r")
            if skipFirstLine:   noise_spikes=f.readlines()[1:]
            else:               noise_spikes=f.readlines()
            
            # --- Formatting the data into [spkid,spkt] pairs
            noise_spikes_1=[]
            for line in noise_spikes:
                line_=line.strip('\n'); 
                line__=line_.split('\t')
                noise_spikes_1.append([float(line__[0]),int(line__[1])])

            edge_gids=[]
            for noise_spike in noise_spikes_1:edge_gids.append(noise_spike[1])
            edge_gids_set=list(set(edge_gids)); edge_gids_set.sort()
            min_set_val=min(edge_gids_set)
            if setIndexToZero: edge_gids_set=[gid-min_set_val for gid in edge_gids_set]

            noise_input_dict={}
            for edge_gid_s    in edge_gids_set:  noise_input_dict.update({edge_gid_s:[]})
            for spkt,edge_gid in noise_spikes_1: noise_input_dict[edge_gid-min_set_val].append(spkt)

            if filterGids is not None:
                noise_input_dict_={}
                for key in noise_input_dict.keys():
                    if key in filterGids: noise_input_dict_.update({key:noise_input_dict[key]})
                del noise_input_dict
                noise_input_dict=noise_input_dict_
        if saveName is not None:
            if '.json' not in saveName: saveName+='.json'
            # Writing to json
            json_object = json.dumps(noise_input_dict, indent=4)
            with open(saveName, "w") as outfile: outfile.write(json_object)

        return noise_input_dict

    def getReportIndexes(circuit_dict, thal_gids):
        ####################################################################################
        # The soma report file returns a matrix with the GIDs for the microcircuit 'mc2'
        # This function converts the original GIDs from 'mc2' into indexes from range [0,len('mc2')]
        ####################################################################################
        
        # --- Creates a dictionary to store {key:value} pairs for {bbp_model_gid:bbp_report_gid}
        ind=0
        mc2_report_index_dict={}
        for key in ['Rt_RC','VPL_TC','VPL_IN']: # this order must not be changed!!!! gids are assigned in this order (Rt_RC: 4909 - VPL_TC: 8952 - VPL_IN: 47)
            for index in list(circuit_dict['mc2'][key].index):
                mc2_report_index_dict.update({index:ind})
                ind+=1
        # --- Selects the bbp_report_gid values based on the bbp_model_gid values in order to select voltage traces
        thal_gid_report_indexes_dict={}
        for thal_gid in thal_gids: 
            thal_gid_report_indexes_dict.update({thal_gid:mc2_report_index_dict[thal_gid]})
        return thal_gid_report_indexes_dict # --- Returns a dictionary with {key:value} pairs of {bbp_model_gid:bbp_report_gid}

    def loadVoltageReport(filePath, cfg_file, gids, timeWindow, dt, microcircuit_number=None, showFig=False):
        circuit_dict = LoadBBPCircuit.getDataFrames_new( cfg_file=cfg_file)
        f = h5py.File(filePath, 'r')
        if type(gids)==int:    selected_cells = [gids]
        elif type(gids)==list: selected_cells = gids
        else:                  return
        if (type(timeWindow)==list) and len(timeWindow)==2: tWindow = [timeWindow[0]/dt,timeWindow[1]/dt]
        thal_gid_report_indexes_dict = LoadSimResults.getReportIndexes(circuit_dict, gids)
        time_trace = np.arange(timeWindow[0],timeWindow[1],dt)
        th_soma_voltage_dict={'t':time_trace,'voltage':{}}
        for gid in thal_gid_report_indexes_dict.keys():
            report_index  = thal_gid_report_indexes_dict[gid]
            voltage_trace = f['report']['VPL_TC']['data'][:,report_index]
            
            th_soma_voltage_dict['voltage'].update({gid:voltage_trace})
            
            #### 
            if showFig:
                import matplotlib.pyplot as plt
                plt.figure(figsize=(20,15))
                plt.plot(time_trace,voltage_trace,'k')
                # plt.plot([5000,5000],[-90,50],'-.r')
                plt.ylim(-90,50)
                plt.ylabel('Membrane voltage (mV)')
                plt.xlabel('Time (ms)')
                # plt.title('Raster plot of BBP model data')
                plt.show()
                plt.savefig('../data/bbp_model_somaVoltages_'+str(gid)+'.png',dpi=500)
            
        return th_soma_voltage_dict