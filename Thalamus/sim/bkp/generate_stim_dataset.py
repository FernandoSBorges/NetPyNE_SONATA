#! /Users/joao/opt/anaconda3/bin/ipython
'''
Class to build the stimulation dataset
Contributors: joao.moreira@downstate.edu, salvadordura@gmail.com
'''
########################################################################################################################
import numpy as np
import json
import math
import pandas as pd
########################################################################################################################
from GenerateStimDataset import LoadPSTH
from GenerateStimDataset import PlotFigs
from GenerateStimDataset import SampleData
from GenerateStimDataset import FitData
########################################################################################################################
NetPyNE_rootFolder              = '/Users/joao/Research/Models/BBP/thalamus_netpyne'
stim_source_code_folder         = '/stims/Minnery_2003/fig1'
stim_dataset_folder             = NetPyNE_rootFolder+stim_source_code_folder+'/stim_datasets/'      # Folder with the digitized PSTH
save_figs_folder                = NetPyNE_rootFolder+stim_source_code_folder+'/figs_stim/'          # Folder to save the output figures
save_deflection_model_folder    = NetPyNE_rootFolder+stim_source_code_folder+'/deflection_model/'   # Folder to save the output deflection model and sampled raster
########################################################################################################################
UpdateDeflectionModel       = False

if UpdateDeflectionModel: 
    SampleRaster            = True
    UpdateSampledRasterData = True
else:
    SampleRaster = True
    UpdateSampledRasterData = False

UpdateSampledRasterFigs  = True
saveSimplifiedRasterData = True # Saves an aditional dictionary dataset with only the spike times organized by GID - lighter, to be synchronized in Github
save_format              = 'pkl'
########################################################################################################################
# deflection_events = [[500,650],[1000,1250],[1500,1850]]
deflection_events = [[2000,2150],[3000,3150],[4000,4250],[5000,5250],[6000,6350],[7000,7350]]
########################################################################################################################
# --- Plot figures
plotFigs = True
########################################################################################################################
# --- scaling probability to adjust for a higher re-sampling rate
rescale_sampling = False            # (True: 0.025 ms) / (False: 1.0 ms)

dt_step_PSTH    = 1.0    # (ms) - extracted from paper methods (measured at 0.1 ms resolution, and binned at 1.0 ms bins)
dt_step_NetPyNE = 0.025  # (ms) - from netpyne simulation cfg file

if rescale_sampling: dt_step_resampling = dt_step_NetPyNE
else:                dt_step_resampling = dt_step_PSTH

########################################################################################################################
# --- Angle to sample rasters from
target_angles = [0,45,90,135,180,225,270,315,None]
# target_angles = [None]
# target_angles = [90,None]
# target_angles = [0,45,135,180,225,270,315]
# target_angles = [0,45,90,135,180,225,270,315]
########################################################################################################################

if UpdateDeflectionModel:
    # Specify the path to your CSV file
    folderPath = stim_dataset_folder
    fileNames = ['PrV.csv','VPM.csv','0stim.csv']
    hist_dict = LoadPSTH.load(folderPath, fileNames)
    if plotFigs: PlotFigs.plotAllTraces(hist_dict, save_path = save_figs_folder)
    ########################################################################################################################
    # --- Skips the samples for the rise time, until I find a better model to fit the curve (poisson, maybe?)
    target_nuclei = 'PrV'

    ON_start  = 62
    OFF_start = 202

    skip_ON  = 3 # samples
    skip_OFF = 2 # samples

    # --- Selecting and formating the data from PrV
    PrV_data = LoadPSTH.formatData(hist_dict,target_nuclei='PrV', ON_start=ON_start, OFF_start=OFF_start, skip_ON=skip_ON, skip_OFF=skip_OFF)

    # --- Calculating the x_data shifted to x=0
    x_data_ON_zero  = [x-min(PrV_data['ON']['x'])  for x in PrV_data['ON']['x']]
    x_data_OFF_zero = [x-min(PrV_data['OFF']['x']) for x in PrV_data['OFF']['x']]

    # --- Updating the dataset with a zero aligned data
    PrV_data.update({'ON_zero' :{'x':x_data_ON_zero, 'y':PrV_data['ON']['y']}})
    PrV_data.update({'OFF_zero':{'x':x_data_OFF_zero,'y':PrV_data['OFF']['y']}})

    #######################################################################################################################################
    # --- Fitting the baseline
    baseline_mean = np.mean(PrV_data['baseline']['y'])
    baseline_std  = np.std( PrV_data['baseline']['y'])

    fitDegree = 25

    # --- DATA ON ORIGINAL TIMESTAMPS
    # --- Fits polynomial to data
    yPoly_ON       = FitData.fitPolynomial(x_data=PrV_data['ON']['x'],       y_data=PrV_data['ON']['y'], degree=fitDegree)
    yPoly_OFF      = FitData.fitPolynomial(x_data=PrV_data['OFF']['x'],      y_data=PrV_data['OFF']['y'], degree=fitDegree)
    # --- Creates polynomial Object to sample from
    yFit_ON       = FitData.getPolynomialObj(yPoly_ON)
    yFit_OFF      = FitData.getPolynomialObj(yPoly_OFF)

    # --- DATA ALIGNED TO ZERO (x0=0)
    # --- Fits polynomial to data
    yPoly_ON_zero  = FitData.fitPolynomial(x_data=PrV_data['ON_zero']['x'],  y_data=PrV_data['ON']['y'], degree=fitDegree)
    yPoly_OFF_zero = FitData.fitPolynomial(x_data=PrV_data['OFF_zero']['x'], y_data=PrV_data['OFF']['y'], degree=fitDegree)
    # --- Creates polynomial Object to sample from
    yFit_ON_zero   = FitData.getPolynomialObj(yPoly_ON_zero)
    yFit_OFF_zero  = FitData.getPolynomialObj(yPoly_OFF_zero)

    if plotFigs: PlotFigs.plotFits(PrV_data,yFit_ON,yFit_OFF,yFit_ON_zero,yFit_OFF_zero,nSamples=250, save_path = save_figs_folder)

    #######################################################################################################################################
    # --- Tunning function - Extracted from Fig 3A
    angular_tuning = {  'x':[  0,      1,                  2,                  3,                   4],
                        'y':[  1.0,    0.6982148208637337, 0.5676691729323308, 0.44418423989406364, 0.4282115869017633]}

    yPoly_angTuning = FitData.fitPolynomial(x_data=angular_tuning['x'], y_data=angular_tuning['y'], degree=3)
    yFit_angTuning  = FitData.getPolynomialObj(yPoly_angTuning)

    if plotFigs: PlotFigs.plotAngularTuning(angular_tuning,nSamples=250,yFit_angTuning=yFit_angTuning, save_path = save_figs_folder)

    #######################################################################################################################################
    # --- Make a histogram of spikes

    # --- Dictionary to store the Sampled Data
    sampledProb={}

    # --- Sampling probabilities from the baseline
    baseline_duration = math.floor(max(PrV_data['baseline']['x']))
    sampledProb.update({'baseline':SampleData.sampleNormal(mean=baseline_mean,std=baseline_std,start=0,stop=baseline_duration,dt=dt_step_resampling)})

    # --- Sampling probabilities from the polynomial fit
    sampledProb.update({'ON' :SampleData.sampleFitFunction(yFit_ON , PrV_data['ON']['x'] , dt_step_resampling)})
    sampledProb.update({'OFF':SampleData.sampleFitFunction(yFit_OFF, PrV_data['OFF']['x'], dt_step_resampling)})

    # --- Sampling probabilities from the polynomial fit - dataset with time set to ZERO
    sampledProb.update({'ON_zero' :SampleData.sampleFitFunction(yFit_ON_zero , PrV_data['ON_zero']['x'] , dt_step_resampling)})
    sampledProb.update({'OFF_zero':SampleData.sampleFitFunction(yFit_OFF_zero, PrV_data['OFF_zero']['x'], dt_step_resampling)})

    # --- Plotting figs
    if plotFigs: PlotFigs.plotSampledDataBar(sampledProb,save_path = save_figs_folder) # Fig 8
    if plotFigs: PlotFigs.plotFitAndSampled(PrV_data, yFit_ON, yFit_OFF, nSamples=250, sampledData=sampledProb, save_path = save_figs_folder) # Fig 81

    #######################################################################################################################################
    # --- Generating whisker deflection dictionary
    for target_angle in target_angles:
        if target_angle is not None:    deflection_times = deflection_events+[[10000,10020]]
        else:                           deflection_times = [[10000,10020]]

        print('Updating whisker deflection Model')
        if target_angle is not None:
            deflection_dict={   
                                # 'deflection_times':             [[500,650],[1000,1250],[1500,1850],[10000,10020]],
                                'deflection_times':             deflection_times,
                                'baseline_prob':                [baseline_mean,baseline_std],
                                'ON_polyFit':                   list(yPoly_ON_zero),
                                'OFF_polyFit':                  list(yPoly_OFF_zero),
                                'resetTime':                    True,
                                'dt':                           dt_step_resampling,
                                'plotFigs':                     True,
                                'target_angle':                 target_angle,
                                'n_cells':                      180,
                                'prob_compensation':            1.25,
                                'yPoly_angTuning':              list(yPoly_angTuning),
                                'maximize_target_angle':        True,
                                'save_stim_data':               True,
                                'analyze_stim_data':            True,
                                'paper_sampling':               dt_step_PSTH,
                                'addMouseHead':                 True,
                                'stim_dataset_folder':          stim_dataset_folder,
                                'save_figs_folder':             save_figs_folder,
                                'save_deflection_model_folder': save_deflection_model_folder,
                                }
            stims_string ='_stims|' + "|".join(str(sublist[0])+'@'+str(sublist[1]) for sublist in deflection_times)
        else:
            # --- Generates a baseline stimulation dataset
            deflection_dict={   
                                'deflection_times':             deflection_times,
                                'baseline_prob':                [baseline_mean,baseline_std],
                                'ON_polyFit':                   list(yPoly_ON_zero),
                                'OFF_polyFit':                  list(yPoly_OFF_zero),
                                'resetTime':                    True,
                                'dt':                           dt_step_resampling,
                                'plotFigs':                     True,
                                'target_angle':                 target_angle,
                                'n_cells':                      180,
                                'prob_compensation':            1,
                                'yPoly_angTuning':              None,
                                'maximize_target_angle':        False,
                                'save_stim_data':               True,
                                'analyze_stim_data':            True,
                                'paper_sampling':               dt_step_PSTH,
                                'addMouseHead':                 True,
                                'stim_dataset_folder':          stim_dataset_folder,
                                'save_figs_folder':             save_figs_folder,
                                'save_deflection_model_folder': save_deflection_model_folder,
                                }
            stims_string = '_stims|None'

        # --- Saves the stim_string into the deflection dictionary for retrieval
        deflection_dict.update({'stims_string':stims_string})

        # --- Saves the sampling interval (dt) into the deflection dictionary for retrieval
        dt_string = '_samplingBins|'+str(int(dt_step_resampling*1000))+'us'

        # --- Save the deflection dictionary
        with open(save_deflection_model_folder+'deflectionModel_'+str(deflection_dict['target_angle'])+'|deg'+dt_string+stims_string+'.json', 'w') as fp: json.dump(deflection_dict, fp, indent=4)

if SampleRaster:
    #######################################################################################################################################
    # --- Sample whisker deflection from deflection model
    for target_angle in target_angles:
        if target_angle is not None:    
            deflection_times = deflection_events+[[10000,10020]]
            stims_string ='_stims|' + "|".join(str(sublist[0])+'@'+str(sublist[1]) for sublist in deflection_times)
        else:                           
            deflection_times = [[10000,10020]]
            stims_string = '_stims|None'
        
        # --- Saves the sampling interval (dt) into the deflection dictionary for retrieval
        dt_string = '_samplingBins|'+str(int(dt_step_resampling*1000))+'us'


        deflection_dict = SampleData.LoadJSON(file_path=save_deflection_model_folder+'deflectionModel_'+str(target_angle)+'|deg'+dt_string+stims_string+'.json')
        # deflection_dict = SampleData.LoadJSON(file_path=save_deflection_model_folder+'deflectionModel_'+str(target_angle)+'|deg.json')
        print(' \n\n\n##################################################################')
        print(' ---> Target Angle: ', target_angle)
        

        # --- Make a new raster
        if UpdateSampledRasterData:
            print('Generating dataset')
            spikes_dict = SampleData.deflectionEvent(deflection_dict)
            if save_format=='json':
                with open(deflection_dict['save_deflection_model_folder']+'spike_dicts/mleSpikes_deflectionAngle|'+str(deflection_dict['target_angle'])+'_samplingBins|'+str(int(dt_step_resampling*1000))+'us'+deflection_dict['stims_string']+'.json', 'w') as fp: json.dump(spikes_dict, fp, indent=4)
            else:
                f_name=deflection_dict['save_deflection_model_folder']+'spike_dicts/mleSpikes_deflectionAngle|'+str(deflection_dict['target_angle'])+'_samplingBins|'+str(int(dt_step_resampling*1000))+'us'+deflection_dict['stims_string']+'.pkl'
                with open(f_name, 'wb') as f: pd.to_pickle(spikes_dict, f)

            if UpdateSampledRasterFigs: SampleData.analyzeDeflectionEvent(spikes_dict)
        
        # --- Load saved raster
        else:
            try:
                print('Loading dataset')
                # obs: dt_step_resampling is redefined, so that when loading you can skip updating the deflection model
                file_path = save_deflection_model_folder+'spike_dicts/mleSpikes_deflectionAngle|'+str(deflection_dict['target_angle'])+'_samplingBins|'+str(int(dt_step_resampling*1000))+'us'+deflection_dict['stims_string']+'.pkl'
                # spikes_dict = SampleData.LoadJSON(file_path)
                spikes_dict = SampleData.LoadPickle(file_path)
                if UpdateSampledRasterFigs: SampleData.analyzeDeflectionEvent(spikes_dict)
                print('Loading dataset - worked')
                print('Loaded spike dict - target angle: ',spikes_dict['target_angle'])
            except:
                print('Loading dataset - failed')

        # Saves an aditional dictionary dataset with only the spike times organized by GID - lighter, to be synchronized in Github
        try:
            if saveSimplifiedRasterData:
                simplified_spikes_dict={'spkts':spikes_dict['spkts']}
                # obs: dt_step_resampling is redefined, so that when loading you can skip updating the deflection model
                f_name=deflection_dict['save_deflection_model_folder']+'spike_dicts/mleSpikes_deflectionAngle|'+str(deflection_dict['target_angle'])+'_samplingBins|'+str(int(dt_step_resampling*1000))+'us'+deflection_dict['stims_string']+'_simplifiedDataset'+'.pkl'
                with open(f_name, 'wb') as f: pd.to_pickle(simplified_spikes_dict, f)
                print('Saving simplified dataset - worked')
        except:
            print('Saving simplified dataset - failed')

    #######################################################################################################################################
