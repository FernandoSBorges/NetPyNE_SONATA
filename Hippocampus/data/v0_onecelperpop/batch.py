"""
batch.py 

Influence if slow potassium and Ca channels in bistable firing patterns using NetPyNE

Contributors: protachevicz@gmail.com, fernandodasilvaborges@gmail.com
"""
from netpyne.batch import Batch
from netpyne import specs
import numpy as np

# ----------------------------------------------------------------------------------------------
# Custom
# ----------------------------------------------------------------------------------------------
def custom():
    params = specs.ODict()
    
    params[('seeds', 'conn')] =  [1] 
    # params[('EEGain')] = [0.5, 1, 2, 4]
    # params[('IClamp1')] = [0.001] 

    b = Batch(params=params, netParamsFile='netParams.py', cfgFile='cfg.py')

    return b

# ----------------------------------------------------------------------------------------------
# Run configurations
# ----------------------------------------------------------------------------------------------
def setRunCfg(b, type='mpi_bulletin'):
    if type=='mpi_bulletin' or type=='mpi':
        b.runCfg = {'type': 'mpi_bulletin', 
            'script': 'init.py', 
            'skip': True}

    elif type=='mpi_direct':
        b.runCfg = {'type': 'mpi_direct',
            'cores': 6,
            'script': 'init.py',
            'mpiCommand': 'mpiexec --use-hwthread-cpus', # i7  --use-hwthread-cpus
            'skip': True}

    elif type=='hpc_slurm_gcp':
        b.runCfg = {'type': 'hpc_slurm',
            'allocation': 'default',
            'walltime': '24:00:00',
            'nodes': 1,
            'coresPerNode': 80,
            'email': 'fernandodasilvaborges@gmail.com',
            'folder': '/home/ext_fernandodasilvaborges_gmail_/S1_mouse/sim/',
            'script': 'init.py',
            'mpiCommand': 'mpirun',
            'skipCustom': '_raster.png'}

# ----------------------------------------------------------------------------------------------
# Main code
# ----------------------------------------------------------------------------------------------
if __name__ == '__main__': 
    b = custom() #

    b.batchLabel = 'v0_batch10'  
    b.saveFolder = '../data/'+b.batchLabel
    b.method = 'grid'
    setRunCfg(b, 'mpi_direct')     # setRunCfg(b, 'mpi_bulletin')
    b.run() # run batch
