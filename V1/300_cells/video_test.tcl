# Visual Neuronal Dynamics
# Visualization State
::NeuronVND::loadFiles /home/fernando/Documents/sonata/examples/300_cells/circuit_config.json false false
::NeuronVND::loadSpikes /home/fernando/Documents/sonata/examples/300_cells/output/spikes.h5
# List of representations
::NeuronVND::createRepArgs show true style morphology_line color blue material Glass3 neurons 50 selection {stride 6}
::NeuronVND::createRepArgs show true style soma color Type material Opaque neurons 300 selection {all}
