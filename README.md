# EMG Simulator for intramuscular multichannel shifting electrodes

## Description

## Prerequisites
No specific prerequisites needed, except for working MEX compilation for C/C++ code. See details in the "Intalling" section below. Code tested on Matlab 2018a. Important parts of code should work in other versions. Otherwise, please raise an issue here.

## Installing
MEX should be set up on your Matlab. Run 'mex -setup', if no C compiler is recognized, follow the Matlab official instructions on setting up MEX. 
A. Peyre's "Geodesic farthest sampling" library should be compiled using MEX. For that, in Matlab, go to ./dependancies/geodesic_farthest_point/toolbox_graph, and run 'compile_mex'. Otherwise, you may also follow this guide: https://fr.mathworks.com/matlabcentral/fileexchange/6110-toolbox-fast-marching .

## Running the tests
The model is pre-defined with quite general parameters. It is sufficient to run main.m and wait for the results (may take some time). For better understanding of the model workflow, you may run main.m in a section-by-section manner.

## Authors
Konstantin Akhmadeev, Tianyi Yu, Eric Le Carpentier, Yannick Aoustin, Dario Farina

## License
To be defined. 

## Acknowledgments

