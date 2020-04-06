# EMG Simulator for intramuscular multichannel shifting electrodes

## Prerequisites
No specific prerequisites needed, except for working MEX compilation for C/C++ code. See details in the "Intalling" section below. Code tested on Matlab 2018a. Important parts of code should work in other versions. Otherwise, please, raise an issue here.

## Installing
MEX should be set up on your Matlab. Run 'mex -setup', if no C compiler is recognized, follow the Mathorks' official instructions on setting up MEX.
MEX is needed to compile A. Peyre's "Geodesic farthest sampling" library (see 'Authors' section below). For that, in Matlab, go to folder './dependancies/geodesic_farthest_point/toolbox_graph', and run 'compile_mex'. In case of troubles with compiling it, please first refer to this guide: https://fr.mathworks.com/matlabcentral/fileexchange/6110-toolbox-fast-marching.

## Running the tests
The model is pre-defined with some general parameters. It is sufficient to run main.m and wait for the results (may take some time). For better understanding of the model workflow, you may run main.m in a section-by-section manner.

## Authors
Simulation strategy behing this code is described in our article: 
K. Akhmadeev, T. Yu, E. Le Carpentier, Y. Aoustin and D. Farina, "Simulation of motor unit action potential recordings from intramuscular multichannel scanning electrodes," in IEEE Transactions on Biomedical Engineering. https://ieeexplore.ieee.org/document/8926536

Third-party tools: 
- Toolbox Fast Marching: Gabriel Peyre (2020). Toolbox Fast Marching (https://www.mathworks.com/matlabcentral/fileexchange/6110-toolbox-fast-marching), MATLAB Central File Exchange. Retrieved April 6, 2020. Provided here for simplicity of configuration, otherwise available under the link.

## License
Packages provided in './dependancies' have their own licences. Please, see the corresponding .txt files.


