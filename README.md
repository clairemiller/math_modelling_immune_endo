# Mathematical modelling of the immune response during endometriosis lesion onset

This repository contains the code required to reproduce the results from our paper "Mathematical modelling of the immune response during endometriosis lesion onset", currently in submission and available on BioRxiv.

## Simulations of ODE system

The ODE system for the model was simulated using Matlab, and the code to do this is available in the matlab folder. The scripts in this folder as as follows:

1. `script_example_systems.m`: This script can be used to run a specific example system (first section) and to run all the timeseries (remaining sections) that are in the manuscript, including the default system (Fig. 1b), attachment examples (Fig. 2), reponse to varying rho0 (Fig. 3), and the examples associated with immune dysfunction section (Fig. 4a and 5b). To change output data folder from `../../data/`, change the directory paths in lines 57, 79, 116, 176.
2. `script_parameter_sweep.m`: This script runs (in parallel) the simulations used to generate the immune dysfunction heatmaps in Fig. 4b and Fig. 5a/b. To change the output folder change the directory name in `dirname`, line 3. 

## Bifurcation analysis

The bifurcation analysis was performed using AUTO-07p, using the python interface, and the code to do this is available in the AUTO folder. To run this code, you will need to have AUTO-07p, available from https://github.com/auto-07p/auto-07p/.

For the one parameter bifuraction analysis, script: `one_parameter/run_1d_bifurcation.py` 

1. Change the location of the auto-07p directory in the `auto_dir` variable in line 12. 
2. Specify the location for the output directory in the `data_dir` variable in line 36. 
3. There are 4 systems that need to be run to output all of the required data for the diagrams in the paper. \[P/N\]\[1/2\] Solve the (P)ositive or (N)egative version of the system to get the maximum and minimum bounds of the limit cycle, with outputs (1) covering states M0 to KA and (2) K0 to EA (as auto only outputs a maximum of 6 states). 

- To run all of the systems for one parameter (omega or beta1) run `python run_1d_bifurcation.py -p [omega/beta1] --all` - To run a single system run `python run_1d_bifurcation.py -p [omega/beta1] -s [P1/P2/N1/N2]`. This will also show the auto plot of the diagram at the end of the analysis.

For the co-dimensional bifurcation analysis, script: `two_parameter/run_2D_bifurcation.py` 

1. Change the location of the auto-07p directory in the `auto_dir` variable in line 12. 
2. Specify the location for the output directory in the `data_dir` variable in line 56. 
3. Run all levels of rho0: `python run_2D_bifurcation.py --all` 
4. To run only one level of rho0 (also shows auto plot at end of analysis): `python run_2D_bifurcation.py -r [low/mid/high]`