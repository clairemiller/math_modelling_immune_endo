# Journal article

This repository contains code to reproduce the results from the paper: **Mathematical modelling of macrophage and natural killer cell immune response during early stages of peritoneal endometriosis lesion onset**, currently in submission and available on BioRxiv (https://doi.org/10.1101/2025.01.20.633967).

## Abstract
The immune system is hypothesised to contribute to the onset of endometriosis lesions. However, the precise mechanisms underlying its role are not yet known. We introduce a novel compartmental model that describes the interactions between innate immune cells, specifically macrophages and natural killer cells, and endometrial cells, occurring within the peritoneal fluid during the early stages of superficial peritoneal) endometriosis lesion onset. Our study focuses on retrograde influx, immune detection, and immune clearance. Results show an increased influx of endometrial cells into peritoneal fluid correlates with heightened pro-inflammatory macrophage activation, but does not lead to an increase in disease. We compare the system's response to changes in immune cytotoxicity and ability to detect ectopic endometrial cells. We predict that reduced cytotoxicity is a key driver of disease. These findings align with the increased immune activation observed clinically. Lastly, we predict that an individual can transition to a diseased state following a reduction in immune system cytotoxicity and/or reduced ability to detect ectopic cells. Due to hysteresis, a significant improvement is then required to restore an individual to the disease-free state. This work provides a valuable framework to explore hypotheses of endometriosis lesion onset and assist in understanding of the disease.


## Citation

If you use this code, please cite the paper:

Miller, C., Germano, D. P. J., Chenoweth, A. M., & Holdsworth-Carson, S. (2025). *Mathematical modelling of the immune response during endometriosis lesion onset*. **bioRxiv**. https://doi.org/10.1101/2025.01.20.633967.

**bibtex:**
```bibtex
@article{miller25,
  title = {Mathematical Modelling of the Immune Response during Endometriosis Lesion Onset},
  author = {Miller, Claire and Germano, Domenic P. J. and Chenoweth, Alicia M. and {Holdsworth-Carson}, Sarah},
  year = {2025},
  pages = {2025.01.20.633967},
  publisher = {bioRxiv},
  doi = {10.1101/2025.01.20.633967},
}
```

# Project requirements

## Repository structure

```bash
├── README.md
├── auto_code
│   ├── one_parameter
│   │   ├── c.immune_model.beta1
│   │   ├── c.immune_model.omega
│   │   ├── c.immune_model_N1
│   │   ├── c.immune_model_N2
│   │   ├── c.immune_model_P2
│   │   ├── immune_model_N1.f90
│   │   ├── immune_model_N2.f90
│   │   ├── immune_model_P1.f90
│   │   ├── immune_model_P2.f90
│   │   └── run_1d_bifurcation.py
│   └── two_parameter
│       ├── c.immune_model
│       ├── immune_model.f90
│       ├── parse_fort9.py
│       └── run_2D_bifurcation.py
└── matlab_code
    ├── default_parameters.m
    ├── fn_determine_condition_satisfied.m
    ├── fn_immune_ode_system.m
    ├── script_example_systems.m
    └── script_parameter_sweep.m
```

## Software and libraries

The following is a list of the software versions and any required packages that were used to run the code in this project:

- **Matlab** R2024B
- **Python** 3.9.12, with required additional packages:
    - numpy (Version: 1.23)
    - pandas (Version: 2.2.3)
- **AUTO-07p**

### AUTO-07p
To perform the bifurcation analyses, you will need to have AUTO-07p installed. This is available from https://github.com/auto-07p/auto-07p/. To run the code in this project, the application directory is assumed to be located in the `auto_code` directory with folder name `auto-07p`.

# Running the code

## Simulation results (Section 3.2, Matlab)
The ODE system for the model was simulated using Matlab R2024B, and the code to do this can be found in the `matlab_code` directory. All output from the scripts is saved to the `output` directory. Results are stored in csv files with an associated `model_info.yml` detailing the default parameters and initial conditions used. 

### Time-series results (Section 3.2)

**Script**: `matlab_code/script_example_systems.m`

This script will generate the csv results for all the timeseries plots in Section 3.2 Hypotheses of disease onset: 

- the default system response (Fig. 1b)
- attachment examples (Fig. 2)
- reponse to varying rho0 (Fig. 3), 
- heatmap examples shown in the immune dysfunction results (Fig. 4a and 5c).

### Parameter sweep results (Section 3.2.2)

**Script**: `matlab_code/script_parameter_sweep.m`

This script runs (in parallel) the results for the heatmaps in Section 3.2.2 Immune Dysfunction. The script will generate a csv containing the parameter values and summary statistics for each simulation. The results in the csv files are for the heatmaps in the paper:

- Endometrial cell attachment for varying $\omega$ and $\beta_1$ (Fig. 4b)
- Immune cell activation for varying $\omega$ and $\beta_1$ (Fig. 5a/b.) 

**Script**: `matlab_code/script_parameter_sweep_rhoF.m`

This script runs (in parallel) the results for the heatmaps in Supplementary S4: varying $\omega$ and $\beta_1$ for different endometrial attachment rates.

## Bifurcation analysis (Section 3.3, AUTO-07p)

The bifurcation analysis was performed using AUTO-07p, using the python interface, and the python code to reproduce the results can be found in the `auto_code` directory. The scripts below will generate the results for Section 3.3 Disease Transition, Supplementary S5 Bifurcation diagrams, and Supplementary S6 Analysis of upregulation of M2-type macrophages by attached endometrial cells. All output from the scripts is saved to the `output` directory.

**auto-07 directory:** The python code assumes the `auto-07p` program directory is a subdirectory in `auto_code`. If this is stored elsewhere, the location must be updated in line 9 of `auto_code/parse_auto.py` (noting this must be passed as an absolute path, which is currently resolved from a relative path in line 10). 

### Running Through Docker

We have included a Dockerfile to run the bifurcation code through a Docker container if desired, removing the requirement to do auto and python package installation natively on your machine. For instructions on installing Docker, see: <https://www.docker.com/get-started/>.

**Build the Docker Image:** The docker image must be built before containers can be run (only required to be done once). We assume here this is being run from within the `auto_code` directory (the folder containing the Dockerfile). Here we label our image `auto_analysis`. 

``` bash
docker build -t auto_analysis .
```

**Start a Docker Container:** The following command will create a container for interactive use. It will also bind the output directory in the container to the local output directory, which is required for the output files to be accessed after the container is exited. We assume here that the docker container is started from within the **main repository directory** (i.e. the relative path to the `output` directory is `./output`). 

``` bash
docker run --volume ./output:/output -it --rm auto_analysis
```

**Run the analysis:** The commands to run the analysis interactively through the docker container are the same as those detailed below. Note, the auto generated plots will not appear at the end of the analysis when run through Docker. To run the entire analysis the following commands would be run:

```bash
python run_1d_bifurcation.py -p omega --all
python run_1d_bifurcation.py -p beta1 --all
python run_2d_bifurcation.py --all
python run_m2_upregulation.py
```

**Exit and Stop the Container When Finished:** Exiting the container will also delete it as we included the `--rm` flag in the run command.

``` bash
exit
```

### One parameter bifurcation

**Script**: `auto_code/run_1d_bifurcation.py` 

There are 4 systems that need to be run to output all of the required data for the diagrams in the paper. \[P/N\]\[1/2\] Solve the (P)ositive or (N)egative version of the system to get the maximum and minimum bounds of the limit cycle, with outputs (1) covering states $M_0$ to $K_A$ and (2) $K_0$ to $E_A$ (as auto only outputs a maximum of 6 states). 

This script will generate the results in Fig. 6b and Supplementary Figures S3 and S4. To run this script, run the following commands in bash within the `auto_code` directory (this will also show the auto plot of the diagram at the end of the analysis):
 
- To run a single parameter (omega, $\omega$, or beta1, $\beta_1$) and system (P1, P2, N1, or N2): 
    ``` bash
    python run_1d_bifurcation.py -p [omega/beta1] -s [P1/P2/N1/N2]
    ```

- Run all of the systems for one parameter (omega, $\omega$, or beta1, $\beta_1$): 
    ```bash
    python run_1d_bifurcation.py -p [omega/beta1] --all
    ```

### Codimensional bifurcation
**Script**: `auto_code/run_2D_bifurcation.py` 

This script will generate the results for Fig. 6b and Supplementary Figure S5(a). To run this script, run the following command in bash within the `auto_code` directory: 
    
- To run only one level of $\rho_0$ (low, mid, or high, also shows auto plot at end of analysis):
    ```bash
    python run_2d_bifurcation.py -r [low/mid/high]
    ```

- To run all levels of $\rho_0$: 
    ```bash
    python run_2d_bifurcation.py --all
    ```

### M2-upregulation (Supplementary S6)
**Script**: `auto_code/run_m2_upregulation.py`

This script will generate the results for Fig. S6. To run this script, run the following command in bash within the `auto_code` directory:

```bash
python run_m2_upregulation.py
```