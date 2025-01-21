# Libraries
import numpy as np
import pandas as pd
from pathlib import Path
import argparse
import parse_fort9
from ast import literal_eval

# Need to add the auto-07p paths to our environment to import the AUTO python modules
import sys, os
# os.chdir('immune_model_final') # If running interactively
auto_dir = "../auto-07p"
os.environ['AUTO_DIR'] = auto_dir
sys.path.append(Path(auto_dir,"bin"))
sys.path.append(Path(auto_dir,"cmds"))
import auto
from auto import run, load

def run_1d_bifurcation(parameter, system_chosen):
    # Dictionary to map the system and parameter to the correct model and parameter names
    systems_list = {"P1": 'immune_model_P1', "P2": 'immune_model_P2', 
                    "N1": 'immune_model_N1', "N2": 'immune_model_N2'}
    constants_file = "immune_model." + parameter
    parameter_lab_auto = parameter + "_scaled"

    # Change the parameter label if we are looking at the rho0 system for output location
    if system_chosen == "lowrho0":
        parameter = parameter + "_lowrho0"

    # Correct name ordering for the parameters
    unames_1 = {1: 'M0', 2: 'M1', 3: 'M2', 4: 'K0', 5: 'KA', 6: 'E0', 7: 'EF', 8: 'EA'}
    unames_2 = {(9-k): v for k, v in unames_1.items()}
    unames_dict = {"P1": unames_1, "P2": unames_2, "N1": unames_1, "N2": unames_2}

    # Data directory
    data_dir = Path("../../data/bifurcation",parameter)
    # Create directory if it doesn't exist
    data_dir.mkdir(parents=True, exist_ok=True)

        # Function to read and write the data from a bifDiagBranch object
    def extract_data(branch, data_label):
        header = branch.keys()
        point_types = [pt['TY'] if pt['TY'] != 'No Label' else None
                    for pt in branch]
        df = pd.DataFrame(branch.toArray(), columns=header)
        df['point_type'] = point_types
        df['system_id'] = system_chosen
        df['branch_id'] = data_label
        return(df)

    # Load the model - files: immune_model.f90 and c.immune_model
    immune_model = load(systems_list[system_chosen], c=constants_file)

    print(rf"{'-'*100}")
    print(system_chosen + " : STABLE POINTS")
    print(rf"{'-'*100}")

    # Run the 1D bifurcation (omega) forwards
    omega_forward = run(immune_model, NMX  = 1000000, NPR  = 100000000,
                        unames = unames_dict[system_chosen],
                        EPSL = 1e-05, EPSU = 1e-05, EPSS = 1e-04,
                        DS = 1e-6, DSMIN = 1e-8, DSMAX = 5e-4,
                        UZSTOP = {'beta1_scaled': [-10, -4], 'omega_scaled':[-7, -3]},
                        JAC = 0)
    # Extract the fort.9 data
    fort9_df_fwd = parse_fort9.parse_fort9(downsample=100)

    # Run the 1D bifurcation backwards
    omega_backward =  run(immune_model, NMX  = 200000, NPR  = 2000000,
                        unames = unames_dict[system_chosen],
                        DS = -1e-6, DSMIN = 1e-8, DSMAX = 5e-5,
                        UZSTOP = {'beta1_scaled': [-10, -4], 'omega_scaled':[-7, -3]}, 
                        JAC = 0)
    # Extract the fort.9 data
    fort9_df_bwd = parse_fort9.parse_fort9(downsample=100)

    # Combine the solutions
    omega = omega_forward + omega_backward
    labels = ['SP_fwd', 'SP_bwd']

    # Extract the data and export to a csv file
    omega_df = [extract_data(omega[i], data_label = labels[i]) for i in range(len(omega))]
    pd.concat(omega_df).to_csv(Path(data_dir,system_chosen+"_sp_data.csv"), index = False)

    # Export the fort.9 data (we don't need to keep the eigenvalues)
    fort9_df_fwd['branch_id'] = 'SP_fwd'
    fort9_df_bwd['branch_id'] = 'SP_bwd'
    fort9_df = pd.concat([fort9_df_fwd, fort9_df_bwd])
    
    # eigenvalues_df = fort9_df[[parameter_lab_auto, 'TY', 'stability', 'eigenvalues']]
    # eigenvalues_df.to_csv(Path(data_dir,system_chosen+"_sp_eigenvalues.csv"), index = False)
    
    fort9_df = fort9_df.drop(columns=['eigenvalues'])
    fort9_df['system_id'] = system_chosen
    fort9_df.to_csv(Path(data_dir,system_chosen+"_sp_fort9.csv"), index = False)


    # Run from first hopf bifurcation point
    #-----------------------------------------------------------------------------------------
    print(rf"{'-'*100}")
    print(system_chosen + " : HOPF BIFURCATION 1")
    print(rf"{'-'*100}")

    hb1_start = load(omega('HB1'), ISW=1)
    # Run forward
    hb1_forward = run(hb1_start, IPS=2, MXBF = 1,
                    ICP=[parameter_lab_auto, 'PERIOD'],
                    THL =  {'PERIOD': 0.10},
                    NTST = 100, NCOL = 5, NMX=400000, NPR=100000,
                    DS=0.001, DSMAX = 1e-1,
                    EPSL = 1e-4, EPSU = 1e-4, EPSS = 1e-2,
                    SP = ['BP1'],
                    JAC = 0)
    fort9_df_fwd = parse_fort9.parse_fort9()

    # Run backwards
    hb1 = hb1_forward 
    hb1_labels = ['HB1_fwd', 'HB1_bwd']

    # Extract the data and export to a csv file
    hb1_df = [extract_data(hb1[i], data_label = hb1_labels[i]) for i in range(len(hb1))]
    pd.concat(hb1_df).to_csv(Path(data_dir,system_chosen+"_hb1_data.csv"), index = False)

    # Export the fort.9 data
    fort9_df_fwd['branch_id'] = 'HB1_fwd'
    fort9_df_bwd['branch_id'] = 'HB1_bwd'
    fort9_df = fort9_df_fwd #pd.concat([fort9_df_fwd, fort9_df_bwd])
    fort9_df = fort9_df.drop(columns=['multipliers'])
    fort9_df['system_id'] = system_chosen
    fort9_df.to_csv(Path(data_dir,system_chosen+"_hb1_fort9.csv"), index = False)


    # Run from the second hopf bifurcation point for beta1 (omega has only 1 hopf point)
    #-----------------------------------------------------------------------------------------
    if (parameter == "omega"):
        return omega + hb1
    
    print(rf"{'-'*100}")
    print(system_chosen + " : HOPF BIFURCATION 2")
    print(rf"{'-'*100}")
    hb2_start = load(omega('HB2'), ISW=1)
    # Run forward
    hb2_forward = run(hb2_start, IPS=2, MXBF = 1,
                    ICP=[parameter_lab_auto, 'PERIOD'],
                    THL =  {'PERIOD': 0.10},
                    NTST = 100, NCOL = 5, NMX=20000, NPR=100000,
                    DS=0.001, DSMAX = 1e-1,
                    EPSL = 1e-4, EPSU = 1e-4, EPSS = 1e-2,
                    SP = ['BP1'],
                    JAC = 0)
    fort9_df_fwd = parse_fort9.parse_fort9()

    # Combine the solutions
    hb2 = hb2_forward 
    hb2_labels = ['HB2_fwd', 'HB2_bwd']

    # Extract the data and export to a csv file
    hb2_df = [extract_data(hb2[i], data_label = hb2_labels[i]) for i in range(len(hb2))]
    pd.concat(hb2_df).to_csv(Path(data_dir, system_chosen+"_hb2_data.csv" ), index = False)

    # Export the fort.9 data (we don't need to output the multipliers data)
    fort9_df_fwd['branch_id'] = 'HB2_fwd'
    fort9_df_bwd['branch_id'] = 'HB2_bwd'
    fort9_df = fort9_df_fwd #pd.concat([fort9_df_fwd, fort9_df_bwd])
    fort9_df = fort9_df.drop(columns=['multipliers'])
    fort9_df['system_id'] = system_chosen
    fort9_df.to_csv(Path(data_dir,system_chosen+"_hb2_fort9.csv"), index = False)

    return omega + hb1 + hb2

def main():
    parser = argparse.ArgumentParser(description="Run 1D bifurcation analysis.")
    parser.add_argument('-p', '--parameter', default='omega', choices=['omega', 'beta1'], 
                        help='Parameter to use for bifurcation analysis (default: omega)')
    parser.add_argument('-s', '--system', default='P1', choices=['P1', 'P2', 'N1', 'N2'],
                        help='System (defined by .f90 files) to use for bifurcation analysis (degault: P1)')
    parser.add_argument('--all', action='store_true', 
                        help='Use all systems for bifurcation analysis for parameter in -p')

    # Get the specified parameters
    args = parser.parse_args()
    parameter = args.parameter
    system_chosen = args.system

    # Run the bifurcation analysis
    if args.all:
        for system in ['P1', 'P2', 'N1', 'N2']:
            print("Running bifurcation analysis for system: ", system)
            # Continue to next system even if there is an error
            try:
                run_1d_bifurcation(parameter, system)
            except:
                print("Error occured for system: ", system)
                continue
    else:
        print("Running bifurcation analysis for system: ", system_chosen)
        bd = run_1d_bifurcation(parameter, system_chosen)
        auto.plot(bd)
        auto.wait()

if __name__ == "__main__":
    # run_1d_bifurcation("beta1", "P1")
    main()