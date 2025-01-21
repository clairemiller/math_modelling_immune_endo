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
from auto import run, load, plot, wait, clean

# Function to read and write the data from a bifDiagBranch object
def extract_data(branch, data_label):
    header = branch.keys()
    point_types = [pt['TY'] if pt['TY'] != 'No Label' else None
                for pt in branch]
    df = pd.DataFrame(branch.toArray(), columns=header)
    df['point_type'] = point_types
    df['branch_id'] = data_label
    return(df)

# Function to write the forwards and backwards data
def write_fwd_bwds(data_dir, sol_fwd, sol_bwd, fort9_df_fwd, fort9_df_bwd, label):
    # Combine the solutions
    sol = sol_fwd + sol_bwd

    # Work out the labelling
    filelab = label.lower()
    sol_labels = [(label+"_fwd"), (label+"_bwd")]

    # Extract the data and export to a csv file
    sol_df = [extract_data(sol[i], data_label = sol_labels[i]) for i in range(len(sol))]
    pd.concat(sol_df).to_csv(Path(data_dir, filelab + "_data.csv"), index = False)

    # Export the fort.9 data
    fort9_df_fwd['branch_id'] = label+'_fwd'
    fort9_df_bwd['branch_id'] = label+'_bwd'
    fort9_df = pd.concat([fort9_df_fwd, fort9_df_bwd])
    if 'multipliers' in fort9_df.columns:
        fort9_df = fort9_df.drop(columns=['multipliers'])
    fort9_df.to_csv(Path(data_dir, filelab + "_fort9.csv"), index = False)

def run_2d_bifurcation(system_chosen):

    # Make sure the system_chosen is one of low/mid/high
    assert system_chosen in ['low', 'mid', 'high'], "system_chosen must be one of 'low', 'mid', 'high'"

    # Data directory
    data_dir = Path("../data/bifurcation/codimensional",system_chosen)
    # Create directory if it doesn't exist
    data_dir.mkdir(parents=True, exist_ok=True)


    # Clean the directory and load the model

    #-----------------------------------------------------------------------------------------

    clean()

    # Load the model - files: immune_model_[low/mid/high]rho0.f90 and c.immune_model
    f90_filename = f"immune_model_{system_chosen}rho0"
    immune_model = load(f90_filename,c="immune_model")


    # Run 1D bifurcation in omega
    #-----------------------------------------------------------------------------------------

    # Output current stem
    print(rf"{'-'*100}")
    print("STABLE POINTS")
    print(rf"{'-'*100}")

    # Run the 1D bifurcation (omega) forwards
    omega_forward = run(immune_model,
                NMX  = 1000000, NPR  = 100000000,
                EPSL = 1e-05, EPSU = 1e-05, EPSS = 1e-04,
                DS = 1e-6, DSMIN = 1e-8, DSMAX = 1e-4,
                UZSTOP = {'beta1_scaled': [-10, -3], 'omega_scaled':[-7, -3]})
    # Extract the fort.9 data
    fort9_df_fwd = parse_fort9.parse_fort9(downsample=300)

    # Run the 1D bifurcation backwards
    omega_backward =  run(immune_model,
                        NMX  = 1000000, NPR  = 100000000,
                        EPSL = 1e-05, EPSU = 1e-05, EPSS = 1e-04,
                        DS = -1e-6, DSMIN = 1e-8, DSMAX = 1e-4,
                        UZSTOP = {'beta1_scaled': [-10, -3], 'omega_scaled':[-7, -3]})
    # Extract the fort.9 data
    fort9_df_bwd = parse_fort9.parse_fort9(downsample=100)

    # Combine the solutions
    omega = omega_forward + omega_backward

    # Write the output
    write_fwd_bwds(data_dir, omega_forward, omega_backward, fort9_df_fwd, fort9_df_bwd, 'SP_omega')


    # Run from first hopf bifurcation point
    #-----------------------------------------------------------------------------------------
    print(rf"{'-'*100}")
    print("HOPF BIFURCATION 1")
    print(rf"{'-'*100}")

    hb1_start = load(omega('HB1'), ISW=1)

    hb1_forward = run(hb1_start, IPS=2,
                    ICP=['omega_scaled', 'PERIOD'],
                    THL =  {'PERIOD': 0.10},
                    NTST = 100, NCOL = 5, NMX=40000, NPR=100000,
                    DS=0.005, DSMAX = 1e-1,  MXBF = 1,
                    EPSL = 1e-4, EPSU = 1e-4, EPSS = 1e-2,
                    SP = ['BP1'])
    fort9_df_fwd = parse_fort9.parse_fort9(downsample=10)

    hb1_backward = run(hb1_start, IPS=2,MXBF = 1,
                    ICP=['omega_scaled', 'PERIOD'],
                    THL =  {'PERIOD': 0.10},
                    NTST = 100, NCOL = 5, NMX=10000, NPR=100000,
                    DS=-0.005, DSMAX = 1e-1,
                    EPSL = 1e-4, EPSU = 1e-4, EPSS = 1e-2,
                    SP = ['BP1'])
    fort9_df_bwd = parse_fort9.parse_fort9()

    # Write the forward and backward solutions
    write_fwd_bwds(data_dir, hb1_forward, hb1_backward, fort9_df_fwd, fort9_df_bwd, 'HB1_omega')

    # Combine the solutions
    hb1 = hb1_forward + hb1_backward


    # Run from second hopf bifurcation point (not actually used for 2D)
    #-----------------------------------------------------------------------------------------
    # print(rf"{'-'*100}")
    # print("HOPF BIFURCATION 2")
    # print(rf"{'-'*100}")
    # hb2_start = load(omega('HB2'), ISW=1)

    # hb2_forward = run(hb2_start, IPS=2, MXBF = 1,
    #                   ICP=['omega_scaled', 'PERIOD'],
    #                   THL =  {'PERIOD': 0.10},
    #                   NTST = 100, NCOL = 5, NMX=20000, NPR=100000,
    #                   DS=0.005, DSMAX = 1e-1,
    #                   EPSL = 1e-4, EPSU = 1e-4, EPSS = 1e-2,
    #                   SP = ['BP1'])
    # fort9_df_fwd = parse_fort9.parse_fort9(downsample=5)

    # hb2_backward = run(hb2_start, IPS=2, MXBF = 1,
    #                   ICP=['omega_scaled', 'PERIOD'], THL =  {'PERIOD': 0.10},
    #                   NTST = 100, NCOL = 5, NMX=10000, NPR=10000,
    #                   DS=-0.005, DSMAX = 1e-1,
    #                   EPSL = 1e-4, EPSU = 1e-4, EPSS = 1e-2,
    #                   SP = ['BP1'])
    # fort9_df_bwd = parse_fort9.parse_fort9()

    # # Combine the solutions
    # hb2 = hb2_forward + hb2_backward

    # # Write
    # write_fwd_bwds(hb2_forward, hb2_backward, fort9_df_fwd, fort9_df_bwd, 'HB2_omega')

    # # plot solution
    # plot(hb2 + hb1 + omega)


    # Run from first limit point of hopf bifurcation 1 in 2D
    #-----------------------------------------------------------------------------------------
    print(rf"{'-'*100}")
    print("FROM LIMIT POINT OF HB1 - 2D")
    print(rf"{'-'*100}")
    # Set the new start label to the first LP label of hopf bifurcation 1
    lpbp1 = load(hb1('LP1'), ISW=2)

    # Continue from this label in two parameters
    lpstart = run(lpbp1,DS= 1e-4 ,
                    ISW =2, IPS=2, IRS =   5, ILP =   1,
                    ICP=['omega_scaled', 'beta1_scaled', 'PERIOD'],
                    THL =  {'PERIOD': 0.0})

    loci_1 = run(lpstart, NMX=100000, NPR=1000000, 
                DS = 1e-4, DSMIN = 1e-9, DSMAX = 1e-4,
                EPSL = 1e-4, EPSU = 1e-4, EPSS = 1e-2,
                SP = ['BP2'])
    fort9_df_fwd = parse_fort9.parse_fort9(downsample=100)

    loci_2 = run(lpstart, NMX=100000, NPR=1000000, 
                DS =-1e-4, DSMIN = 1e-9, DSMAX = 1e-4,
                EPSL = 1e-4, EPSU = 1e-4, EPSS = 1e-2,
                SP = ['BP2'])
    fort9_df_bwd = parse_fort9.parse_fort9(downsample=10)

    # Write
    write_fwd_bwds(data_dir, loci_1, loci_2, fort9_df_fwd, fort9_df_bwd, 'LP1_HB1_2D')

    # Combine the solutions and plot
    lp1 = loci_1 + loci_2
    # p = plot(lp1)
    # p.config(bifurcation_x=['omega_scaled'])
    # p.config(bifurcation_y=['beta1_scaled'])


    # Continuing from limit point 1 of SPs
    #-----------------------------------------------------------------------------------------
    print(rf"{'-'*100}")
    print("FROM LIMIT POINT 1 OF SP - 2D")
    print(rf"{'-'*100}")

    # Set the new start label to the first LP label
    lp1 = load(omega('LP1'), ISW=2)

    # Continue from this label in two parameters
    lp1_twoDim_fwd = run(lp1,DS= 1e-2 )
    fort9_df_fwd = parse_fort9.parse_fort9(downsample=200)
    lp1_twoDim_bwd = run(lp1, DS= -1e-2  )
    fort9_df_bwd = parse_fort9.parse_fort9(downsample=500)

    # Write
    write_fwd_bwds(data_dir, lp1_twoDim_fwd, lp1_twoDim_bwd, fort9_df_fwd, fort9_df_bwd, 'LP1_2D')

    # Combine and plot
    lp1_twoDim = lp1_twoDim_fwd + lp1_twoDim_bwd
    # p = plot(lp1_twoDim)
    # p.config(bifurcation_x=['omega_scaled'])
    # p.config(bifurcation_y=['beta1_scaled'])

    # Continuing from limit point 1 of SPs
    #-----------------------------------------------------------------------------------------
    print(rf"{'-'*100}")
    print("FROM LIMIT POINT 2 OF SP - 2D")
    print(rf"{'-'*100}")

    # Set the new start label to the first LP label in b.mu and s.mu
    lp2 = load(omega('LP2'), ISW=2)

    # Continue from this label in two parameters
    lp2_twoDim_fwd =  run(lp2, DS= 1e-2)
    fort9_df_fwd = parse_fort9.parse_fort9(downsample=200)
    lp2_twoDim_bwd =  run(lp2,DS= -1e-2)
    fort9_df_bwd = parse_fort9.parse_fort9(downsample=500)

    # Write
    write_fwd_bwds(data_dir, lp2_twoDim_fwd, lp2_twoDim_bwd, fort9_df_fwd, fort9_df_bwd, 'LP2_2D')

    # Combine and plot
    lp2_twoDim = lp2_twoDim_fwd + lp2_twoDim_bwd
    # # p = plot(lp2_twoDim)
    # # p.config(bifurcation_x=['omega_scaled'])
    # # p.config(bifurcation_y=['beta1_scaled'])


    # Continuing from HB1 in two dimensions
    #-----------------------------------------------------------------------------------------
    print(rf"{'-'*100}")
    print("FROM HPF BIFURCATION 1 IN 2D")
    print(rf"{'-'*100}")

    # Set the new start label to the first LP label
    hb1 = load(omega('HB1'), ISW=2)

    # Continue from this label in two parameters
    hb1_twoDim_fwd = run(hb1,DS= 1e-2)
    fort9_df_fwd = parse_fort9.parse_fort9(downsample=200)
    hb1_twoDim_bwd = run(hb1, DS= -1e-2)
    fort9_df_bwd = parse_fort9.parse_fort9(downsample=300)

    # Write
    write_fwd_bwds(data_dir, hb1_twoDim_fwd, hb1_twoDim_bwd, fort9_df_fwd, fort9_df_bwd, 'HB1_2D')

    # Combine and plot
    hb1_twoDim = hb1_twoDim_fwd + hb1_twoDim_bwd
    # # p = plot(hb1_twoDim)
    # # p.config(bifurcation_x=['omega_scaled'])
    # # p.config(bifurcation_y=['beta1_scaled'])


    # Continuing from HB1 in two dimensions
    #-----------------------------------------------------------------------------------------
    print(rf"{'-'*100}")
    print("FROM HPF BIFURCATION 2 IN 2D")
    print(rf"{'-'*100}")

    # Set the new start label to the first LP label
    hb2 = load(omega('HB2'), ISW=2)
    # Continue from this label in two parameters
    hb2_twoDim_fwd =  run(hb2,DS= 1e-2)
    fort9_df_fwd = parse_fort9.parse_fort9(downsample=50)
    hb2_twoDim_bwd =  run(hb2, DS= -1e-2)
    fort9_df_bwd = parse_fort9.parse_fort9(downsample=100)

    # Write
    write_fwd_bwds(data_dir, hb2_twoDim_fwd, hb2_twoDim_bwd, fort9_df_fwd, fort9_df_bwd, 'HB2_2D')

    # Combine and plot
    hb2_twoDim = hb2_twoDim_fwd + hb2_twoDim_bwd
    # # p = plot(hb2_twoDim)
    # # p.config(bifurcation_x=['omega_scaled'])
    # # p.config(bifurcation_y=['beta1_scaled'])



    # Combine all the solutions for plotting
    #-----------------------------------------------------------------------------------------
    return hb2_twoDim + hb1_twoDim + lp2_twoDim + lp1_twoDim + loci_1 + loci_2

def main():
    parser = argparse.ArgumentParser(description="Run 2D bifurcation analysis for varying rho0 levels.")
    parser.add_argument('-r', '--rho0', default='mid', choices=['low', 'mid', 'high'],
                        help='Level of rho0 to run (default: mid).')
    parser.add_argument('--all', action='store_true', 
                        help='Use all systems for bifurcation analysis for parameter in -p')

    # Get the specified parameters
    args = parser.parse_args()
    system_chosen = args.rho0

    # Run the bifurcation analysis
    if args.all:
        for system in ['low', 'mid', 'high']:
            print("Running bifurcation analysis for system: ", system)
            # Continue to next system even if there is an error
            try:
                run_2d_bifurcation(system)
            except:
                print("Error occured for system: ", system)
                continue
    else:
        print("Running bifurcation analysis for system: ", system_chosen)
        bd = run_2d_bifurcation(system_chosen)
        plot(bd)
        wait()

if __name__ == "__main__":
    #run_2d_bifurcation("mid")
    main()