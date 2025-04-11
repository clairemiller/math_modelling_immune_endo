# Libraries
import pandas as pd
from pathlib import Path
import os
# Must import parse_auto before we import auto
from parse_auto import parse_fort9, extract_data
import auto
from auto import run, load

def run_1d_bifurcation(system):

    # Data directory
    data_dir = Path("../../output/bifurcation/m2_upregulation")

    # Create directory if it doesn't exist
    data_dir.mkdir(parents=True, exist_ok=True)

    # Choose the correct constants file
    c_file = "immune_model.alpha"
    f90_file = "immune_model_" + system
    # Load the model
    immune_model = load(f90_file, c=c_file)

    # Sort out the names order
    unames_P1 = {1: 'M0', 2: 'M1', 3: 'M2', 4: 'K0', 5: 'KA', 6: 'E0', 7: 'EF', 8: 'EA'}
    unames_P2 = {(9-k): v for k, v in unames_P1.items()}
    unames_order = unames_P1
    if ("P2" in system or "N2" in system):
        unames_order = unames_P2

    print(rf"{'-'*100}")
    print(system + " : STABLE POINTS")
    print(rf"{'-'*100}")

    # Run the stationary points
    #-----------------------------------------------------------------------------------------

    # Run the bifurcation forwards
    fwd = run(immune_model, NMX  = 10000000, NPR  = 100000000,
              unames = unames_order,
              EPSL = 1e-05, EPSU = 1e-05, EPSS = 1e-05,
              DS = 1e-10, DSMIN = 1e-10, DSMAX = 5e-4,
              UZSTOP = {'alpha_log': [-6, 9], 'omega_scaled':[-7, -2]},
              JAC = -1)
    # Extract the fort.9 data
    fort9_df_fwd = parse_fort9(downsample=100)
    
    # Run the 1D bifurcation backwards
    bwd = run(immune_model, 
              unames = unames_order,
              NMX  = 10000000, NPR  = 2000000,
              DS = -1e-10, DSMIN = 1e-8, DSMAX = 5e-5,
              UZSTOP = {'alpha_log': [-6, 9], 'omega_scaled':[-7, -2]},
              JAC = -1)
    # Extract the fort.9 data
    fort9_df_bwd = parse_fort9(downsample=250)

    # Combine the solutions
    bd = fwd + bwd
    labels = ['fwd', 'bwd']

    # Extract the data and export to a csv file
    bd_df = [extract_data(bd[i], data_label = labels[i], system=system) for i in range(len(bd))]
    pd.concat(bd_df).to_csv(Path(data_dir,system+"_sp_data.csv"), index = False)
    # Export the fort.9 data (we don't need to keep the eigenvalues)
    fort9_df_fwd['branch_id'] = 'fwd'
    fort9_df_bwd['branch_id'] = 'bwd'
    fort9_df = pd.concat([fort9_df_fwd, fort9_df_bwd])
    fort9_df = fort9_df.drop(columns=['eigenvalues'])
    fort9_df['system_id'] = system
    fort9_df.to_csv(Path(data_dir,system+"_sp_fort9.csv"), index = False)

    # Return if not the high manifold no hopf bifurcation
    if (not ("high" in system)):
        return(bd)


    # Run from first hopf bifurcation point
    #-----------------------------------------------------------------------------------------

    print(rf"{'-'*100}")
    print(system + " : HOPF BIFURCATION 1")
    print(rf"{'-'*100}")

    hb1_start = load(bd('HB1'), ISW=1)

    # Run forward
    hb1_fwd = run(hb1_start, IPS=2, MXBF = 1,
                    ICP=["alpha_log", 'PERIOD'],
                    THL =  {'PERIOD': 0.0},
                    NTST = 100, NCOL = 5, NMX=400000, NPR=100000,
                    DS=0.001, DSMAX = 1e-1,
                    EPSL = 1e-4, EPSU = 1e-4, EPSS = 1e-2,
                    SP = ['BP1'],
                    JAC = 0)
    fort9_df_hb1_fwd = parse_fort9()

    # Run backward
    hb1_bwd = run(hb1_start, IPS=2,MXBF = 1,
                ICP=["alpha_log", 'PERIOD'],
                THL =  {'PERIOD': 0.0},
                NTST = 100, NCOL = 5, NMX=10000, NPR=100000,
                DS=0.001, DSMAX = 1e-1,
                EPSL = 1e-4, EPSU = 1e-4, EPSS = 1e-2,
                SP = ['BP1'],
                JAC = 0)
    fort9_df_hb2_bwd = parse_fort9()

    # Combine the solutions
    hb1 = hb1_fwd + hb1_bwd
    hb1_labels = ['HB1_fwd', 'HB1_bwd']

    # Extract the data and export to a csv file
    hb1_df = [extract_data(hb1[i], data_label = hb1_labels[i], system=system) for i in range(len(hb1))]
    pd.concat(hb1_df).to_csv(Path(data_dir,system+"_hb1_data.csv"), index = False)

    # Export the fort.9 data
    fort9_df_hb1_fwd['branch_id'] = 'HB1_fwd'
    fort9_df_hb2_bwd['branch_id'] = 'HB1_bwd'
    fort9_df = pd.concat([fort9_df_hb1_fwd, fort9_df_hb2_bwd])
    fort9_df = fort9_df.drop(columns=['multipliers'])
    fort9_df['system_id'] = system
    fort9_df.to_csv(Path(data_dir,system+"_hb1_fort9.csv"), index = False)

    return (bd + hb1)

def main():

    # Change to the correct working directory containing the auto files
    auto_files_dir = "m2_upregulation_auto_files/"
    os.chdir(auto_files_dir)

    # Run P1 and P2 low
    print("*** Running P1 low... ***")
    bdP1L = run_1d_bifurcation("P1_low")

    print("*** Running P2 low... ***")
    bdP2L = run_1d_bifurcation("P2_low")

    # Run P1 and P2_high
    print("*** Running P1 high... ***")
    bdP1H = run_1d_bifurcation("P1_high")
        
    print("*** Running P2 high... ***")
    bdP2H = run_1d_bifurcation("P2_high")
        
    p = auto.plot(bdP1L + bdP2L + bdP1H + bdP2H)
    p.config(bifurcation_y=['EA'])
    auto.wait()

    # Clear the files
    auto.clean()


if __name__ == "__main__":
   main()