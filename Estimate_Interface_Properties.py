import argparse
import numpy as np
import pandas as pd
from Generate_Connectivity import *
from CCRE import calc_CCRE
from TRE import calc_TRE

boundary = {"CCRE": 0.0572991533003641, "TRE": 0.03231143068649855} # (min_phys+max_chem)/2 = (anc+ttc)/2
elecPE_from_CCREPE = {"m_c": 194.9028879086265, "n_c": -9.977871123415241}
elecPE_from_TREPE = {"m_c": 408.884262157611, "n_c": -12.127302374271972}
elecPE_from_TAPE = {"CCRE": elecPE_from_CCREPE, "TRE": elecPE_from_TREPE}
prep_from_elec = {"m_c": -0.1883204174466927, "n_c": 46.12762704240741}
disp_from_Atoms = {"m_p": -9.592982832618027, "n_p": -3.246660944205928,
                   "m_p*": -9.648233292831103, "n_p*": -0.199992709599087,
                   "m_c": -10.793432432432434, "n_c": 8.268567567567686}
elec_from_disp = {"m_p": -0.12158262530535541, "n_p": -6.765094069822705,
                  "m_p*": -0.12153499860935289, "n_p*": -6.648080124226723}
Dads_from_elec = {"m_p": 0.006773379411837688, "n_p": 2.575442438617571,
                  "m_p*": 0.006493091894837407, "n_p*": 2.5985700158562475,
                  "m_c": 0.0012623590585077085, "n_c": 2.425472864012389}
θbend_from_Dads = {"m_p": -10.694521327961631, "n_p": 30.528310836923023,
                   "m_p*": -9.468338480064869, "n_p*": 27.119180385036653,
                   "m_c": -31.034604913688945, "n_c": 81.93414783704351}
dq_from_int = {"m_p": 0.00032359660226307394, "n_p": 0.053797801762578326,
               "m_p*": 0.00032500165488610834, "n_p*": 0.05952360996769873,
               "m_c": 0.0011726722987925175, "n_c": 0.10112937736742111}
dΦ_from_disp = {"m_p": 0.0011280709524252785, "n_p": -0.15455878707300463,
                "m_p*": 0.0011280709524252785, "n_p*": -0.15455878707300463,
                "m_c": 0.0011280709524252785, "n_c": -0.15455878707300463}
                
def Estimate_Interface_Properties(con_mat, methods, include_anc_in_fit):
    if include_anc_in_fit == True: add = ""
    else: add = "*"
    PE = len(con_mat)
    Atoms = PE + (PE - np.sum([1 for row in con_mat if np.sum(row) == 3]))
    TAPE = {}
    adsorption_regime = {}
    disp = {}
    elecPE = {}
    elec = {}
    prep = {}
    int = {}
    ads = {}
    Dads = {}
    θbend = {}
    dq = {}
    dΦ = {}
    for method in methods:
        if method == "CCRE":
            print()
            print("################################################")
            print(" Results of Conjugated Circuit Resonance Energy ")
            print("################################################")
            print()
            CCRE, CCREPE = calc_CCRE(con_mat)
            print(f"CCRE   = {np.round(CCRE,4)} β")
            print(f"CCREPE = {np.round(CCREPE,5)} β")
            TAPE[method] = CCREPE
        elif method == "TRE":
            print()
            print("################################################")
            print("    Results of Topological Resonance Energy     ")
            print("################################################")
            print()
            TRE, TREPE = calc_TRE(con_mat)
            print(f"TRE   = {np.round(float(TRE),4)} β")
            print(f"TREPE = {np.round(float(TREPE),5)} β")
            TAPE[method] = TREPE
        else:
            raise SystemExit("Error, your specified method is unknown!")
            
        if TAPE[method] < boundary[method]:
            adsorption_regime[method] = "chemisorption"
            m = "m_c"; n = "n_c"
        else:
            adsorption_regime[method] = "physisorption"
            m = "m_p"+add; n = "n_p"+add
    
        disp[method] = disp_from_Atoms[m] * Atoms + disp_from_Atoms[n]
        if adsorption_regime[method] == "chemisorption":
            elecPE[method] = elecPE_from_TAPE[method][m] * TAPE[method] + elecPE_from_TAPE[method][n]
            elec[method] = elecPE[method] * PE
            if elec[method] > 0:
                print()
                print("### WARNING ###")
                print(f"Although {method} predicted chemisorption, there is a positive electronic interaction with the surface!")
            prep[method] = prep_from_elec[m] * elec[method] + prep_from_elec[n]
        else:
            elec[method] = elec_from_disp[m] * disp[method] + elec_from_disp[n]
            if include_anc_in_fit == True: prep[method] = 6.2
            else: prep[method] = 5.4
    
        int[method] = elec[method] + disp[method]
        ads[method] = prep[method] + elec[method] + disp[method]
        Dads[method] = Dads_from_elec[m] * elec[method] + Dads_from_elec[n]
        θbend[method] = θbend_from_Dads[m] * Dads[method] + θbend_from_Dads[n]
        dq[method] = dq_from_int[m] * int[method] + dq_from_int[n]
        dΦ[method] = dΦ_from_disp[m] * disp[method] + dΦ_from_disp[n]

    if "TRE" in methods and "CCRE" in methods:
        if adsorption_regime["TRE"] != adsorption_regime["CCRE"]:
            print()
            print("### WARNING ###")
            print("The two topological aromaticity methods predict different adsorption regimes!")
            print("CCRE:", adsorption_regime["CCRE"])
            print("TRE: ", adsorption_regime["TRE"])

    return adsorption_regime, disp, elec, prep, int, ads, Dads, θbend, dq, dΦ
    
def print_results_table(methods, adsorption_regime, disp, elec, prep, int, ads, Dads, θbend, dq, dΦ):
    if_props = ["adsorption regime",
            "E_ads", "  E_prep", "  E_int", "    E_int(elec)", "    E_int(disp)",
            "D_ads", "θ_bend", "Δq", "ΔΦ"]
    units = ["", "[kJ/mol]", "[kJ/mol]", "[kJ/mol]", "[kJ/mol]", "[kJ/mol]",
             "[Å]", "[°]", "[e]", "[eV]"]
    results = pd.DataFrame()
    results["Interface Properties"] = if_props
    results.set_index("Interface Properties",inplace=True)
    results[""] = units
    for method in methods:
        results[method] = [adsorption_regime[method],
                           np.round(float(ads[method]),1), np.round(float(prep[method]),1), np.round(float(int[method]),1),
                           np.round(float(elec[method]),1), np.round(float(disp[method]),1),
                           np.round(float(Dads[method]),3), np.round(float(θbend[method]),2),
                           np.round(float(dq[method]),3), np.round(float(dΦ[method]),4)]
    print()
    print("################################################")
    print("      Estimated Interface Properties are:       ")
    print("################################################")
    print()
    print(results)
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Estimation of Cu(111)-Organic Interface Properties using Topological Aromaticity Indices derived solely from the Structure of the Molecule')
    parser.add_argument('Input', help='Structure Input from .xyz-Filename OR SMILES-String', type=str, default='input.xyz')
    parser.add_argument('--method', help='Topological Aromaticity Index which is used', nargs="+", default=["CCRE", "TRE"])
    parser.add_argument('--include_anc_in_fit', help='Include Anthracene (anc) into the fitted values of the data set', action='store_true', default=False)
    args = parser.parse_args()
    
    if args.Input.split(".")[-1] == "xyz": con_mat = con_mat_from_struc(args.Input)
    elif args.Input[0] == "[" and args.Input[-1] == "]": con_mat = con_mat_from_mat(args.Input)
    else: con_mat = con_mat_from_SMILES(args.Input)
    adsorption_regime, disp, elec, prep, int, ads, Dads, θbend, dq, dΦ = Estimate_Interface_Properties(con_mat, args.method, args.include_anc_in_fit)
    print_results_table(args.method, adsorption_regime, disp, elec, prep, int, ads, Dads, θbend, dq, dΦ)
