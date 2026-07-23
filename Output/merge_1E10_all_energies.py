"""
Purpose: Merge (sum) ComptCam data and ground truth histograms across 10 PBT simulation root files,
repeated for each proton beam energy (120, 121, 122, 123 MeV),
save the combined histograms to a new root file per energy,
and save the combined ComptCam data to a new csv file.
"""
import uproot4
import numpy as np
import pandas as pd
from ComptCamFunctions import (get_true_information,
                               get_detector_information,
                               get_compton_scatters,
                               )

# Names of the histograms to merge
HIST_NAMES = {
    "BraggPeak":    1,  # 1D
    "Gamma2MeV":    2,  # 2D
    "Gamma4MeV":    2,
    "Gamma5MeV":    2,
    "Gamma6MeV":    2,
    "TotalGammas":  2,
}

BEAM_ENERGIES = [120, 121, 122, 123]  # MeV
N_FILES = 10
N_PRIMARIES = "1E9"

def get_hist_data(rootfile, hist_names):
    """Reads to_numpy() output for each named histogram in one file."""
    file = uproot4.open(rootfile)
    data = {}
    for name in hist_names:
        h = file[name]
        data[name] = h.to_numpy()  # (values, edges...) tuple
    return data

def merge_hists(base_path, hist_names, n_files=N_FILES):
    """Sums each named histogram across n_files, keeping the first file's edges."""
    summed = {}   # name -> [values_sum, edge1, edge2, ...]

    for i in range(1, n_files + 1):
        rootfile = f"{base_path}_{i}.root"
        try:
            data = get_hist_data(rootfile, hist_names)
        except Exception as e:
            print(f"Warning: could not read {rootfile} ({e}), skipping.")
            continue

        for name, tup in data.items():
            vals, *edges = tup
            if name not in summed:
                summed[name] = [np.zeros_like(vals), *edges]  # init with edges from first file
            summed[name][0] += vals

    return summed

def get_comptcam_data(base_path, low_energy=2, high_energy=7, n_files=N_FILES):
    """Reads ComptCam data from all files into a single dataframe."""
    all_ideal_dfs = []
    all_realistic_dfs = []
    n_read = 0
    for i in range(1, n_files + 1):
        rootfile = f"{base_path}_{i}.root"
        try:
            # ------- Ideal data -------
            df_true = get_true_information(rootfile, drop_faulty=True, stages_hit=3)
            df_true_cc = get_compton_scatters(df_true, on='pos').sort_values(by='EventID')
            df_true_ene = df_true_cc[(df_true_cc['initEnergy'] < high_energy) & (df_true_cc['initEnergy'] > low_energy)]
            all_ideal_dfs.append(df_true_ene)

            # ------- Realistic data -------
            df_measured = get_detector_information(rootfile, Npix=120, pitch=0.5, stages_hit=3, pixlimit=5)
            df_measured['Zcm_1'] = 2.5
            df_measured['Zcm_2'] = 31.5
            df_measured['Zcm_3'] = 60.5
            df_measured_cc = get_compton_scatters(df_measured, on='cm').sort_values(by='EventID')
            df_measured_ene = df_measured_cc[(df_measured_cc['initEnergy'] < high_energy) & (df_measured_cc['initEnergy'] > low_energy)]
            all_realistic_dfs.append(df_measured_ene)
            n_read += 1
        except Exception as e:
            print(f"Warning: could not read ComptCam data from {rootfile} ({e}), skipping.")
            continue

    print(f"  Merged ComptCam data from {n_read}/{n_files} files.")
    combined_df_ideal = pd.concat(all_ideal_dfs, ignore_index=True) if all_ideal_dfs else pd.DataFrame()
    combined_df_realistic = pd.concat(all_realistic_dfs, ignore_index=True) if all_realistic_dfs else pd.DataFrame()

    return combined_df_ideal, combined_df_realistic

def save_hists(summed, out_file):
    """Writes summed histograms to a new root file using uproot4's write support."""
    with uproot4.recreate(out_file) as f:
        for name, tup in summed.items():
            vals, *edges = tup
            # (values, edges...) matches np.histogram / np.histogram2d convention,
            # which uproot4 recognizes and writes as TH1D/TH2D.
            f[name] = tuple([vals, *edges])
    print(f"  Wrote {len(summed)} histograms to {out_file}")

def save_comptcam_data(df_ideal, df_realistic, out_file):
    """Writes ComptCam data to a new csv file."""
    df_ideal.to_csv("idealData_" + out_file, index=False)
    df_realistic.to_csv("realData_" + out_file, index=False)
    print(f"  Wrote ComptCam data to {out_file}")

# Loop over each proton beam energy
for energy in BEAM_ENERGIES:
    print(f"Processing {energy} MeV...")
    base_path = f"pbt_{energy}MeV_{N_PRIMARIES}"
    summed = merge_hists(base_path, HIST_NAMES, n_files=N_FILES)
    out_file = f"pbt_merged_{energy}MeV"
    save_hists(summed, out_file+"_hists.root")

    df_ideal, df_realistic = get_comptcam_data(base_path, low_energy=2, high_energy=7, n_files=N_FILES)
    save_comptcam_data(df_ideal, df_realistic, out_file+"_comptcam.csv")