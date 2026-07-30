"""
Purpose: Merge (sum) ComptCam data and ground truth histograms across 10 PBT simulation root files,
repeated for each proton beam energy (120, 121, 122, 123 MeV),
save the combined histograms to a new root file per energy,
and save the combined ComptCam data to a new csv file.
"""
import uproot4
import numpy as np
import pandas as pd
from ComptCamFunctions import (
                               merge_hists,
                               save_hists,
                               get_comptcam_data,
                               save_comptcam_data
                               )

# Names of the histograms to merge
HIST_NAMES = {
    "BraggPeak": 1,  # 1D
    "Gamma2MeV": 2,  # 2D
    "Gamma4MeV": 2,
    "Gamma5MeV": 2,
    "Gamma6MeV": 2,
    "TotalGammas": 2,
}

BEAM_ENERGIES = [120, 121, 122, 123]  # MeV
N_FILES = 10
N_PRIMARIES = "1E9"

# Loop over each proton beam energy
for energy in BEAM_ENERGIES:
    print(f"Processing {energy} MeV...")
    base_path = f"pbt_{energy}MeV_{N_PRIMARIES}"
    summed = merge_hists(base_path, HIST_NAMES, n_files=N_FILES)
    out_file = f"pbt_merged_{energy}MeV"
    save_hists(summed, out_file + "_hists.root")

    df_ideal, df_realistic = get_comptcam_data(base_path, low_energy=2, high_energy=7, n_files=N_FILES)
    save_comptcam_data(df_ideal, df_realistic, out_file + "_comptcam.csv")
