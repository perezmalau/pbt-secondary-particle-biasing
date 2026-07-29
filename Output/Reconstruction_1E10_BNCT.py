"""
Last modified on Wed 24 Jun 2026
@author: Maria Perez-Lara, University College London, University of Warwick

Purpose: Execute the analysis and reconstruction of the Compton scattering data with high stats (10 sims 1E9 each)

We define the true data to be CZT camera output having:
- Ability to find coincident hits coming from the same gamma track
- Perfect energy and spatial resolution
- If drop_faulty is set to true, also the ability to drop backscatters (e.g. compt-phot-compt)

The realistic camera output therefore has:
- Coincident hits corresponding to the same event ID (but no individual track or particle info)
- Finite spatial resolution (500 micron pitch, no depth resolution)
- Discrimination from events producing more than one cluster in a single stage

NOTE: The number of faulty events due to backscattering drops significantly when you use solely triple coincidences.
For the future: either only use 3-stage coincidences or have time resolution between stages to build a corrected cone.
"""
from ComptCamFunctions import (get_true_information,
                               get_detector_information,
                               get_compton_scatters,
                               get_position_matrix,
                               get_cones,
                               simple_back_projection, stochastic_origin_ensemble,
                               gaussian_prior
                               )
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

k = 1 # empirical
Xbins = 70
Ybins = 100
Zbins = 1 # only 2D for now
xmin = -70 # mm
xmax = 70 # mm
ymin = -100  # mm
ymax = 100  # mm
zmin = -200
zmax = -200
low_energy = 0.46 # MeV
high_energy = 0.49 # MeV
r = get_position_matrix(Xbins, Ybins, Zbins, xmin, xmax, ymin, ymax, zmin, zmax)
slice_num = int(Zbins / 2)
gauss = gaussian_prior(Xbins, Ybins, Zbins, sigma_y=0.5, sigma_z=0.5)

all_ideal_dfs = []
all_realistic_dfs = []
for i in range(1, 51):
    rootfile = f"bnct_test_5E10/bnct_test_1E9_{i}.root"

    # --------- True data ---------
    # df_true = get_true_information(rootfile, drop_faulty=True, stages_hit=2)
    # df_true_cc = get_compton_scatters(df_true, on='pos').sort_values(by='EventID')
    # df_true_ene = df_true_cc[(df_true_cc['initEnergy'] < high_energy) & (df_true_cc['initEnergy'] > low_energy)]
    # all_ideal_dfs.append(df_true_ene)
    # ------- Realistic data -------
    df_measured = get_detector_information(rootfile, Npix=240, pitch=0.25, stages_hit=2, pixlimit=5)
    df_measured['Zcm_1'] = 2.5
    df_measured['Zcm_2'] = 31.5
    df_measured['Zcm_3'] = 60.5
    df_measured_cc = get_compton_scatters(df_measured, on='cm').sort_values(by='EventID')
    df_measured_ene = df_measured_cc[(df_measured_cc['initEnergy'] < high_energy) & (df_measured_cc['initEnergy'] > low_energy)]
    all_realistic_dfs.append(df_measured_ene)

#combined_df_ideal = pd.concat(all_ideal_dfs, ignore_index=True)
combined_df_realistic = pd.concat(all_realistic_dfs, ignore_index=True)

#combined_df_ideal.to_csv("bnct_test_5E10/fullBNCT_5E10_idealComptCam_478keV.csv", index=False)
combined_df_realistic.to_csv("bnct_test_5E10/fullBNCT_5E10_realisticComptCam_478keV.csv", index=False)

#cones_ideal = get_cones(combined_df_ideal, r, k, on='pos')
cones_realistic = get_cones(combined_df_realistic, r, k, on='cm')

# ---- Reconstructions -------
SBP = simple_back_projection(Xbins, Ybins, Zbins, cones_realistic)

plt.figure()
plt.imshow(SBP[:, :, slice_num].T, cmap='gnuplot', extent=[-xmax, -xmin, ymin, ymax], origin='lower')
cb = plt.colorbar()
cb.set_label(r'$N_\gamma$ (counts)', labelpad=15)
plt.xlabel("X (mm)")
plt.ylabel("Y (mm)")
plt.tight_layout()
plt.title("Ideal camera - SBP", fontweight='bold')
plt.show()

SOE_chain, probs = stochastic_origin_ensemble(cones_realistic,
                                              Xbins, Ybins, Zbins,
                                              N_events=200,
                                              N_soe=3000,
                                              #weights=gauss,
                                              #percent_convergence=3,
                                              #alpha=0.001
                                              )
SOE = np.mean(SOE_chain, axis=0) # or to look at the last one only, use SOE_chain[-1]

plt.figure()
plt.imshow(SOE[:, :, slice_num].T, cmap='inferno', extent=[-xmax, -xmin, ymin, ymax], origin='lower')
cb = plt.colorbar(orientation='horizontal')
cb.set_label(r'$N_\gamma$ (counts)', labelpad=15)
plt.xlabel("X (mm)")
plt.ylabel("Y (mm)")
plt.tight_layout()
plt.title("Realistic camera - SOE", fontweight='bold')
plt.show()
