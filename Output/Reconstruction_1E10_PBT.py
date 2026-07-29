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
                               gaussian_prior, get_range_sigmoid
                               )
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os

k = 1 # empirical
Xbins = 100
Ybins = 50
Zbins = 1 # only 2D for now
xmin = -100 # mm
xmax = 100 # mm
ymin = -50  # mm
ymax = 50  # mm
zmin = -200
zmax = -200
low_energy = 1 # MeV
high_energy = 7 # MeV
r = get_position_matrix(Xbins, Ybins, Zbins, xmin, xmax, ymin, ymax, zmin, zmax)
slice_num = int(Zbins / 2)
gauss = gaussian_prior(Xbins, Ybins, Zbins, sigma_y=0.5, sigma_z=0.5)
process_realistic = True
proton_energy = 120

# -------------------------------------------------------------

if process_realistic:
    file_path = f"pt_test_120MeV_1E10/fullPBT_{proton_energy}MeV_1E10_realisticComptCam_{low_energy}-{high_energy}MeV.csv"
    on = 'cm'
else:
    file_path = f"pt_test_120MeV_1E10/fullPBT_{proton_energy}MeV_1E10_idealComptCam_{low_energy}-{high_energy}MeV.csv"
    on = 'pos'

if os.path.isfile(file_path):
    combined_df = pd.read_csv(file_path)
else:
    all_dfs = []
    for i in range(1, 11):
        rootfile = f"pt_test_{proton_energy}MeV_1E10/pt_test_{proton_energy}MeV_1E9_{i}.root"
        if process_realistic:
            df = get_detector_information(rootfile, Npix=120, pitch=0.5, stages_hit=3, pixlimit=5)
            df['Zcm_1'] = 2.5
            df['Zcm_2'] = 31.5
            df['Zcm_3'] = 60.5
            df_cc = get_compton_scatters(df, on='cm').sort_values(by='EventID')
            df_ene = df_cc[(df_cc['initEnergy'] < high_energy) & (df_cc['initEnergy'] > low_energy)]
            all_dfs.append(df_ene)
        else:
            df = get_true_information(rootfile, drop_faulty=True, stages_hit=3)
            df_cc = get_compton_scatters(df, on='pos').sort_values(by='EventID')
            df_ene = df_cc[(df_cc['initEnergy'] < high_energy) & (df_cc['initEnergy'] > low_energy)]

    combined_df = pd.concat(all_dfs, ignore_index=True)
    combined_df.to_csv(file_path, index=False)


# ---- Reconstructions -------
cones = get_cones(combined_df, r, k, on=on)

# SBP = simple_back_projection(Xbins, Ybins, Zbins, cones)
#
# plt.figure()
# plt.imshow(SBP[:, :, slice_num].T, cmap='gnuplot', extent=[-xmax, -xmin, ymin, ymax], origin='lower')
# cb = plt.colorbar(orientation='horizontal')
# cb.set_label(r'$N_\gamma$ (counts)', labelpad=15)
# plt.xlabel("X (mm)")
# plt.ylabel("Y (mm)")
# plt.tight_layout()
# plt.title("SBP", fontweight='bold')
# plt.show()

all_chains = []
for i in range(1):
    SOE_chain = stochastic_origin_ensemble(cones,
                                                  Xbins, Ybins, Zbins,
                                                  N_events=200,
                                                  N_soe=20000,
                                                  weights=gauss,
                                                  #percent_convergence=3,
                                                  #alpha=0.001
                                                  )
    SOE = np.mean(SOE_chain, axis=0) # or to look at the last one only, use SOE_chain[-1]
    all_chains.append(SOE)

MoM_SOE = np.median(all_chains, axis=0)


x = np.linspace(xmin, xmax, Xbins)
y = np.sum(MoM_SOE[:, 15:35, slice_num].T, axis=0)

f, ax = plt.subplots(2, 1)

im = ax[0].imshow(MoM_SOE[:, :, slice_num].T, cmap='inferno',
                  extent=[xmin, xmax, ymin, ymax], origin='lower')
cb = f.colorbar(im, ax=ax[0])
cb.set_label(r'Average $N_\gamma$', labelpad=15)
ax[0].set_xlabel("X (mm)")
ax[0].set_ylabel("Y (mm)")

ax[1].plot(x, y)
ax[1].set_xlabel("X (mm)")
ax[1].set_ylabel(r'Summed value over Y in ROI')

f.tight_layout()
plt.show()
