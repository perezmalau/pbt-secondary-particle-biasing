"""
Last modified on Fri 20 Mar 2026
@author: Maria Perez-Lara, University College London, University of Warwick

Purpose: Execute the analysis and reconstruction of the Compton scattering data

We define the true data to be CZT camera output having:
- Ability to find coincident hits coming from the same gamma track
- Perfect energy and spatial resolution
- If drop_faulty is set to true, also the ability to drop backscatters

The realistic camera output therefore has:
- Coincident hits corresponding to the same event ID (but no track info)
- Finite spatial resolution (500 micron pitch)
- Discrimination from events producing more than one cluster in a single stage
"""
from ComptCamFunctions import (get_true_information,
                               get_detector_information,
                               get_compton_scatters,
                               get_position_matrix,
                               get_simple_backprojection
                               )
import matplotlib.pyplot as plt


rootfile = "biasing_test_5E7.root"
k = 1.5 # empirical
Xbins = 100
Ybins = 50
Zbins = 1
xmin = -100 # mm
xmax = 100 # mm
ymin = -50 # mm
ymax = 50 # mm
zref = -200 # mm
low_energy = 1 # MeV
high_energy = 7 # MeV
r = get_position_matrix(Xbins, Ybins, Zbins, xmin, xmax, ymin, ymax, zref, zref)

# --------- True data ---------
df_true = get_true_information(rootfile, drop_faulty=True)
df_true_cc = get_compton_scatters(df_true, on='pos').sort_values(by='EventID')
df_true_ene = df_true_cc[(df_true_cc['initEnergy'] < high_energy) & (df_true_cc['initEnergy'] > low_energy)]
SBP_t, df_sbp_t = get_simple_backprojection(df_true_ene, Xbins, Ybins, Zbins, r, k, zref, on='pos')

# ------- Realistic data -------
df_measured = get_detector_information(rootfile, Npix=120, pitch=0.5)
df_measured['Zcm_1'] = 2.5
df_measured['Zcm_2'] = 31.5
df_measured['Zcm_3'] = 60.5
df_measured_cc = get_compton_scatters(df_measured, on='cm').sort_values(by='EventID')
df_measured_ene = df_measured_cc[(df_measured_cc['initEnergy'] < high_energy) & (df_measured_cc['initEnergy'] > low_energy)]
SBP_m, df_sbp_m = get_simple_backprojection(df_measured_ene, Xbins, Ybins, Zbins, r, k, zref, on='pos')

# ------- Reconstructed images ---------
plt.figure()
plt.imshow(SBP_t[:, :, 0].T, cmap='gnuplot', extent=[-xmax, -xmin, ymin, ymax], origin='lower')
cb = plt.colorbar(orientation='horizontal')
cb.set_label(r'$N_\gamma$ (counts)', labelpad=15)
plt.xlabel("X (mm)")
plt.ylabel("Y (mm)")
plt.tight_layout()
plt.title("Ideal camera", fontweight='bold')
plt.show()

plt.figure()
plt.imshow(SBP_m[:, :, 0].T, cmap='gnuplot', extent=[-xmax, -xmin, ymin, ymax], origin='lower')
cb = plt.colorbar(orientation='horizontal')
cb.set_label(r'$N_\gamma$ (counts)', labelpad=15)
plt.xlabel("X (mm)")
plt.ylabel("Y (mm)")
plt.tight_layout()
plt.title("Realistic camera", fontweight='bold')
plt.show()
