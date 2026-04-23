"""
Last modified on Fri 20 Mar 2026
@author: Maria Perez-Lara, University College London, University of Warwick

Purpose: Execute the analysis and reconstruction of the Compton scattering data

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
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable


rootfile = "biasing_120MeV_1E8.root"
k = 1.5 # empirical
Xbins = 100
Ybins = 50
Zbins = 50
xmin = -100 # mm
xmax = 100 # mm
ymin = -50  # mm
ymax = 50  # mm
zmin = -250
zmax = -150
low_energy = 1 # MeV
high_energy = 7 # MeV
r = get_position_matrix(Xbins, Ybins, Zbins, xmin, xmax, ymin, ymax, zmin, zmax)
slice_num = int(Zbins / 2)
# --------- True data ---------
df_true = get_true_information(rootfile, drop_faulty=True, stages_hit=2)
df_true_cc = get_compton_scatters(df_true, on='pos').sort_values(by='EventID')
df_true_ene = df_true_cc[(df_true_cc['initEnergy'] < high_energy) & (df_true_cc['initEnergy'] > low_energy)]

gauss = gaussian_prior(Xbins, Ybins, Zbins, sigma_y=0.5, sigma_z=0.5)
# plt.figure()
# plt.imshow(gauss[:, :, int(Zbins/2)].T, cmap='hot', extent=[-xmax, -xmin, ymin, ymax], origin='lower')
# cb = plt.colorbar()
# cb.set_label('P (a.u.)', labelpad=15)
# plt.title(r'Gaussian prior ($\bf{\sigma_y}$ = 0.5, $\bf{\sigma_z}$ = 0.5)', fontweight='bold')
# plt.xlabel("X (mm)")
# plt.ylabel("Y (mm)")
# plt.tight_layout()
# plt.show()

# # ----SBP-------
cones = get_cones(df_true_ene, r, k, on='pos')
# SBP = simple_back_projection(Xbins, Ybins, Zbins, cones)
#
# plt.figure()
# plt.imshow(SBP[:, :, slice_num].T, cmap='gnuplot', extent=[-xmax, -xmin, ymin, ymax], origin='lower')
# cb = plt.colorbar()
# cb.set_label(r'$N_\gamma$ (counts)', labelpad=15)
# plt.xlabel("X (mm)")
# plt.ylabel("Y (mm)")
# plt.tight_layout()
# plt.title("Ideal camera - SBP", fontweight='bold')
# plt.show()

# ----SOE-------
SOE_chain, probs = stochastic_origin_ensemble(cones,
                                              Xbins, Ybins, Zbins,
                                              N_events=200,
                                              N_soe=3000,
                                              weights=gauss,
                                              percent_convergence=3
                                              )

f, ax = plt.subplots(2,2)
im1 = ax[0,0].imshow(SOE_chain[0][:,:,slice_num].T, cmap='gnuplot', extent=[-xmax, -xmin, ymin, ymax], origin='lower')
ax[0,0].set_title("$Y_0$")
divider1 = make_axes_locatable(ax[0,0])
cax1 = divider1.append_axes("right", size="5%", pad=0.05)
cb1 = plt.colorbar(im1, cax=cax1, orientation='vertical')
im2 = ax[1].imshow(SOE_chain[-1][:, :, slice_num].T, cmap='gnuplot', extent=[-xmax, -xmin, ymin, ymax], origin='lower')
ax[1].set_title("$Y_{N_{itr}}$")
divider2 = make_axes_locatable(ax[1])
cax2 = divider2.append_axes("right", size="5%", pad=0.05)
cb2 = plt.colorbar(im2, cax=cax2, orientation='vertical')
#plt.subplots_adjust(wspace=0.3)
plt.suptitle(r"SOE ($\bf{N_{itr}}$=1000, $\bf{N_{evt}}$=200, Gaussian prior)", fontweight='bold')
plt.tight_layout()
plt.show()
# final_state = SOE_chain[-1]
# plt.figure()
# plt.imshow(final_state[:, :, slice_num].T, cmap='gnuplot', extent=[-xmax, -xmin, ymin, ymax], origin='lower')
# cb = plt.colorbar(orientation='horizontal')
# cb.set_label(r'$N_\gamma$ (counts)', labelpad=15)
# plt.xlabel("X (mm)")
# plt.ylabel("Y (mm)")
# plt.tight_layout()
# plt.title("Ideal camera - SOE", fontweight='bold')
# plt.show()


# ------- Realistic data -------
# df_measured = get_detector_information(rootfile, Npix=120, pitch=0.5, stages_hit=3, pixlimit=5)
# df_measured['Zcm_1'] = 2.5
# df_measured['Zcm_2'] = 31.5
# df_measured['Zcm_3'] = 60.5
# df_measured_cc = get_compton_scatters(df_measured, on='cm').sort_values(by='EventID')
# df_measured_ene = df_measured_cc[(df_measured_cc['initEnergy'] < high_energy) & (df_measured_cc['initEnergy'] > low_energy)]
# SBP_m, df_sbp_m = get_simple_backprojection(df_measured_ene, Xbins, Ybins, Zbins, r, k, zref, on='cm')
#
# plt.figure()
# plt.imshow(SBP_m[:, :, 0].T, cmap='gnuplot', extent=[-xmax, -xmin, ymin, ymax], origin='lower')
# cb = plt.colorbar(orientation='horizontal')
# cb.set_label(r'$N_\gamma$ (counts)', labelpad=15)
# plt.xlabel("X (mm)")
# plt.ylabel("Y (mm)")
# plt.tight_layout()
# plt.title("Realistic camera - SBP", fontweight='bold')
# plt.show()
