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
from ComptCamFunctions import (
    get_position_matrix,
    get_cones,
    simple_back_projection,
    stochastic_origin_ensemble,
    gaussian_prior,
    get_range_sigmoid,
    SOE_median_of_means,
    get_hist_data
)
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable

k = 1  # empirical
Xbins = 100
Ybins = 50
Zbins = 1  # only 2D for now
xmin = -100  # mm
xmax = 100  # mm
ymin = -50  # mm
ymax = 50  # mm
zmin = -200
zmax = -200
r = get_position_matrix(Xbins, Ybins, Zbins, xmin, xmax, ymin, ymax, zmin, zmax)
slice_num = int(Zbins / 2)
gauss = gaussian_prior(Xbins, Ybins, Zbins, sigma_y=0.5, sigma_z=0.5)
process_realistic = False
use_median_of_means = False
proton_energy = 120
n_events_to_move = 200
n_soe_iterations = 3000
output_subdir = "pbt_merged/"

# ------------------ Read ground truth histogram
rootfile_path = output_subdir + f"pbt_merged_{proton_energy}MeV_hists.root"
hist_names = {
    "BraggPeak": 1,  # 1D
    "Gamma2MeV": 2,  # 2D
    "Gamma4MeV": 2,
    "Gamma5MeV": 2,
    "Gamma6MeV": 2,
    "TotalGammas": 2,
}
hist_data = get_hist_data(rootfile_path, hist_names)
values, xedges, yedges = hist_data["TotalGammas"]

# ------------------- Read compt cam data
if process_realistic:
    file_path = output_subdir + f"realData_pbt_merged_{proton_energy}MeV_comptcam.csv"
    on = 'cm'
else:
    file_path = output_subdir + f"idealData_pbt_merged_{proton_energy}MeV_comptcam.csv"
    on = 'pos'

df = pd.read_csv(file_path)

# ---- Reconstructions -------
cones = get_cones(df, r, k, on=on)

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

if use_median_of_means:
    SOE = SOE_median_of_means(cones, Xbins, Ybins, Zbins, 10, n_events_to_move, n_soe_iterations, gauss)
else:
    SOE_chain = stochastic_origin_ensemble(cones, Xbins, Ybins, Zbins, n_events_to_move, n_soe_iterations, gauss)
    SOE = np.mean(SOE_chain, axis=0)

# plt.figure()
# plt.imshow(SOE[:, :, slice_num].T, cmap='inferno',
#                   extent=[xmin, xmax, ymin, ymax], origin='lower')
# cb = plt.colorbar(orientation='horizontal')
# cb.set_label(r'Average $N_\gamma$', labelpad=15)
# plt.xlabel("X (mm)")
# plt.ylabel("Y (mm)")
# plt.tight_layout()
# plt.title("SOE", fontweight='bold')
# plt.show()

# --------------------- Results ---------------------------------------
x = np.linspace(xmin + 100, xmax + 100, Xbins)

f, ax = plt.subplots(2, 2, figsize=(10, 8))

# --- Top row: images with tightly-fitted colorbars ---
im1 = ax[0, 0].imshow(
    values.T,
    extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
    origin='lower',
    cmap='inferno'
)
div1 = make_axes_locatable(ax[0, 0])
cax1 = div1.append_axes("right", size="5%", pad=0.05)
cb1 = f.colorbar(im1, cax=cax1)
cb1.set_label(r'$N_\gamma$', labelpad=15)
ax[0, 0].set_xlabel("X (mm)")
ax[0, 0].set_ylabel("Y (mm)")
ax[0, 0].set_title("Ground truth", fontweight='bold')

im2 = ax[0, 1].imshow(
    SOE[:, :, slice_num].T,
    cmap='inferno',
    extent=[xmin, xmax, ymin, ymax],
    origin='lower'
)
div2 = make_axes_locatable(ax[0, 1])
cax2 = div2.append_axes("right", size="5%", pad=0.05)
cb2 = f.colorbar(im2, cax=cax2)
cb2.set_label(r'Average $N_\gamma$', labelpad=15)
ax[0, 1].set_xlabel("X (mm)")
ax[0, 1].set_ylabel("Y (mm)")
ax[0, 1].set_title("SOE", fontweight='bold')

# --- Bottom row: line plots, with invisible spacer axes to match top-row width ---
ax[1, 0].plot(x, np.sum(values.T, axis=0))
ax[1, 0].set_xlabel("Depth (mm)")
ax[1, 0].set_ylabel(r'Summed value over Y in ROI')
div3 = make_axes_locatable(ax[1, 0])
cax3 = div3.append_axes("right", size="5%", pad=0.05)
cax3.axis("off")

ax[1, 1].plot(x, np.sum(SOE[:, 15:35, slice_num].T, axis=0))
ax[1, 1].set_xlabel("Depth (mm)")
ax[1, 1].set_ylabel(r'Summed value over Y in ROI')
div4 = make_axes_locatable(ax[1, 1])
cax4 = div4.append_axes("right", size="5%", pad=0.05)
cax4.axis("off")

f.tight_layout()
plt.show()
