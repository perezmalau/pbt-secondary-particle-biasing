"""
Last modified on Thu 2 Jul 2026
@author: Maria Perez-Lara, University College London, University of Warwick

Purpose: Obtain dose and gamma production histograms for the ground truth data, which are saved in the root outputs.
"""
import uproot4
import matplotlib.pyplot as plt
import numpy as np

# Loads histogram data from root output file
# Voxel sizes are 2 mm
def get_hists_bnct(rootfile):
    file = uproot4.open(rootfile)
    h_dose = file["Dose3D"]
    h_478 = file["Gamma478"]
    h_2200 = file["Gamma2200"]

    # to_numpy() gives (values, x_edges, y_edges, z_edges)
    dose_vals, xedges, yedges, zedges = h_dose.to_numpy()
    gamma478_vals, _, _, _ = h_478.to_numpy()
    gamma2200_vals, _, _, _ = h_2200.to_numpy()

    return dose_vals, gamma478_vals, gamma2200_vals

# Plots slices that pass through boronated tumour
# Specify title of the figure, quantity (e.g. Dose or N gammas) and unit (e.g. a.u. or Gy)
def plot_slices(arr, title, quantity, unit):
    f, ax = plt.subplots(1, 3, figsize=(12, 5))
    f.suptitle(title, fontsize=14, fontweight='bold')
    # Sagittal plane (YZ)
    im1 = ax[0].imshow(arr[18], cmap='gnuplot', aspect='auto', extent=[-85, 85, -100, 100], vmin=0, vmax=np.max(arr))
    ax[0].set_title("Sagittal plane", fontsize=12)
    ax[0].set_aspect('equal', adjustable='box')
    ax[0].set_xlabel("Z (mm)")
    ax[0].set_ylabel("Y (mm)")
    # Coronal plane (XZ)
    im2 = ax[1].imshow(np.rot90(arr, axes=(1, 0))[50], cmap='gnuplot', aspect='auto', extent=[-85, 85, -70, 70],
                       vmin=0, vmax=np.max(arr))
    ax[1].set_title("Coronal plane", fontsize=12)
    ax[1].set_aspect('equal', adjustable='box')
    ax[1].set_xlabel("Z (mm)")
    ax[1].set_ylabel("X (mm)")
    # Axial plane (XY)
    im3 = ax[2].imshow(np.rot90(arr, axes=(0, 2))[37], cmap='gnuplot', aspect='auto', extent=[-70, 70, -100, 100],
                       vmin=0, vmax=np.max(arr))
    cb = plt.colorbar(im3, ax=ax.ravel().tolist(), orientation='horizontal', shrink=0.8)
    cb.set_label(quantity+" ("+unit+")")
    ax[2].set_title("Axial plane", fontsize=12)
    ax[2].set_aspect('equal', adjustable='box')
    ax[2].set_xlabel("X (mm)")
    ax[2].set_ylabel("Y (mm)")
    plt.show()

d, g1, g2 = get_hists_bnct("bnct_test_1E7.root")
plot_slices(d, "Dose distribution", "Dose", "Gy")