"""
Last modified on Mon 6 Jul 2026
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

# Sums histograms across multiple simulation files
def get_summed_hists_bnct(base_path, n_files=50):
    dose_sum = None
    g478_sum = None
    g2200_sum = None

    for i in range(1, n_files + 1):
        rootfile = f"{base_path}_{i}.root"
        try:
            d, g1, g2 = get_hists_bnct(rootfile)
        except Exception as e:
            print(f"Warning: could not read {rootfile} ({e}), skipping.")
            continue

        if dose_sum is None:
            dose_sum = np.zeros_like(d)
            g478_sum = np.zeros_like(g1)
            g2200_sum = np.zeros_like(g2)

        dose_sum += d
        g478_sum += g1
        g2200_sum += g2

    return dose_sum, g478_sum, g2200_sum

# Plots slices that pass through boronated tumour
# Specify title of the figure, quantity (e.g. Dose or N gammas) and unit (e.g. a.u. or Gy)
def plot_slices(arr, title, quantity, unit):
    f, ax = plt.subplots(1, 3, figsize=(12, 5))
    f.suptitle(title, fontsize=14, fontweight='bold')
    # Sagittal plane (YZ)
    im1 = ax[0].imshow(arr[16], cmap='gnuplot', aspect='auto', extent=[-85, 85, -100, 100], vmin=0, vmax=np.max(arr))
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

base_path = "bnct_test_5E10/bnct_test_1E9"
d, g1, g2 = get_summed_hists_bnct(base_path, n_files=50)
plot_slices(d, "Total dose (5E10 total)", "Dose", "Gy")
plot_slices(g1, "478 keV Gamma Distribution (5E10 total)", "N gammas", "a.u.")
plot_slices(g2, "2.2 MeV Gamma Distribution (5E10 total)", "N gammas", "a.u.")
