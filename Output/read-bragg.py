import numpy as np
import uproot4
import matplotlib.pyplot as plt
from ComptCamFunctions import get_hist_data

# --- physical constants ---
MEV_TO_JOULE = 1.602176634e-13
WATER_DENSITY = 1.0  # g/cm^3

# --- scoring geometry ---
SCORING_AREA_MM2 = 200.0 * 200.0  # 20x20 cm phantom cross-section

def energy_to_dose(values_mev, edges_mm, area_mm2, density_g_cm3=WATER_DENSITY):
    """
    Convert an energy-deposit histogram (MeV per bin, summed over all
    simulated primaries) into dose (Gy per bin).
    """
    bin_widths_mm = np.diff(edges_mm)                 # mm, will be 1.0 here
    volume_cm3 = area_mm2 * bin_widths_mm / 1000.0     # mm^3 -> cm^3
    mass_kg = density_g_cm3 * volume_cm3 / 1000.0      # g -> kg

    energy_j = values_mev * MEV_TO_JOULE
    dose_gy = energy_j / mass_kg
    return dose_gy


def range_finder(x, dose, percentage_falloff):
    peak_idx = np.argmax(dose)
    peak_dose = dose[peak_idx]
    dose_after_peak = dose[peak_idx:]
    x_after_peak = x[peak_idx:]
    range_dose = percentage_falloff * peak_dose
    closest_idx = np.argmin(np.abs(dose_after_peak - range_dose))
    proton_range = x_after_peak[closest_idx]
    return proton_range


energies = [120, 121, 122, 123]
plt.figure()
for e in energies:

    data = get_hist_data(f"pbt_merged/pbt_merged_120MeV_hists.root", ["BraggPeak"])
    values, edges = data["BraggPeak"]
    centers = 0.5 * (edges[:-1] + edges[1:]) + 100

    dose = energy_to_dose(values, edges, SCORING_AREA_MM2)

    R = range_finder(centers, dose, 0.9)
    plt.plot(centers, dose, label=f"{e} MeV, R = {R:.2f} mm")

plt.xlabel("Depth (mm)")
plt.ylabel("Dose (Gy)")
plt.legend()
plt.grid(alpha=0.5)
plt.show()