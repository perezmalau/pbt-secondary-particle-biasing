import uproot4
import numpy as np
import matplotlib.pyplot as plt


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
    file = uproot4.open(f"pbt_{e}MeV_1E6.root")
    hist = file["BraggPeak"]
    values, edges = hist.to_numpy()
    centers = 0.5 * (edges[:-1] + edges[1:]) + 100
    R = range_finder(centers, values, 0.9)
    plt.plot(centers, values, label=f"{e} MeV, R = {R:.2f} mm")
plt.xlabel("Depth (mm)")
plt.ylabel("Dose (Gy)")
plt.legend()
plt.grid(alpha=0.5)
plt.show()