from ComptCamFunctions import (get_true_information,
                               get_detector_information,
                               get_triple_scatters,
                               get_position_matrix,
                               get_conical_surface,
                               get_real_dataframe,
                               get_simple_backprojection
                               )
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


rootfile = "biasing_test_2E7.root"
df_true = get_true_information(rootfile, drop_double_coincidences=False)
df_cc = get_triple_scatters(df_true, on='pos')

k=1
Xbins = 100
Ybins = 100
Zbins = 1
xmin = -100
xmax = 100
ymin = -100
ymax = 100
zref = -200

r = get_position_matrix(Xbins, Ybins, Zbins, xmin, xmax, ymin, ymax, zref, zref)
df_ene = df_cc[(df_cc['initEnergy'] < 6.1) & (df_cc['initEnergy'] > 2.1)]
SBP, df_sbp = get_simple_backprojection(df_ene, Xbins, Ybins, Zbins, r, k, zref, on='pos')

plt.figure()
plt.imshow(np.rot90(SBP[:, :, 0], k=3), cmap='gnuplot', extent=[-xmax, -xmin, ymin, ymax])
cb = plt.colorbar()
cb.set_label(r'$N_\gamma$ (counts)', labelpad=15)
plt.xlabel("X (mm)")
plt.ylabel("Y (mm)")
plt.tight_layout()
plt.show()

# df_measured = get_detector_information(rootfile, Npix=120, pitch=0.5)
# z1 = np.nanmean(df_true['PosZ_1'].values)
# z2 = np.nanmean(df_true['PosZ_2'].values)
# z3 = np.nanmean(df_true['PosZ_3'].values)
# df_measured['Zcm_1'] = z1
# df_measured['Zcm_2'] = z2
# df_measured['Zcm_3'] = z3
# df_cc = get_triple_scatters(df_measured)
# df_merged = pd.merge(df_cc, df_true, on='EventID', how="left")
#
