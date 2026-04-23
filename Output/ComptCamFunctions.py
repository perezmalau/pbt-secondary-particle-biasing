"""
Last modified on Fri 20 Mar 2026
@author: Maria Perez-Lara, University College London, University of Warwick

Purpose: Read out 2 stage and 3 stage compton cam data,
taking into account different geometries for scatterer and absorber
"""

import pandas as pd
import numpy as np
import progressbar
import random
import time
import uproot4
from sklearn.cluster import DBSCAN
# import sys

# ------------------------------------------------------
# ------------ PREPROCESSING FUNCTIONS -----------------
# ------------------------------------------------------

# careful after implementing this!! Works for the detector that is facing the source (absorber)
# As the scatterer is backwards, x axis is mirrored, then change the flip flag to true
# Needs N pixels and pitch of the detector in mm
# Default values based on 3x3 HEXITEC detector array (40x40 pixels and 500 micron pitch each)
def pixel_to_coordinate(i, j, Npix, pitch, flip=False):
    if flip:
        x = pitch * ((Npix-1)/2 - i)
    else:
        x = pitch * (i - (Npix-1)/2)
    y = pitch * (j - (Npix-1)/2)
    return x, y

# Perform a DBSCAN clustering algorithm in a given hexitec file
def cluster_df(df_hex, on='EventID'):
    df = df_hex.copy()
    scale = np.array([0.5, 0.5, 2.0])
    X = scale * df[['Xpix', 'Ypix', on]].values
    scan = DBSCAN(eps=1.5, min_samples=1).fit(X)
    df['ClusterID'] = scan.labels_
    df_clustered = []
    groups = df.groupby('ClusterID')
    # print("Clustering...")
    # bar = progressbar.ProgressBar(maxval=len(groups), widgets=[progressbar.Bar('=', '[', ']'), ' ',
    #                                                            progressbar.Percentage()])
    # bar.start()

    for i, (_, group) in enumerate(groups):
        Edep_group = group['Edep'].sum()
        weighted_coords = np.average(group[['Xpix', 'Ypix']].values, axis=0, weights=group['Edep'])
        Xpix_cm, Ypix_cm = np.round(weighted_coords).astype(int)
        group['Xcm'] = Xpix_cm
        group['Ycm'] = Ypix_cm
        group['Etot'] = Edep_group
        df_clustered.append(group)
    #     bar.update(i + 1)
    # bar.finish()
    df_clustered = pd.concat(df_clustered, ignore_index=True)
    # df_clustered = df_clustered.drop(columns=['ClusterID'])
    return df_clustered


# Filter only events that show E1+E2 = Esource
def energy_filter_mono(df1, df2, Esource, epsilon, category='EventID'):
    df1_first = df1.groupby(category).first()
    df2_first = df2.groupby(category).first()
    E0 = df1_first['Etot'] + df2_first['Etot']
    isPrompt = ((Esource - epsilon <= E0) & (E0 <= Esource + epsilon))
    useful_events = isPrompt[isPrompt].index.values
    filtered_hex1 = df1[df1[category].isin(useful_events)]
    filtered_hex2 = df2[df2[category].isin(useful_events)]
    return filtered_hex1, filtered_hex2


# Energy resolution function, returns resolution as a percentage of energy
def resolution(a, b, E):
    return a + (b*np.sqrt(E))

# Blurring function to account for energy resolution in geant4 data
# For S017: a=-4.348, b=1.097
#For S018: a=-2.836, b=0.853
def blur_energies(df_merged):
    dfb = df_merged.copy()
    energies1 = dfb['Etot_1'].values*1000
    resolutions1 = resolution(-4.348, 1.097, energies1)
    resolutions1[resolutions1 < 0] = 0
    blurred_energies1 = np.random.normal(energies1, resolutions1/2.355)
    energies2 = dfb['Etot_2'].values*1000
    resolutions2 = resolution(-2.836, 0.853, energies2)
    resolutions2[resolutions2 < 0] = 0
    blurred_energies2 = np.random.normal(energies2, resolutions2/2.355)
    dfb['Etot_1'] = blurred_energies1/1000
    dfb['Etot_2'] = blurred_energies2/ 1000
    return dfb

# Only keep events that have cluster sizes of less than pixlimit and only show a single cluster per event
def filter_by_clusters(df, pixlimit=10, groupby='EventID'):
    cluster_sizes = df.groupby(groupby)['ClusterID'].count()
    nclusters = df.groupby(groupby)['ClusterID'].nunique()
    eventInfo = pd.DataFrame({'ClusterSize': cluster_sizes, 'Nclusters': nclusters})
    keeper_events = eventInfo[(eventInfo['ClusterSize'] < pixlimit) & (eventInfo['Nclusters'] == 1)].index
    df_filt = df[df[groupby].isin(keeper_events)]
    df_first = df_filt.groupby('EventID').first()[['Xcm', 'Ycm', 'Etot']]
    return df_first


# ------------------------------------------------------
# ------------ 3-STAGE CAMERA FUNCTIONS ----------------
# ------------------------------------------------------

# For ideal data

# Remove events with more than one hit in ideal data
def filter_single_hits(df):   # Discard multiple interactions in a single laye
    counts = df[['EventID', 'TrackID']].value_counts()
    single = counts[counts == 1].index
    return df.set_index(['EventID', 'TrackID']).loc[single].reset_index()

# Label an event as genuine double/triple coincidence or faulty
# A faulty event in a 3-stage coincidence could be due to backscatter (e.g. compt-phot-compt)
# A faulty event in a 2-stage coincidence could be due to double compton or backscatter (phot-compt)
def classify_event(row):
    has_1 = pd.notna(row.get('Xpos_1'))
    has_2 = pd.notna(row.get('Xpos_2'))
    has_3 = pd.notna(row.get('Xpos_3'))

    stages_hit = has_1 + has_2 + has_3

    if stages_hit == 3:
        # If the first stages show absorption, there was a backscatter, not useful
        if row['Process_1'] == 'phot' or row['Process_2'] == 'phot':
            return 'Faulty'
        return 'TripleCoinc'

    elif stages_hit == 2:
        if has_1 and has_2:
            valid = row['Process_1'] == 'compt' and row['Process_2'] == 'phot'
        elif has_1 and has_3:
            valid = row['Process_1'] == 'compt' and row['Process_3'] == 'phot'
        elif has_2 and has_3:
            valid = row['Process_2'] == 'compt' and row['Process_3'] == 'phot'
        return 'DoubleCoinc' if valid else 'Faulty'

    return 'Faulty'  # single stage, shouldn't appear but just in case

# Obtain dataframe from ideal information from root file
def get_real_dataframe(root_filename, tree_name):
    with uproot4.open(root_filename) as f:
        tree = f[tree_name]
        df = tree.arrays(['EventID', 'TrackID', 'Xpos', 'Ypos', 'Zpos', 'Etot', 'Process'], library="pd")
        df['Process'] = df['Process'].astype(str)
        return df

# Merge ideal information from all 3 stages in a single dataframe using only single hit events
# Classifies coincidences as double, triple or faulty
def get_true_information(root_filename, drop_faulty=True, stages_hit=2):
    df1 = get_real_dataframe(root_filename, "G4S1TrueGammaInfo")
    df2 = get_real_dataframe(root_filename, "G4S2TrueGammaInfo")
    df3 = get_real_dataframe(root_filename, "G4S3TrueGammaInfo")

    df1_f = filter_single_hits(df1)
    df2_f = filter_single_hits(df2)
    df3_f = filter_single_hits(df3)

    merged = df1_f.merge(df2_f, on=['EventID', 'TrackID'], suffixes=('_1', '_2'), how='outer') \
        .merge(df3_f, on=['EventID', 'TrackID'], how='outer')
    merged.rename(columns={c: f'{c}_3' for c in df3_f.columns if c not in ['EventID', 'TrackID']}, inplace=True)
    # Count how many stages hit for each (EventID, TrackID)
    merged['NHits'] = (merged[[c for c in merged.columns if 'Xpos' in c]].notna().sum(axis=1))
    # Filter by number of stages hit
    merged = merged[merged['NHits'] >= stages_hit].reset_index(drop=True)
    merged['Classification'] = merged.apply(classify_event, axis=1)

    if drop_faulty:
        merged = merged[merged['Classification'] != 'Faulty']

    return merged.drop(['Classification', 'Process_1', 'Process_2', 'Process_3'], axis=1)

# For realistic detector data
# Obtain dataframe from measured information from root file
def get_measured_dataframe(root_filename, tree_name):
    with uproot4.open(root_filename) as f:
        tree = f[tree_name]
        df = tree.arrays(['EventID', 'Line', 'Col', 'Edep'], library="pd")
        df['Xpix'] = df['Col']
        df['Ypix'] = df['Line']
        df = df.drop(['Line', 'Col'], axis=1)
        return df.sort_values(by=['EventID'])

# Merge detector information from all 3 stages in a single dataframe
# performing cluster filtering and converting pixel to coordinates
def get_detector_information(rootfile, Npix=120, pitch=0.5, stages_hit=2, pixlimit=10):
    df1 = get_measured_dataframe(rootfile, "G4Sensor1Hits")
    df2 = get_measured_dataframe(rootfile, "G4Sensor2Hits")
    df3 = get_measured_dataframe(rootfile, "G4Sensor3Hits")
    df_hex1_clustered = cluster_df(df1)
    df_hex2_clustered = cluster_df(df2)
    df_hex3_clustered = cluster_df(df3)
    df1_filt = filter_by_clusters(df_hex1_clustered, pixlimit=pixlimit)
    df2_filt = filter_by_clusters(df_hex2_clustered, pixlimit=pixlimit)
    df3_filt = filter_by_clusters(df_hex3_clustered, pixlimit=pixlimit)
    merged = df1_filt.merge(df2_filt, on='EventID', suffixes=('_1', '_2'), how='outer') \
        .merge(df3_filt, on='EventID', how='outer')
    merged.rename(columns={c: f'{c}_3' for c in df3_filt.columns if c not in ['EventID']}, inplace=True)
    merged['stages_hit'] = (merged[[c for c in merged.columns if 'Xcm' in c]].notna().sum(axis=1))

    # Filter by number of stages hit
    df_all = merged[merged['stages_hit'] >= stages_hit].reset_index()
    df_all['Xcm_1'], df_all['Ycm_1'] = pixel_to_coordinate(df_all['Xcm_1'], df_all['Ycm_1'], Npix, pitch)
    df_all['Xcm_2'], df_all['Ycm_2'] = pixel_to_coordinate(df_all['Xcm_2'], df_all['Ycm_2'], Npix, pitch)
    df_all['Xcm_3'], df_all['Ycm_3'] = pixel_to_coordinate(df_all['Xcm_3'], df_all['Ycm_3'], Npix, pitch)
    return df_all.drop(['stages_hit'], axis=1)

# ------------------------------------------------------
# ------------ RECONSTRUCTION FUNCTIONS ----------------
# ------------------------------------------------------

# Get a position vector for each voxel in the imaging volume
def get_position_matrix(Xbins, Ybins, Zbins, xmin, xmax, ymin, ymax, zmin, zmax):
    x_positions = np.linspace(xmin, xmax, Xbins)
    y_positions = np.linspace(ymax, ymin, Ybins)
    z_positions = np.linspace(zmin, zmax, Zbins)
    X, Y, Z = np.meshgrid(x_positions, y_positions, z_positions, indexing='ij')
    position_matrix = np.stack([X, Y, Z], axis=-1)
    return position_matrix

# Build a gaussian prior in 3D space. X axis is constant, Y and Z axes are gaussian.
def gaussian_prior(Xbins, Ybins, Zbins, sigma_y, sigma_z, mean_y=0, mean_z=0):
    y = np.linspace(-3, 3, Ybins)
    z = np.linspace(-3, 3, Zbins)

    gaussian_y = (1 / np.sqrt(2 * np.pi * sigma_y**2)) * np.exp(-0.5 * ((y-mean_y)/sigma_y)**2)
    gaussian_z = (1 / np.sqrt(2 * np.pi * sigma_z**2)) * np.exp(-0.5 * ((z-mean_z)/sigma_z)**2)

    gaussian_2d = np.outer(gaussian_y, gaussian_z)

    array = np.tile(gaussian_2d[np.newaxis, :, :], (Xbins, 1, 1))
    return array

# get compton angles from detector information
def get_compton_angle(delta_E1, delta_E2, E0):
    m_electron=0.511
    E0 = np.where(E0 is None, delta_E1 + delta_E2, E0) # In case E0 is not provided assume full energy deposition
    Emax = E0 / (1 + (m_electron / (2 * E0)))
    cos_angle = 1 - m_electron * ((1 / (E0 - delta_E1)) - (1 / E0))
    # Mask out invalid rows (where delta_E1 >= Emax, or cos_angle out of [-1, 1])
    invalid = (delta_E1 >= Emax) | (cos_angle < -1) | (cos_angle > 1)
    angle = np.arccos(np.clip(cos_angle, -1, 1))
    return np.where(invalid, np.nan, angle)

# Function that finds the incoming photon energy in a 3 stage CC
# Find the initial kinetic energy of the particle from the first two energy deposits and second scattering angle
# Formula taken from S W Peterson et al 2010 Phys. Med. Biol. 55 6841
def initial_energy(delta_E1, delta_E2, angle_2):
    return delta_E1 + 0.5 * (delta_E2 + np.sqrt(delta_E2 ** 2 + (4 * (delta_E2 * 0.511) / (1 - np.cos(angle_2)))))

# Get the initial Compton angle for double or triple coincident data
# To get E0 from triple coincidences, I use the process described in S W Peterson et al 2010 Phys. Med. Biol. 55 6841
def get_compton_scatters(df, on='cm'):
    df_cc = df.copy()
    # Set masks to know the type of coincidences
    m1 = df_cc['Etot_1'].notna()
    m2 = df_cc['Etot_2'].notna()
    m3 = df_cc['Etot_3'].notna()
    # Classify as double or triple coincidence
    triple = m1 & m2 & m3
    double = (m1 & m2) | (m1 & m3) | (m2 & m3)
    # Initialise outputs
    df_cc['initEnergy'] = np.nan
    df_cc['comptAngle'] = np.nan
    # --- Case 1: Triple coincidences ---
    if triple.any():
        dE1 = df_cc.loc[triple, 'Etot_1']
        dE2 = df_cc.loc[triple, 'Etot_2']

        n1 = get_unit_vector(
            df_cc.loc[triple, 'X'+on+'_1'],
            df_cc.loc[triple, 'Y'+on+'_1'],
            df_cc.loc[triple, 'Z'+on+'_1'],
            df_cc.loc[triple, 'X'+on+'_2'],
            df_cc.loc[triple, 'Y'+on+'_2'],
            df_cc.loc[triple, 'Z'+on+'_2'],
        )
        n2 = get_unit_vector(
            df_cc.loc[triple, 'X'+on+'_2'],
            df_cc.loc[triple, 'Y'+on+'_2'],
            df_cc.loc[triple, 'Z'+on+'_2'],
            df_cc.loc[triple, 'X'+on+'_3'],
            df_cc.loc[triple, 'Y'+on+'_3'],
            df_cc.loc[triple, 'Z'+on+'_3'],
        )
        angles_2 = np.arccos(np.sum(n1 * n2, axis=0))
        E0 = initial_energy(dE1, dE2, angles_2)
        df_cc.loc[triple, 'initEnergy'] = E0
        df_cc.loc[triple, 'comptAngle'] = get_compton_angle(dE1, dE2, E0)
    # --- Case 2: Double coincidences ---
    if double.any():
        # Map energies in interaction order
        dE1 = df_cc['Etot_1'].where(m1, df_cc['Etot_2'])
        dE2 = df_cc['Etot_2'].where(m1 & m2, df_cc['Etot_3'])
        # Restrict to doubles only (exclude triples)
        dE1 = dE1[double & ~triple]
        dE2 = dE2[double & ~triple]
        # Assume full energy deposition
        E0 = dE1 + dE2
        df_cc.loc[double & ~triple, 'initEnergy'] = E0
        df_cc.loc[double & ~triple, 'comptAngle'] = get_compton_angle(dE1, dE2, E0)
    return df_cc.dropna(subset=['comptAngle'])


# Obtain a unit vector that points in the direction of two points (x1, y1, z1) and (x2, y2, z2)
def get_unit_vector(x1, y1, z1, x2, y2, z2):
    delta_x = x2 - x1
    delta_y = y2 - y1
    delta_z = z2 - z1
    magnitude = np.sqrt(delta_x ** 2 + delta_y ** 2 + delta_z ** 2)
    nx = delta_x / magnitude
    ny = delta_y / magnitude
    nz = delta_z / magnitude
    n = np.array([nx, ny, nz])
    return n

# x1, y1, z1 hit position in scatterer
# x2, y2, z2 hit position in absorber
# r is a set of positions in the image volume, already converted to cartesian
# k is a constant that determines the threshold (empirical)
# Returns the positions where the cone and the image volume intersect (allowing for a threshold)
# This is defined as proposed by Mundy and Herman (2011) https://aapm.onlinelibrary.wiley.com/doi/10.1118/1.3519873
def calculate_conical_surface(x1, y1, z1, x2, y2, z2, theta, r, k):
    # vectors and matrices
    r1 = np.array([x1, y1, z1])
    n = get_unit_vector(x1, y1, z1, x2, y2, z2)
    #  calculations
    diff = r - r1  # broadcast → (Y, X, Z, 3)
    dist = np.sqrt(np.sum(diff**2, axis=3))
    threshold = k * dist *np.sin(2*theta)
    # Cone equation terms
    dot = np.sum(diff * n, axis=3)
    norm_sq = np.sum(diff**2, axis=3)
    a = dot**2
    b = (np.cos(theta)**2) * norm_sq
    Strue = np.abs(b - a)  # solution matrix
    cone_positions = np.argwhere(Strue <= threshold)
    return  cone_positions

# Gets the cone positions for every event in a dataframe
def get_cones(df, r, k, on='pos'):
    t0 = time.time()

    good_events = 0
    cones = []

    bar = progressbar.ProgressBar(maxval=len(df), widgets=[progressbar.Bar('=', '[', ']'), ' ',
                                                           progressbar.Percentage()])
    bar.start()

    for i in range(len(df)):
        row = df.iloc[i]
        if np.isnan(row['X'+on+'_1']):
            x1, y1, z1 = row['X'+on+'_2'], row['Y'+on+'_2'], row['Z'+on+'_2']
            x2, y2, z2 = row['X'+on+'_3'], row['Y'+on+'_3'], row['Z'+on+'_3']
        else:
            x1, y1, z1 = row['X'+on+'_1'], row['Y'+on+'_1'], row['Z'+on+'_1']

            if np.isnan(row['X'+on+'_2']):
                x2, y2, z2 = row['X'+on+'_3'], row['Y'+on+'_3'], row['Z'+on+'_3']
            else:
                x2, y2, z2 = row['X'+on+'_2'], row['Y'+on+'_2'], row['Z'+on+'_2']
        theta = row['comptAngle']
        cone = calculate_conical_surface(x1, y1, z1, x2, y2, z2, theta, r, k)
        if cone.size == 0:
            continue
        else:
            cones.append(cone)
            good_events += 1
        bar.update(i + 1)
    bar.finish()

    print(f"Completion time {time.time() - t0:.2f} s")
    print(f"{good_events} out of {len(df)} events were successfully reconstructed.")
    return cones

# Obtain simple, unweighted backprojection image from cone calculation
def simple_back_projection(Xbins, Ybins, Zbins, cones):
    Stot = np.zeros([Xbins, Ybins, Zbins])
    for cone in cones:
        for pos in cone:
            Stot[pos[0], pos[1], pos[2]] += 1
    return Stot

# Performs stochastic origin ensemble reconstruction based on cone information
# Allows for incorporation of weighting using prior (i.e. gaussian or SBP)
# If percent_convergence is provided, the function may stop when the movement probability is below this value
def stochastic_origin_ensemble(cones, Xbins, Ybins, Zbins, N_events, N_soe, weights=None, percent_convergence=None):
    t0 = time.time()
    probabilities = []
    all_SOEs = []

    # Set initial state of the SOE
    D = np.zeros([Xbins, Ybins, Zbins])
    old_positions = []
    for cone in cones:
        # we choose our first position within the current cone
        optimumPos = random.choice(cone)  # returns tuple for 3D position
        D[optimumPos[0], optimumPos[1], optimumPos[2]] += 1
        old_positions.append(optimumPos)
    all_SOEs.append(D.copy())
    # Main loop
    bar = progressbar.ProgressBar(maxval=N_soe, widgets=[progressbar.Bar('=', '[', ']'), ' ',
                                                         progressbar.Percentage()])
    bar.start()
    for j in range(N_soe):
        moves_performed = 0
        for i in range(N_events): # Loop over N events that are subject to movement per iteration
            k = random.randint(0, len(cones)-1)
            old_pos = old_positions[k]
            actual_cone = cones[k]
            new_pos = random.choice(actual_cone)  # get a random position within the cone of event k
            old_density = D[old_pos[0], old_pos[1], old_pos[2]]  # actual density at old position
            new_density = D[new_pos[0], new_pos[1], new_pos[2]]  # actual density at new position

            # Definition of priors
            if weights is not None: # If weights are provided, use them to define the prior
                old_w = weights[old_pos[0], old_pos[1], old_pos[2]]
                new_w = weights[new_pos[0], new_pos[1], new_pos[2]]
            else: # Otherwise, use flat prior
                old_w = 1
                new_w = 1
            if old_density > 0:
                ratio = ((new_density + 1)/ old_density) * (new_w / old_w)
            else:
                ratio = 1
            # acceptance_probability = min(1, ratio)
            if ratio >= random.uniform(0, 1):
                moves_performed += 1
                D[new_pos[0], new_pos[1], new_pos[2]] += 1  # update density value at new position
                D[old_pos[0], old_pos[1], old_pos[2]] -= 1  # update density value at old position
                old_positions[k] = new_pos  # update definition of old position

        # Evaluate movement probability
        # TODO: incorporate notion of entropy as a means for convergence testing
        move_probability = 100*(moves_performed/N_events)
        probabilities.append(move_probability)
        if percent_convergence is not None:
            # break the loop if we reach a stable probability of acceptance, defined as a mean of <(percent_convergence)%
            # and a std < 1 for the last 500 iterations
            if j > 500:
                mean = np.mean(probabilities[j-500:j])
                std = np.std(probabilities[j-500:j])
                if (mean < percent_convergence) and (std < 1):
                    break

        all_SOEs.append(D.copy())
        bar.update(j + 1)
    bar.finish()

    print(f"Completion time {time.time() - t0:.2f} s")
    return all_SOEs, probabilities