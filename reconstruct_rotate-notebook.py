# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.15.2
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# # IMPORTANT
# # Disable autosave for Jupytext version control with a paired .py script
# ## But manually saving the notebook frequently is still good

# %autosave 0
# # Acquisition system constants

# nm per camera pixel (XY)

nm_per_px = 133.3

# # Filtering parameters

# max localisation precision in x and y to include

filter_params = {'loc_prec_xy_nm_max': 7}

# ***Filter for Z error as well?***

# # Imports

# +
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# from matplotlib import colors
# -


# # Set path for the localisation here
# Replace the string.

# +
locdatapath = 'C://Janelia_Z-disk_2024//240117-3_Run1-640_c123_sum_X11_processed_IDL_purged_IDL_beadsremoved_IDL_ASCII_xy-in-nm.txt'

locdatapath = Path(locdatapath)
# -

# # Load localisations into a dataframe
# ## Use localisations exported from PeakSelector with X and Y in nm
# The precision estimates (`Sigma X [or Y] Pos rt Nph`) are then also in nm

locs_df = pd.read_table(locdatapath)
locs_df

# ## Print column names
#
# `Sigma X [or Y] Pos rt Nph` is a localisation precision estimate in units of pixels

locs_df.columns

# # Filter XYZ data

# ## Remind ourselves of the filtering parameters

filter_params

# ## Filter and select only XYZ columns

locs_filtered_df = locs_df[(locs_df['Sigma X Pos rtNph'] < filter_params['loc_prec_xy_nm_max']) & (locs_df['Sigma Y Pos rtNph'] < filter_params['loc_prec_xy_nm_max'])]
locs_filtered_df = locs_filtered_df[['X Position', 'Y Position', 'Unwrapped Z', 'Sigma X Pos rtNph', 'Sigma Y Pos rtNph', 'Unwrapped Z Error']]
locs_filtered_df

locs_filtered_df['Sigma Y Pos rtNph'].max()

locs_df['Sigma X Pos rtNph'].max()

# ## Put XYZ coordinates in arrays for easy reference

xyz_coords_nm = locs_filtered_df[['X Position', 'Y Position', 'Unwrapped Z']].to_numpy()
xyz_coords_nm.shape
x_nm = xyz_coords_nm[:, 0]
y_nm = xyz_coords_nm[:, 1]
z_nm = xyz_coords_nm[:, 2]

# # Plot FOV

# ## Set plotting parameters
# * bin-size (nm)
# * max displayed histogram bin value (image saturation point)

binsize = 100
hist_im_sat_value = 50

# ## min and max XYZ for reference while plotting

print('\tx\ty\tz')
print('min\t{}\t{}\t{}'.format(xyz_coords_nm[:, 0].min(), xyz_coords_nm[:, 1].min(), xyz_coords_nm[:, 2].min()))
print('max\t{}\t{}\t{}'.format(xyz_coords_nm[:, 0].max(), xyz_coords_nm[:, 1].max(), xyz_coords_nm[:, 2].max()))

# ## Plot

fig, ax = plt.subplots()
ax.set_aspect(1)
counts, xedges, yedges, im = \
    ax.hist2d(x_nm, y_nm,
              bins=[int(np.ceil(np.ptp(x_nm) / binsize)), int(np.ceil(np.ptp(y_nm) / binsize))],
              cmap='inferno',
              vmin=0, vmax=hist_im_sat_value
              )
fig.colorbar(im, ax=ax)

# # Crop FOV

# ## Set min and max x and y to use

# +
x_min = 10000
x_max = 30000

y_min = 13000
y_max = 41000
# -

# ## Get cropped data

locs_cropped_df = locs_filtered_df[
    (x_nm > x_min) & (x_nm < x_max)
    & (y_nm > y_min) & (y_nm < y_max)
    ]
locs_cropped_df

# ### Check

print(locs_cropped_df['X Position'].min())
print(locs_cropped_df['X Position'].max())
print(locs_cropped_df['Y Position'].min())
print(locs_cropped_df['Y Position'].max())

x_cropped = locs_cropped_df['X Position']
y_cropped = locs_cropped_df['Y Position']
z_cropped = locs_cropped_df['Unwrapped Z']

# ## Plot

fig, ax = plt.subplots()
ax.set_aspect(1)
counts, xedges, yedges, im = \
    ax.hist2d(x_cropped, y_cropped,
              bins=[int(np.ceil(np.ptp(x_cropped) / binsize)), int(np.ceil(np.ptp(y_cropped) / binsize))],
              cmap='inferno',
              vmin=0, vmax=hist_im_sat_value
              )
fig.colorbar(im, ax=ax)

# # Rotate cropped coordinates to make X the axis of the transformed myofibril

# ## Set rotation angle, from original, in degrees

rotation_angle_degrees = 63

# ## Rotate localisation data
# The rotation happens about the (x = 0, y = 0) axis

# +
# Polar coordinates of locs
r = np.sqrt(x_cropped ** 2 + y_cropped ** 2)
theta = np.arctan2(y_cropped, x_cropped)

# Rotate the coordinates
rotation_angle_rad = rotation_angle_degrees * np.pi / 180

newtheta = theta + rotation_angle_rad
x_rotated = r * np.cos(newtheta)
y_rotated = r * np.sin(newtheta) # * -1 # * -1 to avoid flipping with
                                 # this angle convention
z_rotated = z_cropped
# -

# ## Plot

fig, ax = plt.subplots()
ax.set_aspect(1)
counts, xedges, yedges, im = \
    ax.hist2d(x_rotated, y_rotated,
              bins=[int(np.ceil(np.ptp(x_rotated) / binsize)), int(np.ceil(np.ptp(y_rotated) / binsize))],
              cmap='inferno',
              vmin=0, vmax=hist_im_sat_value
              )
fig.colorbar(im, ax=ax)

# # Crop again

# ## Set min and max x and y to use

# +
x_min = -30000
x_max = 0

y_min = 30000
y_max = 34000
# -

# ## Get cropped data

x_rotated_cropped = x_rotated[
    (x_rotated > x_min) & (x_rotated < x_max)
    & (y_rotated > y_min) & (y_rotated < y_max)
    ]
y_rotated_cropped = y_rotated[
    (x_rotated > x_min) & (x_rotated < x_max)
    & (y_rotated > y_min) & (y_rotated < y_max)
    ]
z_rotated_cropped = z_rotated[
    (x_rotated > x_min) & (x_rotated < x_max)
    & (y_rotated > y_min) & (y_rotated < y_max)
    ]

# ## Plot

fig, ax = plt.subplots()
ax.set_aspect(1)
counts, xedges, yedges, im = \
    ax.hist2d(x_rotated_cropped, y_rotated_cropped,
              bins=[int(np.ceil(np.ptp(x_rotated_cropped) / binsize)), int(np.ceil(np.ptp(y_rotated_cropped) / binsize))],
              cmap='inferno',
              vmin=0, vmax=hist_im_sat_value
              )
fig.colorbar(im, ax=ax)

# # Save

xyz_rotated_cropped = np.vstack((x_rotated_cropped, y_rotated_cropped, z_rotated_cropped)).T

x_rotated_cropped.shape

xyz_rotated_cropped.shape

outfilename = locdatapath.stem + '_MFalongX.csv'
outpath = locdatapath.parent / outfilename

np.savetxt(outpath, xyz_rotated_cropped, delimiter=',')


