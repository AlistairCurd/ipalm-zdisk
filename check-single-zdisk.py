# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.16.1
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

# # nm per camera pixel (XY)

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


# # Check one Z-disk

# ## Set path for the localisations here
# Replace the string.

# +
# locdatapath = 'C://Janelia_Z-disk_2024//240117-3_Run1-640_c123_sum_X11_processed_IDL_purged_IDL_beadsremoved_IDL_ASCII_xy-in-nm.txt'
locdatapath = 'T://Data//Janelia-Z-disk-2024//segmentations for analysis//z_disk//output//segmented_z_disks//18_Run1-640_c123_sum_X10_processed_purged_beads_removed_Grouped_peaks_nm_IDL_ASCII_zdisk_1_aligned.csv'

locdatapath = Path(locdatapath)
# -

# ## Load localisations into a dataframe
# ### Use localisations exported from PeakSelector with X and Y in nm
# The precision estimates (`Sigma X [or Y] Pos rt Nph`) are then also in nm

locs_df = pd.read_table(locdatapath, delimiter=',')
locs_df

# ### Print column names
#
# `Sigma X [or Y] Pos rt Nph` is a localisation precision estimate in units of pixels

locs_df.columns

# ## Filter XYZ data

# ### Remind ourselves of the filtering parameters

filter_params

# ### Filter and select only XYZ columns
# x, y, z are after alignmnt by PCA

locs_filtered_df = locs_df[(locs_df['Sigma X Pos rtNph'] < filter_params['loc_prec_xy_nm_max']) & (locs_df['Sigma Y Pos rtNph'] < filter_params['loc_prec_xy_nm_max'])]
locs_filtered_df = locs_filtered_df[['X Position', 'Y Position', 'Unwrapped Z', 'Sigma X Pos rtNph', 'Sigma Y Pos rtNph', 'Unwrapped Z Error', 'x', 'y', 'z']]
locs_filtered_df

locs_filtered_df['Sigma Y Pos rtNph'].max()

locs_df['Sigma X Pos rtNph'].max()

# ### Put XYZ coordinates in arrays for easy reference

# xyz_coords_nm = locs_filtered_df[['X Position', 'Y Position', 'Unwrapped Z']].to_numpy()
xyz_coords_nm = locs_filtered_df[['x', 'y', 'z']].to_numpy()
xyz_coords_nm.shape
x_nm = xyz_coords_nm[:, 0]
y_nm = xyz_coords_nm[:, 1]
z_nm = xyz_coords_nm[:, 2]

# ## Plot FOV

# ### Set plotting parameters
# * bin-size (nm)
# * max displayed histogram bin value (image saturation point)

binsize = 5
hist_im_sat_value = 2

# ### min and max XYZ for reference while plotting

print('\tx\ty\tz')
print('min\t{}\t{}\t{}'.format(xyz_coords_nm[:, 0].min(), xyz_coords_nm[:, 1].min(), xyz_coords_nm[:, 2].min()))
print('max\t{}\t{}\t{}'.format(xyz_coords_nm[:, 0].max(), xyz_coords_nm[:, 1].max(), xyz_coords_nm[:, 2].max()))

# ### Plot

fig, ax = plt.subplots()
ax.set_aspect(1)
counts, xedges, yedges, im = \
    ax.hist2d(x_nm, y_nm,
              bins=[int(np.ceil(np.ptp(x_nm) / binsize)), int(np.ceil(np.ptp(y_nm) / binsize))],
              cmap='inferno',
              vmin=0, vmax=hist_im_sat_value
              )
fig.colorbar(im, ax=ax)


