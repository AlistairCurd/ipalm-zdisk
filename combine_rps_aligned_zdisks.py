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

from scipy.spatial import KDTree

# from matplotlib import colors
# -


# # Find Z-disk localisations

# ## Set directory here
# Replace the string.

# +
# locdatapath = 'C://Janelia_Z-disk_2024//240117-3_Run1-640_c123_sum_X11_processed_IDL_purged_IDL_beadsremoved_IDL_ASCII_xy-in-nm.txt'
locdatapath = 'T://Data//Janelia-Z-disk-2024//segmentations for analysis//z_disk//output//segmented_z_disks'

locdatapath = Path(locdatapath)


# -

# ## Function to get localisations from aligned Z-disks

def get_aligned_locs(input_table_path, filter_params):
    """Extract localisations from aligned Z-disk data.

    Args:
        input_table_path (pathlib Path):
            Path to a csv table of localisationsm including in aligned z-disks after PCA.
            csv table will have desired coordinates in 'x', 'y' and 'z'.
        filter_params (dict):
            Parameters (keys) and values for filtering the localisations.
    Returns:
        locs_xyz (numpy Array):
            Array of localisations, with 3 columns (xyz coordinates)
            x and y are swapped compared with the input table, because PCA gives the short XY axis as X,
            but I want it to be Y.
    """
    locs_df = pd.read_table(input_table_path, delimiter=',')
    locs_df = locs_df[(locs_df['Sigma X Pos rtNph'] < filter_params['loc_prec_xy_nm_max']) \
        & (locs_df['Sigma Y Pos rtNph'] < filter_params['loc_prec_xy_nm_max'])]
    locs_xyz = locs_df[['y', 'x', 'z']].to_numpy()
    return locs_xyz


for file in locdatapath.iterdir():
    if 'aligned' in file.name:
        print(file.name)
        locs_xyz = get_aligned_locs(file, filter_params)
        print(locs_xyz.shape)
        kd_tree = KDTree(locs_xyz)
        pairs = kd_tree.query_pairs(r=150, output_type='ndarray')
        print(pairs.shape)

pairs[0:5]











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


