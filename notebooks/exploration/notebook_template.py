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

# # Title
# short description of what the notebook does

# ## User input
# define variables/lists to quickly change inputs to the notebook

# +
work_dir       = '.' # location of cloned repository
data_dir       = work_dir + '/data' # location of data files
fig_dir        = work_dir + '/figures' # location of figure files

exp_list       = ['exp1', 'exp2'] # list of data sets to loop over

save_figures   = True # flag whether to save figures to disk or not
# -

# ## Load packages

# + colab={"base_uri": "https://localhost:8080/"} id="WcMt0hmK_38Q" outputId="1ca62e3a-a818-48f8-aea4-d2827f558cdb"
### some standard packages
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

### some optional packages that I always have to google ...

#### cartopy maps
#from cartopy import config
#import cartopy.crs as ccrs
#from cartopy.util import add_cyclic_point

#### colormaps
#import cmocean

#### csv parser
#import csv


# -

# ## Main code







# #### appendix: some code snippets I regularly use

# + colab={"base_uri": "https://localhost:8080/", "height": 713} id="055047f3" outputId="8b79a52d-0d53-4631-91b9-e4fb17cfc984"
#### loop analysis over data sets   
# for expCount, exp in enumerate(exp_list):

#### load netcdf data set
# ds = xr.open_dataset(work_dir + '/data/file_name.nc')

#### new multi-panel figure
# fig, axes = plt.subplots(nrows, ncols, constrained_layout=True, figsize=(width, height) ) # figsize in inches

#### map plot with cartopy
# ax = fig.add_subplot(nrows, ncols, index, projection=ccrs.Robinson()) # or e.g. ccrs.PlateCarree()
# ax.set_extent([minlon,maxlon, minlat,maxlat], ccrs.PlateCarree()) # or ax.set_global()
# ax.coastlines()
# ax.contourf(ds['variable_name'], transform=ccrs.PlateCarree(), levels=21, 
#             vmin=..., vmax=..., cmap=cmocean.cm.topo, add_colorbar=False)

#### add cyclic longitude to field and coordinate (from cartopy.util import add_cyclic_point)
# variable_cyclic, longitude_cyclic = add_cyclic_point(variable, coord=longitude)

#### save figure
# if save_figures:
#      plt.savefig(work_dir + '/figures/figure_name.pdf')  
#      plt.savefig(work_dir + '/figures/figure_name.png', dpi=200)  

