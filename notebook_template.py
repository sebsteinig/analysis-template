# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.13.8
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# + [markdown] colab_type="text" id="view-in-github" tags=[]
# <a href="https://colab.research.google.com/github/sebsteinig/analysis_template/blob/main/notebook_template.ipynb" target="_parent"><img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/></a>
# -

# # Title
# short description of what the notebook does

# ## Prelude (only necesseary when running on Google Colab)
# If running on Google Colab, execute the following two cells seperately to download the files and install the necessary conda environment on the virtual machine. This will take several minutes and involves restarting the kernel.
#
# If running somewhere else, you can execute the whole notebooks and this part will be skipped automatically.

# +
# configure environment on Colab
try:
    import google.colab

    # if on Colab, clone the repository to access the data locally
    import os
    repo = "analysis_template"

    # clone repo if it does not already exist
    if not os.path.exists(repo):
        print('cloning GitHub repository: ' + repo)
        # !git clone https://github.com/sebsteinig/{repo}.git
        
    # %cd /content/{repo}

    # install condacolab to easily handle conda environments on Colab
    # !pip install -q condacolab
    import condacolab
    condacolab.install()
    
except:
    print('not running on Google Colab')

# + colab={"base_uri": "https://localhost:8080/"} id="DmiZzoIu_5lO" outputId="46885ef8-d49b-49ad-a1d7-844393356d5b"
try:
    import google.colab
    
    # install packages from environment.yml file
    # !conda env update -n base -f environment.yml
    
except:
    print('not running on Google Colab')
# -

# ## User input
# define variables/lists to quickly change inputs to the notebook

# +
work_dir       = '.' # location of cloned repository
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

#### suppress warnings
#import warnings
#warnings.filterwarnings('ignore')

# -

# ## Main code







# #### appendix: some code snippets I regularly use

# + colab={"base_uri": "https://localhost:8080/", "height": 713} id="055047f3" outputId="8b79a52d-0d53-4631-91b9-e4fb17cfc984" tags=[]
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

