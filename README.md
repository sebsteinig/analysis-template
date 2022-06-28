# boilerplate python notebook template
This repo should contain all documentation, scripts and data to recreate the final figure/analysis from the source data. Final analysis and plotting is done in Jupyter notebooks with a python kernel for interactive analysis of the source or preprocessed data and inline documentation. 

## purpose
What do I want to do?

## source data
Some external, static data that builds the starting point of the analysis (e.g model output, CMIP archive, OPeNDAP server, ...) 

## preprocessing
Large source data should be preprocessed before importing into the jupyter notebook (e.g., regridding, averaging, ...) to speed up interactive analysis and to allow distribution via github (100 MB file size limit). Put the preprocessing script into this repo and describe what you have done with the source data

## running the notebooks
Notebooks can either be run on [Google Colab](https://colab.research.google.com/) (online, Google account required) or locally. Notebooks should include a button to open directly in Colab, otherwise you can also directly load a GitHub repo. Easiest way to run locally is to install [conda](https://conda.io/projects/conda/en/latest/index.html) 
and create an environment `env_name` with 

```
conda env create --name env_name --file=environment.yml
``` 

using the `environment.yml` file from this repository to install all necessary packages. The notebooks can then be run 
interactively by typing

```
jupyter lab
```

or directly from the command line with

```
jupyter nbconvert --to notebook --inplace --execute notebook_name.ipynb
```

