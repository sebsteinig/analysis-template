# boilerplate python notebook template
This repo should contain all documentation, scripts and data to recreate the final figure/analysis from the source data. Final analysis and plotting is done in jupyter notebooks with a python kernel for interactive analysis of the source or preprocessed data and inline documentation. 

## purpose
What do I want to do?

## source data
Some external, static data that builds the starting point of the analysis (e.g model output, CMIP archive, OPeNDAP server, ...) 

## preprocessing
Large source data should be preprocessed before importing into the jupyter notebook (e.g., regridding, averaging, ...) to speed up interactive analysis and to allow distribution via github (100 MB file size limit). Put the preprocessing script into this repo and describe what you have done with the source data

## running the notebooks
Notebooks can either be run on [Google Colab](https://colab.research.google.com/) (online, Google account required) or locally. Notebooks should include a button to open it in Colab, otherwise you can also directly load a GitHub repo within Colab. Easiest way to run locally is to first download the repo with


```
git clone https://github.com/USERNAME/REPOSITORY
``` 

and then install [conda](https://conda.io/projects/conda/en/latest/index.html) (if not installed already). Then create an environment `env_name` with 

```
conda env create --name env_name --file=environment.yml
``` 

using the `environment.yml` file from this repository to install all necessary python packages. The notebooks can then be run interactively by typing

```
jupyter lab
```

We further use the [jupytext](https://jupytext.readthedocs.io/en/latest/index.html) package to automatically save a pure python script (*.py) each time the notebook is saved. This is very helpful for clean diffs in the version control and allows you to run the analysis in your local terminal with:

```
python notebook_name.py
```
The python file can also be shared with others to work on the code together using all the version control benefits (branches, pull requests, ...). You can edit it with any tex editor/IDE and it can also be converted back to a jupyter notebook (with no output) via
```
jupytext --to notebook notebook_name.py
```

