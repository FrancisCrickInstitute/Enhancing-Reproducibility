# Overview

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/djpbarry/Enhancing-Reproducibility/main?urlpath=%2Fdoc%2Ftree%2Fnotebooks%2Fcompanion_notebook.ipynb) [![Python 3.11](https://img.shields.io/badge/python-3.11-blue.svg)](https://www.python.org/downloads/release/python-3100/) ![Commit activity](https://img.shields.io/github/commit-activity/y/djpbarry/Enhancing-Reproducibility?style=plastic) ![GitHub](https://img.shields.io/github/license/djpbarry/Enhancing-Reproducibility?color=green&style=plastic)

This repository contains the Python code associated with the following paper:

- Barry DJ, Marcotti S, Gerontogianni L and Kelly G (2025). Enhancing Reproducibility Through Bioimage Analysis: The Significance of Effect Sizes and Controls.

# Get Started

The quickest and easiest way to try this code is to [try it on Binder](https://mybinder.org/v2/gh/djpbarry/Enhancing-Reproducibility/main?urlpath=%2Fdoc%2Ftree%2Fnotebooks%2Fcompanion_notebook.ipynb). This will allow you to reproduce the plots in the associated publication.

# Run On Your Own Data

To test our code on your own data, the easiest thing to do is download this repo and run the [Nuclear_Fascin.cppipe](Nuclear_Fascin.cppipe) [CellProfiler](https://cellprofiler.org/) pipeline on your own images and replace the files in [the cell_profiler_outputs](inputs/cell_profiler_outputs) folder. You can then use the [Jupyter Notebook](notebooks/companion_notebook.ipynb) to generate plots for your own images.

A step-by-step guide is presented below. **You only need to perform steps 1 and 2 once.** Every subsequent time you want to run the code, skip straight to step 3.

## Step 1: Install a Python Distribution

We recommend using conda as it's relatively straightforward and makes the management of different Python environments simple. You can install conda from [here](https://conda.io/projects/conda/en/latest/user-guide/install/index.html#regular-installation) (miniconda will suffice).

## Step 2: Set Up Environment

Once conda is installed, open Anaconda Prompt and run the following series of commands:

```
conda create --name enhancing-reproducibility pip
conda activate enhancing-reproducibility
python -m pip install -r <path to this repo>/requirements.txt
```
where you need to replace `<path to this repo>` with the location on your file system where you downloaded this repo. You will be presented with a list of packages to be downloaded and installed. The following prompt will appear:
```
Proceed ([y]/n)?
```
Hit Enter and all necessary packages will be downloaded and installed - this may take some time. When complete, you can deactivate the environment you have created with the following command.

```
conda deactivate
```
You have successfully set up the necessary conda environment!

## Step 3: Run The Code!

The following commands will launch a Jupyter notebook allowing you to run the code on your own data:
```
conda activate enhancing-reproducibility
jupyter notebook <path to this repo>/notebooks/companion_notebook.ipynb
```

The Jupyter Notebook should open in your browser - follow the step-by-step instructions in the notebook to run the code. If you are not familiar with Jupyter Notebooks, you can find a detailed introduction [here](https://jupyter-notebook.readthedocs.io/en/latest/notebook.html#introduction).
