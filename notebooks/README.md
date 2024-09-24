# How to setup the companion notebook

## Download the material from GitHub
1. Navigate to the [GitHub repo](https://github.com/djpbarry/dont-trust-p-values)
2. Locate the blue `<> Code` button on the top right of the page
3. Copy the material locally by either cloning the repository or hitting on `Download ZIP`

## Setup a Python environment
This step implies you already have a Python environment manager (such as Miniconda) installed on your machine.
If this is not the case, please find more instructions [here](https://github.com/RMS-DAIM/introduction-to-image-analysis/blob/main/Pages/Installation-Instructions.md#installing-conda).
1. Open a terminal 
2. Create a new Python environment with the following command `conda create -n p-val-env`
3. Access the newly created environment by typing `conda activate p-val-env`
4. Install the required packages by typing each of the following commands one at a time:
    * `conda install pip`
    * `pip install jupyterlab`
    * `pip install matplotlib`
    * `pip install pandas`
    * `pip install scikit_posthocs`
5. Navigate to the folder where you saved the material in the previous step with `cd /path/to/folder/`
6. Launch Jupyter lab with the following command `jupyter-lab`

## Open the notebook
1. In Jupyter lab, use the lefthand tab to navigate to the folder `notebooks`
2. Open the file name `companion_notebook.ipynb`
3. Run cells sequentially to reproduce the figures in the paper
