# ipalm-zdisk
Process Z-disk localisations from the Janelia iPALM

Developed on Windows 10 in Python with Mamba, a cross-platform manager designed as a drop-in replacement for conda.

Visualise, crop and rotate localisation distribution from iPALM data.

## Installation

### Environment

#### Before first use: create environment

1. Clone or download this repository.

2. Install Miniforge with Mamba using Mambaforge (optionally stick to Miniconda or Anaconda if you have it and want to).

3. In Miniforge Prompt (installed in the previous step):

  * Navigate to the cloned repository
  * Input `mamba env create -f environment.yml` (or `conda ...` with Miniconda or Anaconda).

This will create an environment called ipalm-zdisk with the required packages.

4. Activate the environment

  * Input `mamba activate ipalm-z-disk` (or `conda ...`)

*Test?*

#### For subsequent use: activate the environment

In Miniforge Prompt:

* Input `mamba activate ipalm-zdisk` (or `conda ...`)

## Usage

* Navigate to the local copy of this repository.
* Input `jupyter notebook`.
* Open `reconstruct_rotate-notebook.py`. Open it as a notebook (sometimes *Open With...*), not in a code editor.
* Run the first cell to prevent autosaving, which can cause conflicts with the paired notebooks that are created. Save manually if desired.

## Developing notebooks

Notebooks are paired with .py scripts using Jupytext, for convenient version control.

**Open the required .py script ending in "-notebook" with Jupyter.**
Saving it will generate .ipynb notebooks with the same name. .gitignore is set to ignore the .ipynb files.

If developing and switching between Git branches, work by using the .py versions of the notebooks in Jupyter that are version controlled. Jupytext enables this working well. Open as a notebook (from *Open With* in Jupyter in some versions), not in a code editor. 

Run the first cell that stops Jupyter autosaving, and save progress manually to do avoid accidental conflicts. If you checkout a different branch and open an .ipynb notebook, it will not have changed with the checkout - you need to open and save the .py version first, which does change with version control.

### When creating a new notebook:
Include the first cell that deactivates autosave.
In Jupyter File > Jupytext, select Pair Notebook with light Script. This will cause a paired .py script to be generated when you save the notebook.
