#!/bin/bash


# Install Miniconda if it's not already installed
if [ ! -d "$HOME/miniconda3" ]; then
    echo "== Miniconda not found, installing =="
    # Download the latest Miniconda installer for macOS
    curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
    # Install Miniconda in batch mode (no prompts)
    bash Miniconda3-latest-MacOSX-x86_64.sh -b -p $HOME/miniconda3
    # Initialize conda for bash
    $HOME/miniconda3/bin/conda init bash
    # Clean up the installer
    rm Miniconda3-latest-MacOSX-x86_64.sh
    echo "== Miniconda installed =="
    # Source bashrc to make conda available in current shell
    # Source appropriate profile based on OS
    if [[ "$OSTYPE" == "darwin"* ]]; then
        # macOS uses bash_profile
        source ~/.bash_profile
    else
        # Linux typically uses bashrc
        source ~/.bashrc
    fi
    
else
    echo "== Miniconda already installed =="
fi

export CONDA_ALWAYS_YES="true"

if [[ "$OSTYPE" == "darwin"* ]]; then
    # macOS uses bash_profile
    source ~/.bash_profile
else
    # Linux typically uses bashrc
    source ~/.bashrc
fi

# check if already exists
conda list -n dask

if [ $? != 0 ] ; then
    echo "== conda env not found, installing =="
    conda create -n dask python=3.10
fi

echo "== activating conda env dask =="
conda activate dask

conda info --env

# Install Dask and its dependencies
conda install dask==2024.6.1

# Install Dask-Jobqueue
conda install dask-jobqueue==0.8.5 -c conda-forge

# Install JupyterLab and NodeJS
conda install nodejs==6.13.1 jupyterlab==4.2.3 -c conda-forge -y

# Install Dask-JupyterLab extension
pip install dask_labextension==7.0.0
# Integrate Dask-Labextension with Jupyter (requires NodeJS)
jupyter labextension install dask-labextension

unset CONDA_ALWAYS_YES