# **Stack Exchange Analysis**

We begin by providing instructions for running the code locally [[Run Locally](#run-locally)], on a Windows machine (though it can be easily adapted to run on any platforms). Next, we explain how to set up the environment and execute the code in an HPC environment which uses a SLURM scheduler [[HPC with SLURM](#hpc-with-slurm)].


# Table of Contents

- [Run Locally](#run-locally)
- [HPC with SLURM](#hpc-with-slurm)

<br>


# Run Locally

## Download the Stack Exchange Dataset

The first step is to get access to the Stack Exchange dataset, which is publicly available at (https://archive.org/details/stackexchange). In the `7Z FILES` tab, select for example `math.stackexchange.com.7z` ('math') and save it into the `data` folder. 

*Note:* StackOverflow ('cs', i.e. Computer Science in the paper) comprises several `.7z` files, download them in the `data` folder. `7z_extraction.py` will perform a custom extraction and standardize the data format within the `data` folder.

## Install

To run the analysis locally, you need to have Python available. Create two Python virtual environments (using Python version 3.10.11) with the dependencies listed in `requirements_load.txt` and `requirements_dask.txt`.

- The first environment (e.g., `.load`) is used for loading and preprocessing the data, as well as producing the plots presented in the paper.
- The second environment (e.g., `.dask`) is used to set up a local Dask cluster and compute the time the *best* answer needs to become stably the best (*time to emerge*).


## Data Loading and Preprocessing

Activate the `.load` virtual environment and, on Windows, run the Powershell Script:

```powershell
./_win_load.ps1
```

It will extract, load and preprocess the `math` dataset. To load the `cs` (Computer Science) one, just modify the environment variables (see the comments within). Make sure to have all the necessary `.7z` files in the `data` folder.


## Dask Local Processing

Activate the `.dask` virtual environment and instantiate a Jupyter server by running:

```powershell
jupyter lab
```

JupyterLab serves as an interactive development environment for notebook-based workflows. Specifically, you can now engage with the `process_dataset_dask.ipynb`, which by default has its `mode` set to "local_debug".

For monitoring your local cluster, a convenient dashboard is accessible via your browser (see [Browser](#browser)).

## Plotting

The information needed for the plots is stored in `distrib` and `results` folder. To obtain figure 8 (a) and (b), run:

```powershell
python plot_answers_ccdfs.py
python plot_answers_ccdfs.py --dataset math
```

While to obtain Figure 9 run:

```powershell
python plot_best_emerge.py
```

The figures (included in the paper) will be stored in `figs\paper` directory. Please be aware that using an updated Stack Exchange dataset with more recent information may lead to minor differences in the generated plots.

<br>

# HPC with SLURM

The data preparation and preprocessing steps are identical to the local procedure; see [Data Loading and Preprocessing](#data-loading-and-preprocessing) for details. Similarly, the plotting process remains the same (refer to [Plotting](#plotting)). The Dask configuration, however, has some minor differences, which are outlined below.

## Dask Distributed

Below guidelines refer to HPC environment with SLURM scheduler. Refer to the [official documentation](#additional-resources) for usage in different  environment. 

*Note:* change the `mode` from "local_debug" (default) to "hpc" in the `process_dataset_dask.ipynb` notebook.

## Install

We provide an `_dask_install.sh` script to automate installation. It will create a `conda` environment and install:

* **Dask Job Queue**
* **JupyterLab Dask extension** 

## Run instructions

### Step 1. Dask configuration

Create a configuration file for Dask by making a directory and adding your configuration:

```bash
mkdir -p ~/.config/dask
```

Create `~/.config/dask/jobqueue.yaml` with the following content (adjust settings as needed for your cluster):

* `network interface name` is usually `ib0` in HPC with Infiniband interconnect. 
* `queue` can be found by looking at the value of `PARTITION` column in one of your jobs, listed by `squeue`



```yaml
jobqueue:
  slurm:
    name: dask-worker
    # replace below with your custom configurations
    interface: <network interface name>            
    queue: <slurm queue>
    account: <youruser>
    
    scheduler-options: {}

    # modify dashboard link in jupyter notebook
    distributed:
      dashboard:
        link: "/proxy/{port}/status"
```


### Step 2: Sanity check

This code will start a Dask cluster over Slurm and suspend for 4 minutes. Run as a Python script or Jupyter notebook:

```python
from dask_jobqueue import SLURMCluster
from dask.distributed import Client
import time

# Create a SLURM cluster
cluster = SLURMCluster(
    cores=4,                # Number of cores per job
    memory="8GB",           # Memory per job
    processes=1,            # Number of processes per job
)

# Scale the cluster (start 10 worker jobs)
cluster.scale(4)

# Connect to the cluster
client = Client(cluster)

# At this point, you can use regular Dask patterns
# For example:
import dask.array as da
x = da.random.random((10000, 10000), chunks=(1000, 1000))
result = x.mean().compute()

time.sleep(240)

# When done, close the client and cluster
client.close()
cluster.close()
```

You can check the output of `squeue` command after the code suspends and verify 4 `dask-worker` jobs have started.

```
$ squeue -u <youruser>
        JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
    1196319    global dask-wor xxx  R       1:02      1 compute-2-23
    1196320    global dask-wor xxx  R       1:02      1 compute-2-23
    1196321    global dask-wor xxx  R       1:02      1 compute-2-23
    1196322    global dask-wor xxx  R       1:02      1 compute-2-23
```


## Troubleshooting

1. **Network Interface**: Make sure you specify the correct network interface (e.g., `ib0` for InfiniBand) in your configuration.
2. **Job Submission Failures**: Check SLURM logs if workers aren't starting. Common issues include:
    - Incorrect account/partition names
    - Resource constraints (memory, cores)
    - Walltime limits

## Monitoring Your Cluster with Dashboard

You can monitor execution via the interactive dashboard. 

**Step 1.** Set a password for the dashboard for convenience `jupyter lab password` (this avoid the need to copy and paste a randomly generated token at every execution). 

**Step 2.** You can access the dashboard in two ways:

#### Jupyter notebook
When running through JupyterLab with the Dask extension installed, you'll see a Dask tab that allows you to monitor your cluster.

#### Browser
Alternatively, you can access the dashboard directly. Print the dashboard link with:

```python
print(client.dashboard_link)
# Example output: http://127.0.0.1:8787/status
```

Then, forward the port to access it on your local machine. On a terminal from the remote HPC node run:
```
ssh -L 8787:127.0.0.1:8787 username@hpc-cluster-address
```

Once connected via SSH with port forwarding, open a browser on your local machine and go to:
```
http://localhost:8787/status
```

## Additional Resources

- [Dask-JobQueue Documentation](https://jobqueue.dask.org/)
- [Dask on HPC Documentation](https://docs.dask.org/en/stable/deploying-hpc.html)
- [Dask Distributed Documentation](https://distributed.dask.org/)
- [JupyerLab Dask extension](https://jobqueue.dask.org/en/latest/clusters-interactive.html)

