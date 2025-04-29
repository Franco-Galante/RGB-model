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