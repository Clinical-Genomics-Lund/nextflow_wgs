### Getting the source code

Clone the source code from GitHub.

```
git clone https://github.com/Clinical-Genomics-Lund/nextflow_wgs.git
```

### Dependencies

#### Nextflow

In order to run the pipeline, you will need to have [Nextflow](https://www.nextflow.io/) installed. At the moment, the pipeline is implemented using the DSL1 syntax. This is replaced by DSL2-syntax in recent versions of Nextflow.

To run this pipeline, you thus need to use an older version of Nextflow (latest version 21). As mentioned on the [Nextflow page](https://www.nextflow.io/docs/latest/getstarted.html), it is possible to run an older version as such:

```
NXF_VER=20.04.0 nextflow run hello
```

#### Singularity

[Singularity](https://docs.sylabs.io/guides/latest/user-guide/) is used to containerize dependencies for processes. This means that if you have the containers, you can run the workflow on any computer without requiring any further installations.

#### Running on a cluster

If you are running this on a computational cluster, the required dependencies can be loaded as such (if working on Hopper in Lund, or a computational cluster with nextflow and singularity set up).

```
module load Java
module load nextflow
module load singularity
```

### Running the pipeline

#### Stub-runs

Stub-runs allow you to test-run the pipeline without actually performing the full analysis. This creates dummy files for each processing steps, and normally completes in a matter of minutes on any computer. This is a useful way to test whether you have the pipeline set up correctly.

```
nextflow run main.nf -stub-run
```

#### Running real data

The following command can be used to run a full dataset, assuming that the input CSV is properly set up (in this case in the path `/path/to/input.csv`), that the `nextflow.config` file is configured (see more in the [input files section](input_files.md)) and that you are running using the `wgs` profile.

```
nextflow run main.nf \
    -profile wgs \
    --csv path/to/input.csv
```

Note that the difference between single dash arguments (`-profile`) and double dash arguments (`--csv`) is that single-dashed arguments are provided directly to Nextflow while double-dashed arguments are provided as params to the script itself, which overrides any options specified in the `nextflow.config` file.

Additional useful arguments:

* `-resume` If possible, continue a previous run and only rerun required steps.
* `-w` To control the location of the work folder (default is the same folder as `main.nf`)

### Running as a batch script

Note that if you are running jobs on Hopper in CMD, Lund, then the recommended way to execute jobs is using the `start_nextflow_analysis.pl` script.

Otherwise you will typically execute this using SLURM. A minimal example of a SLURM run script is shown below.

This can be executed by running `sbatch jobfile.run`.

```
# SBATCH --job-name=job_name
# SBATCH --output=slurm_log_%j.log
# SBATCH --ntasks=2
# SBATCH --mem=4gb
# SBATCH --time=2-00:00:00

module load Java
module load nextflow
module load singularity

nextflow run main.nf \
    -profile wgs \
    --csv path/to/input.csv
```

