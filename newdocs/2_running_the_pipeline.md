##### Dependencies

In order to run the pipeline, you will need to have [Nextflow](https://www.nextflow.io/) installed. At the moment, the pipeline is implemented in DSL-1 - an older version of Nextflow. This is no longer supported by the latest version of Nextflow. 

To run an older version (21 or earlier is needed), you can follow the instructions outlined [here](https://www.nextflow.io/docs/latest/getstarted.html) and run it as such:

```
NXF_VER=20.04.0 nextflow run hello
```

You will also need to have [Singularity](https://docs.sylabs.io/guides/latest/user-guide/) available. Singularity is a tool similar to Docker, which allows running specific sets of dependencies in so-called containers. This allows retaining frozen environments for individual processing steps.

If you are running this on the cluster in Lund (or a similar environment), these dependencies can be loaded as such:

```
module load Java
module load nextflow
module load singularity
```

##### Getting the source code

To get the source code, simply clone it from GitHub.

```
git clone https://github.com/Clinical-Genomics-Lund/nextflow_wgs.git
```

To actually run it, only a subset of the repository is needed. This is further discussed in the [deployment](3_how_to_deploy.md) section.

##### Running the pipeline

To run a so-called "stub-run", you can navigate into the repository and execute:

```
nextflow run main.nf -stub-run
```

This will run through the workflow only producing minimal files for each step. This is a good way to test-run the pipeline, and should complete in a matter of minutes.

To run real data, you will typically specify the input `--csv` file (further discussed in the section on [input files](4_1_input_files.md)). Furthermore, if rerunning a previous run, you will typically add the `-resume` flag to not recreate files unnecessarily.

Note that the difference between single dash arguments (`-argument`) and double dash arguments (`--argument`) is that single-dashed arguments are provided directly to Nextflow (`-resume`, `-stub-run`, `-profile`) while double-dashed arguments are provided as params to the script itself (`--csv`, `--resultsdir`).

If you have prepared the input data and is ready to start a run using the `wgs` profile, then you would execute the workflow as such:

```
nextflow run main.nf \
    -profile wgs \
    --csv path/to/input.csv \
    -resume
```

Additional useful arguments:

* `-w` To control the location of the work folder (default is the same folder as `main.nf`)
* `--resultsdir` Location of the output dir (FIXME: Is this correct?)

##### Running as a batch script

If you are running on a cluster environment, you will typically execute this using SLURM. 

If you happen to work at CMD in Lund, there is a dedicated system for this.

If not, below is an example SLURM script:

```
# SBATCH ...
# FIXME: First lines

nextflow run main.nf \
    -profile wgs \
    --csv path/to/input.csv \
    -resume
```

