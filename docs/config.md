## this is the documentation for nextflow.config

The nextflow.config file defines most of the things needed to run the pipeline. This includes paths to [references](references.md), singularity configuation, environmental configs, output-folders, profiles and executor (slurm) configurations.

Example [config](https://github.com/Clinical-Genomics-Lund/nextflow_wgs/blob/master/configs/nextflow.hopper.config). Do note that this file is setup for a very specific cluster in Lund. You will need to configure it for your specific enironment.

### Config headers

#### singularity

* enabled = true (if you are using containers)
* runOptions, what commands should be run with chosen container. --bind /fs1/ (our fileserver) and --bind /local/ (scratch-folders on our nodes) adds these paths to the container so that it may read/write to these paths.

#### params

All these parameters can be accessed by all other sections of config and by main.nf. Here we define all global parameters that will be used by all profiles (see below section). Can be accessed in main.nf by params.YOURPARAMETER. Set your reference-paths here and later in the profile sections.

#### process
* executor, what executor will be run. We use 'slurm', can use 'local' or most other executors (sge etc)
* queue, this is that priority or rather what slurm partition the sample will be run in. The default is defined in the params section. Use your own slurm partitions!!!
* time, default time per process in main.nf unless time is defined in process header
* container, what container will be run. Default is defined in params section. Set your own path to the [main container](installing.md)

#### profiles

This is the analysis specific section. References that are specific per analysis goes here, as well as output directories. Also, definitions of what processes that are to be run. Profiles are chosen by adding -profile [chosen profile] to the nextflow run command. 

Below are some of the important things needed to be set for new systems

* params.outdir, where published files go, added together with..
* params.subdir
* params.crondir, directory used for symlinking between university/regional computers (might not be needed for general user)
* params.intersect_bed, bedfile defining what regions of genome is of interest for SNV calling.
* params.svdb, [svdb](references.md) vcf for artefact filtering in SV vcfs
* params.rank_model_s, path to single SNV rank model
* params.rank_model, path to family SNV rank model
* params.svrank_model_s, path to single SNV rank model
* params.svrank_model, path to family SNV rank model
* params.vcfanno, path to definition file for vcfanno ADD TO REFERENCE.md!!!