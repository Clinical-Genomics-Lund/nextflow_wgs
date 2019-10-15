## Nextflow WGS pipeline

Still under development!

## Deploy
To run pipeline first deploy to hopper using scp or the nfs-mount(although not reliable atm)
Use a new folder on hopper to avoid deploying and changing another person's current work.

NEW!! Need to deploy bin-folder to runfolder on hopper!

The most important files to deploy are:
* main.nf 
* configs/nextflow.hopper.config (to be named nextflow.config and be located in the same folder as main.nf)
* shards.csv

To run nextflow with singularity on hopper you need first module load.

`module load Java singularity nextflow`

the run command should look like this.

`nextflow run main.nf --csv input.csv -with-singularity /fs1/resources/containers/wgs_active.sif -resume`

the script has a few options
--align (include if you want align)
--annotate (include if you want annotations)


## CSV-file
needs to have the following headers
`clarity_sample_id,id,assay,sex,diagnosis,phenotype,group,father,mother,clarity_pool_id,platform,read1,read2,analysis_dir`

In diagnosis field, use; validering, this will make sure it ends up in the wgsvalidering institute on SCOUT

example-files can be located /fs1/pipelines/wgs_germline/annotation/ any csv that is not shards.csv should work.


## Known Issues

process bqsr sometimes fail, without any obvious reason, rerun pipeline with -resume. Possible solution, add errorstrategy to nextflow process.
