## Nextflow WGS pipeline

Still under development!

## Deploy
To run pipeline first deploy to hopper using scp or the nfs-mount(although not reliable atm)
Use a new folder on hopper to avoid deploying and changing another person's current work.

The most important files to deploy are:
* main.nf 
* configs/nextflow.hopper.config (to be named nextflow.config and be located in the same folder as main.nf)
* shards.csv/shards_38.csv

To run nextflow with singularity on hopper you need first module load.

`module load Java singularity nextflow`

the run command should look like this.

`nextflow run main.nf --csv input.csv -with-singularity /fs1/resources/containers/wgs_active.sif -resume -with-trace name.trace`

the script has a few options
--align (include if you want align)
    --shardbwa (sharded bwa, if running one sample that needs to be done fast)
--varcall (include if you want variant calling)
--annotate (include if you want annotations)
--gatkcov (include if you want coverage reports for sv-tools)


## CSV-file
needs to have the following headers
`clarity_sample_id,id,assay,sex,diagnosis,phenotype,group,father,mother,clarity_pool_id,platform,read1,read2,analysis_dir`

read1 and read2 can be fq.gz files or bam or vcf. Depending on which used different parts of pipeline will be called.
* fq.gz will start at mapping and then proceed to variant calling if flag is present and so on
* bam will start att variant calling and then proceed to annotation
* vcf will start at annotation

example-files can be located /fs1/bjorn/wgs_validation_csv/*.csv


## Known Issues

process bqsr sometimes fail, without any obvious reason, rerun pipeline with -resume. Possible solution, add errorstrategy to nextflow process.
