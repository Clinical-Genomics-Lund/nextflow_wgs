## Nextflow WGS pipeline

## Deploy
To run pipeline first deploy to hopper using scp.

The most important files to deploy are:
* main.nf 
* configs/nextflow.hopper.config (to be named nextflow.config and be located in the same folder as main.nf)
* shards.csv/shards_38.csv
* git.hash
* bin-folder

All these files a deployed through deploy_hopper.sh/deploy_hopper_dev.sh


## Executing
Make sure nextflow, Java and singularity is loaded on hopper

nextflow run main.nf --csv path/to/csv.csv -profile see-below

## Profiles
The current pipeline can be run with three different panels. These are chosen by adding -profile to the nextflow run command. These are the three and their corresponding profiles:
* whole-genome sequencing -profile wgs
* oncogenetics -profile onco
* exome -profile exome
