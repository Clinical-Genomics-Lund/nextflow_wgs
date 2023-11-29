## Nextflow WGS pipeline

## Deploy
To run pipeline first deploy to hopper using scp.

The most important files to deploy are:
* `main.nf`
* `configs/nextflow.hopper.config` (to be named nextflow.config and be located in the same folder as main.nf)
* `shards.csv` / `shards_38.csv`
* `git.hash`
* `bin-folder`

All these files a deployed through deploy_hopper.sh/deploy_hopper_dev.sh


## Executing
Make sure nextflow, Java and singularity is loaded on hopper

```
nextflow run main.nf --csv path/to/csv.csv -profile see-below
```

## Profiles
The current pipeline can be run with six different panels. These are chosen by adding -profile to the nextflow run command. These are the three and their corresponding profiles:
* whole-genome sequencing `-profile wgs`
* oncogenetics `-profile onco` or `-profile oncov1`
* exome `-profile exome`
* other panels `-profile myeloid` or `-profile modycf`

# NextflowWgs

### OS-requirements

This is a pipeline written in Nextflow. Versions prior to 3.0 (not released) only works with Nextflow v19.04 or earlier. This is due to some deprecated nomenclature. Besides Nextflow, it requires Java (to run nextflow) and singularity to execute all software. 

### Installing the pipeline

Clone the repo. 

### Containers

#### Main

In the container directory install the main container, containing most software. 

`sudo -E singularity build name_of_your_container.sif Singularity`

Used when starting nextflow. 

`nextflow run main.nf -profile wgs -with-singularity name_of_your_container.sif`

#### Madeline2

Due to dependency issues madeline2 is installed separately:

`sudo -E singularity build madeline2.sif Singularity_madeline2`

#### VEP

VEP-container is built by 

`singularity pull docker://ensemblorg/ensembl-vep:release_103`. 

This container also needs VEP-references.

* VEP cache
* VEP-fasta
* gnomAD (exomes, genomes, mitochondria) VCFs
* phyloP
* phastcons
* CADD
* maxentscan plugin. 

Required by 3 VEP annotation steps, edit main.nf. use ctrl + f `ensembl-vep_release_103.sif` and replace with your container

#### CADD

`git clone https://github.com/Clinical-Genomics-Lund/CADD-container/tree/v1.6` 

Download all references

GRCh38 / hg38 (really big file! 200GB zipped)

`wget -c https://krishna.gs.washington.edu/download/CADD/v1.6/GRCh38/annotationsGRCh38_v1.6.tar.gz`

`tar -xzvf annotationsGRCh38_v1.6.tar.gz` 

Make sure to edit the following lines in the file `Singularity` to match where you saved the annotations-folder 

`ln -sr /fs1/resources/ref/hg38/annotation_dbs/CADD_v1.6/data/annotations/GRCh38_v1.6 GRCh38_v1.6` `sudo -E singularity build cadd_v1.6.sif Singularity` 

Required by CADD indels, edit main.nf. use ctrl + f `cadd_v1.6.sif` and replace with your container 

#### GATK 4.19 

`singularity pull docker://broadinstitute/gatk`

Required by CNV-calling, edit main.nf. use ctrl + f `gatk_4.1.9.0.sif` and replace with your container




