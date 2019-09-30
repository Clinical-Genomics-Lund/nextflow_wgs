## Nextflow WGS pipeline

Still under development!


To run pipeline first deploy to hopper using scp or the nfs-mount(although not reliable atm)
Use a new folder on hopper to avoid deploying and changing another person's current work.
The most important files to deploy are main.nf and configs/nextflow.hopper.config, the last to be name nextflow.config and be located in the same folder as main.nf

To run nextflow with singularity on hopper you need first module load.

`module load Java singulairty nextflow`

the run command should look like this.

`nextflow run main.nf --csv input.csv -with-singularity /fs1/resources/containers/YOURCONTAINER.sif -resume`

the script has a few options
--align (include if you want align)
--annotate (include if you want annotations)
