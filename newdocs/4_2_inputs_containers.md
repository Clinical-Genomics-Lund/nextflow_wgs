Several containers are used within the wgs pipeline. There is a "default" container used for many processes (defined in `nextflow.config`), and some more narrow containers used for specific pipelines.

The base container can be built by navigating into the `container` directory and running the following command:

```
sudo -E singularity build name_of_your_container.sif Singularity
```

The `-E` flag means that environment variables will be retained also within the singularity container (FIXME: Why do we want this? Does this assume that certain environment variables are set, which for instance are present on Lennart?)

Beyond the base container, a set of other containers are also used.

##### Madeline2

`madeline2` is installed separately. The Singularity recipe for building the dependencies is also present within the `container` folder in the repository.

```
sudo -E singularity build madeline2.sif Singularity_madeline2
```

##### VEP

The VEP container is built from a publicly available docker container. Currently, the pipeline uses version 103.

```
singularity pull docker://ensemblorg/ensembl-vep:release_103
```

When running this container, additional reference files needs to be provided. These are outlined in the [[4.2 Inputs - Annotation files]] page.

##### CADD

The CADD container is hosted in a separate repository: `git clone https://github.com/Clinical-Genomics-Lund/CADD-container/tree/v1.6`

In addition, reference files needs to be downloaded. Note that the file is really large - 200GB zipped. Use `wget` with the `-c` flag such to allow continuing the download if it is interrupted.

* `wget -c https://krishna.gs.washington.edu/download/CADD/v1.6/GRCh38/annotationsGRCh38_v1.6.tar.gz`
* `tar -xzvf annotationsGRCh38_v1.6.tar.gz`

Make sure that the lines in the `Singularity` recipe matches to where you have saved the annotations folder.

`ln -sr /fs1/resources/ref/hg38/annotation_dbs/CADD_v1.6/data/annotations/GRCh38_v1.6 GRCh38_v1.6`

Finally you build the container:

```
sudo -E singularity build cadd_v1.6.sif Singularity
```

##### GATK

This container is used for CNV-calling. It is retrieved from a Docker container provided by Broad Institute.

`singularity pull docker://broadinstitute/gatk`

##### Sentieon

The Sentieon container can be retrieved from the Galaxy Project: https://depot.galaxyproject.org/singularity/

This pipeline currently uses version 202112.

* `sentieon_202112.sif`
##### Stranger

* `stranger_0.8.sif`

FIXME: Where is this container retrieved from?

##### Twist myeloid

* `twistmyeloid_active.sif`

FIXME: Where is this container retrieved from?

##### AnnotSV

A repository containing the recipe to build the AnnotSV container can be retrieved as such:

```
git clone https://github.com/Clinical-Genomics-Lund/annotsv_container.git
```

To build the container, `cd` into the repository and run the following command:

```
sudo singularity build annotsv.v2.3.sif Singularity
```
##### POD

FIXME

`/fs1/resources/containers/POD_2020-05-19.sif`

Used in `plot_pod` process.
##### Genmod

`/fs1/resources/containers/genmod.sif`

##### Reviewer

Used in `reviewer` processs

`/fs1/resources/containers/REViewer_2021-06-07.sif`

##### Melt

`/fs1/resources/containers/melt_2.2.2.sif`
