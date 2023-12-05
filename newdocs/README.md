The wgs pipeline is used both to process full genome sequencing, exomes and panels. It is designed to find rare variants. The expected input data are FASTQ files with raw sequencing data, but can also work with pre-aligned data (BAM-files) or perform annotation only (VCF-files).

The output from this workflow is closely related to the variant selection software [Scout](https://github.com/Clinical-Genomics/scout), an open source software for identifying variants with clinical significance. It has been developed on the local computational cluster in Lund, CMD, and is primarily designed to fulfill the needs of CMD. This it will likely not run "out of the box" in a different environment.

The pipeline uses singularity containers to allow executing software with various dependencies without the need to install these on the server. It is designed to work with the resource management system SLURM.

If you run into any issues using this pipeline, feel free to open an issue in the GitHub repository, or reach out to the responsible bioinformaticians.

### Table of content

* [Running the pipeline](2_running_the_pipeline.md)
    * Profiles
    * Singles and trios
* [How to deploy](3_how_to_deploy.md)
* [Inputs](4_0_inputs.md)
    * [Input files](4_1_input_files.md)
    * [Annotation files](4_2_annotation_files.md)
    * [Containers](4_3_inputs_containers.md)
    * [Sentieon license](4_4_inputs_sentieon_license.md)
* [Outputs](5_outputs.md)
* [List of used software](6_list_of_all_used_software.md)

##### Overview and the profiles

Below you see an overview of the workflow and the different run modes. 

The five profiles are:

* wgs
* onco / oncov1
* exome
* myeloid
* mobycf

The typical workflow consists of the following steps:

1. Optionally preprocess, align and dedup reads
2. Calculate QC and coverage information
3. Perform SNV and SV calling using separate workflows
4. Generate annotations
5. Calculate "severity" scores and merge variants

In addition, different profiles will influence the execution of the workflow.

The `wgs` run mode will in addition run a mitochondria variant calling workflow. If run in a trio (family mode), it will in addition attempt to identify UPD and ROH. For structural variants, it will run a different set of callers adapted to call variants across the full genome (FIXME?), including a dedicated workflow to identify copy-number variants. If run in family mode, it will also attempt to use this information to find compound variants where variants from parents together influence both alleles of a certain gene.

For the `onco` run mode, FreeBayes will be used in addition to DNAScope to perform SNV calling (FIXME: why?).

For panels, a different set of callers is used for SV calling.

![overview_img](img/wgs_overview_200.drawio.png)

