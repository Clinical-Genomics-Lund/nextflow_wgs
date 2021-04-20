## NextflowWgs

OS-requirements
This is a pipeline written in Nextflow. Versions prior to 3.0 (not released) only works with Nextflow v19.04 or earlier. This is due to some deprecated nomenclature. Besides Nextflow, it requires Java (to run nextflow) and singularity to execute all software. 

### Contents
  * [Installing](installing.md)
  * [nextflow.config](config.md)
  * [References](references.md)
  * [Input-csv and meta-information](input_meta_csv.md)
  * Pipeline content
    * Alignment
      * Sharded
      * Normal
    * Distributed dedup bqsr and SNV variant calling
    * Annotation
    * SNP-analyses
      * UDP
      * Peddy
      * Madeline2
    * Mitochondria calling
    * Short tandem repeats calling
    * Structural variation calling
    * Scout-yaml creation

