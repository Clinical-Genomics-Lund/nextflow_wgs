## NextflowWgs

OS-requirements
This is a pipeline written in Nextflow. Versions prior to 3.0 (not released) only works with Nextflow v19.04 or earlier. This is due to some deprecated nomenclature. Besides Nextflow, it requires Java (to run nextflow) and singularity to execute all software. 

### Contents
  * [Installing](installing.md)
  * [nextflow.config](config.md)
  * [References](references.md)
  * [Input-csv and meta-information](input_meta_csv.md)
  * Pipeline content
    * [Alignment](alignment.md)
      * Sharded
      * Normal
    * [Post alignment Quality Control](quality.md)
    * [Distributed dedup bqsr and SNV variant calling](snv-calling.md)
      * [Annotation](annotation.md)
    * SNP-analyses
      * UDP
      * Peddy
      * Madeline2
    * [Short tandem repeats calling](str-calling.md)
    * [Structural variation calling](sv-calling.md)
      * [Annotation](svannotation.md)
    * Scout-yaml creation

