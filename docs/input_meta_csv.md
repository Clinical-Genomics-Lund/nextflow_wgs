## This is the CSV and meta data docs


### Required fields, can be in any order.

#### ID headers
* id (ID of sample)
* clarity_sample_id (could be the same as Sample-ID, used for scout-db matches between clarity-lims)
* group (group/family name. Changed when samples are rerun)
* mother (id of mother, empty if single or not proband)
* father (id of father, empty if single or not proband)

#### Sample-information
* type (proband/relative)
* sex (F/M)
* phenotype (affected/unaffected)
* diagnosis (chosen gene-panel for sample)
* read1 (fastq read1) *
* read2 (fastq read2) *
* assay (what analysis, automated choice through start_nextflow_analysis.pl -profile wgs, and choses right institute in scout for sample)
* analysis (subanalysis, ie wgsdev would put samples in validation institute instead of clinical)

"*" read1 and read2 can be supplemented bam and bqsr for the samples to bypass alignment. Also somewhat works for vcfs (vcf + idx) for SNV annotation

### Other useful fields

* priority (changes the priority of sample in slurm, default is defined in nextflow.config)

#### Example CSV-file

[demo-csv](demo/1999-20.csv)
