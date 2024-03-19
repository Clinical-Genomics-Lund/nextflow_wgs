## Short Tandem Repeat variant calling

Short tandem repeats (STR) is called by [expansionhunter](https://github.com/Illumina/ExpansionHunter) and annotated with [stranger](https://github.com/moonso/stranger)

### ExpansionHunter
* input - deduped BAM, variant catalogue
  * only run for proband!
  * also generate [graph alignment viewer](https://github.com/Illumina/GraphAlignmentViewer) images per sample

### Stranger
* input - expansionhunter VCF
  * fix INFO-field spaces with grep and sed

### vcfbreakmulti_expansionhunter
* input - EH + Stranger VCF
  * fix sample names of VCF with Picard RenameSampleInVcf
  * vcfbreakmulti, split multi-allele variants
  * if family, create a fake VCF with 0/0 parents with familify_str.pl