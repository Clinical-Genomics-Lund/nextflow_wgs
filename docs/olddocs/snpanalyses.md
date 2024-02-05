## SNP analyses

A variety of processes creates non-vcf outputs for visualization and database management of variants.


### add to loqusdb
* input - intersected SNV VCF
  * creates a load command for loading sample variants into local observation database

### peddy
* input - SNV VCF + pedigree
  * quality control
    * sex
    * parent matching
    * ancestry

### fastgnomad
* input intersected SNV VCF
  * software annotate
  * selects common SNPs from sample
  * used for various downstream analyses

### uniparental disomy (upd)
* input - fastgnomad VCF
  * for trios only
  * finds stretches of uniparental disomies

### upd table
* input - data from upd calculations
  * creates a table with informative SNPs as support for plots

### regions of heterozygozity (roh)
* input - fastgnomad VCF
  * bcftools roh --rec-rate 1e-9 --AF-tag GNOMADAF

### gatkcov
* input - deduped BAM
  * gatk CollectReadCounts, DenoiseReadCounts and PlotDenoisedCopyRatios 
  * generates inputs for various downstream analyses

### overview_plot
* input - updplot, roh_plot and gatk_plot
  * creates a plot with all three types of data

### generate gens data
* input - gatk cov data
  * generate_gens_data.pl, creates input for GENS software. Bins of cov depths over whole genome

### svvcf to bed
* input - SV VCF
  * creates a BED file from CNV segments, using cnv2bed.pl

### plot pod (parental origin of duplication)
* input - fastgnomad VCF and SV VCF, pedigree file and metainfo
  * parental_origin_of_duplication.pl, requires specific container. Not documented!! 
