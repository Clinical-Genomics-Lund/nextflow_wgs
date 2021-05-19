## Annotation and normalization

Several annotation steps are done after variant calling. Here all annotations of SNVs are documentated. If not specified, all SNV annotation is relevant for both mitochondrial and nuclear variants.

### Mitochondrial variants

#### split/normalize/filter
* input - mutect2 vcf
  * vcfbreakmulti
  * bcftools sort
  * tabix
  * bcftools norm 
  * bcftools remove below 5% variants, 0 genotypes and variants missing in proband

#### hmtnote
* input - split/normalized/filtered vcf
  * fix info-field remove spaces with sed

#### haplogrep

* input mutect2 vcf
  * plot haplotypes, join trio images with montage

#### eklipse
* input - mito-bam
  * plot coverage in circle plots

### Single nucleotide variants

#### split/normalize/filter/intersect
* input - dnascope merged vcf
  * vcfbreakmulti
  * bcftools norm
  * bcftools sort
  * wgs_DPAF_filter.pl
  * bedtools intersect (coding + clinvar introns)
  * combine with mito variants with picard MergeVcfs
  * rename M -> MT (genmod requires)

#### Variant Effect Predictor
* input - split/normalized/filtered/intersected
  * --offline
  * --everything
  * --merged
  * --vcf
  * --no_stats
  * --synonyms $params.SYNONYMS
  * --force_overwrite
  * --plugin CADD,$params.CADD
  * --plugin LoFtool
  * --plugin MaxEntScan,$params.MAXENTSCAN,SWA,NCSS
  * --fasta $params.VEP_FASTA
  * --dir_cache $params.VEP_CACHE
  * --dir_plugins $params.VEP_CACHE/Plugins
  * --distance 200
  * -cache
  * -custom $params.GNOMAD_EXOMES,gnomADe,vcf,exact,0,AF_popmax,AF,popmax
  * -custom $params.GNOMAD_GENOMES,gnomADg,vcf,exact,0,AF_popmax,AF,popmax
  * -custom $params.GNOMAD_MT,gnomAD_mt,vcf,exact,0,AF_hom,AF_het
  * -custom $params.PHYLOP
  * -custom $params.PHASTCONS

#### vcfanno
* input - VEP vcf
  * clinvar variants
  * loqusdb variants (local observations)
  * gene annotations from ensembl

#### inheritance models
* input - vcfanno, pedigree from meta info
  * genmod models

#### modify vcf
* input - inheritance vcf
  * modify_vcf_scout.pl
  * adjusts, counts and presents information for ranking and visualization in scout

#### mark splice
* input - modified vcf
  * mark_spliceindels.pl
  * self-explanatory

#### indels special route
* input - split/normalize/filter/intersect vcf
  * bcftools extracts indels
  * minimal VEP annotation
  * calculate CADD score for indels in sample using CADD script
  * re-add indel CADD scores to other SNVs using genmod annotate

#### genmod score
* input - marksplice + cadd scores vcf + rankmodel
  * score according to rankmodel one for family and one for single, genmod score
  * calculate compounding variants with genmod compound for families

### Structural variants
* input 

### Short tandem repeats
