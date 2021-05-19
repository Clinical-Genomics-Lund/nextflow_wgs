## Structural variant calling

Structural variant calling is done by three different callers. These are then merged by type (del, dup, ins etc) if they overlap by 70% or more. See sv-annotation for information of manipulation and annotation of vcfs.


### Manta
* input - deduped BAM
  * outputs diploid VCF with del dup ins and BNDs.
  * includes mitochondrial variants

### gatk
* input - deduped BAM
  * several steps
    * gatk CollectReadCounts (coverage)
    * gatk DetermineGermlineContigPloidy 
    * GermlineCNVCaller, subdivided into 8 parts
    * PostprocessGermlineCNVCalls joins 8 parts and outputs genotyped segments (del and dup)
    * filtered with filter_gatk.pl which removed 0 genotypes
    * mergeGATK.pl join neighbouring segments of the same type. large del or large dups

### tiddit
* input - deduped BAM
  * filtered output by variant quality = PASS

### svdb merge
* input - all three called SV-VCFs for all samples
  * merges 9 vcfs if trio
  * 70% overlap required
  * no within-sample merging
  * breakpoint priority from manta>tiddit>gatk
  * remove all BNDs, not used and causes VEP to fail if first variant is BND
  * merge_callsets.pl, joins callers into neatly named callers for scout to present per variant

  

  