## Distributed Locus Collector, Deduplication, Base Recalibration Scoring and Variant calling with dnascope

To speed up bam post alignment quality processing we use Sentieon distributed mode. This again splits the genome into smaller parts. We split the genome into 5 equally sized shards. The shards are defined in nextflow.config under params as two parameters; genomic_shards_file = "$baseDir/shards_5_38.csv" and genomic_shards_num = 5. The first is a file located in the main repository, and has to be in the same folder as main.nf. The basic process consists of 4 steps (see image below):


### LocusCollector

An algorithm that collects read information used later by the deduplication algorithm.
* input - full BAM
* output - duplication likelihoods for all shards

`sentieon driver --temp_dir /local/scratch/ -t ${task.cpus(threads)} -i $bam(full bam) $shard(1-5) --algo LocusCollector --fun score_info ${shard_name(1-5)}_${id(id of sample)}.score`

### Deduplication

Uses above scores to perform deduplication (removes duplicated reads in bam)
* input - full BAM + scores per shard
* output - dedup bam + dedup metrics, both per shard

`sentieon driver --temp_dir /local/scratch/ -t ${task.cpus} -i $bam $shard --algo Dedup --score_info $scores --metrics ${shard_name}_${id}_dedup_metrics.txt --rmdup ${shard_name}_${id}.bam`

In a separate process, the 5 shards are merged into a final deduped bam. This is then used by all callers other than dnascope. Dedup metrics are merged and used for [post-alignment QC](quality.md) presentation.

### Base recalibration score

Calculates BQSR-tables which is needed by dnascope.
* input - bam of present shard and its neighbouring shards. (ie shard 1 needs bam from 1 and 2, shard 2 needs bam from 1 2 and 3, and so on. )
* references - Mills_and_1000G_gold_standard.indels.hg38.vcf.gz and genome fasta used in alignment
* output - bqsr tables per shard

`sentieon driver --temp_dir /local/scratch/ -t ${task.cpus} -r $genome_file -i $bam_neigh $shard --algo QualCal -k $params.KNOWN ${shard_name}_${id}.bqsr.table`

In a separate process the bqsr-table shards are merged into a single bqsr table. This is used both for dnascope in distributed mode and for dnascope run in non-distributed mode.

### dnascope, snv calling

Single nucleotide variant and small insertions and deletions caller. Uses machine learning to make accurate calls
* input - deduped bam shard (with neighbours as in bqsr) and full bqsr-table
* references - genome fasta used in alignment
* output - GVCF per shard

`sentieon driver --temp_dir /local/scratch/ -t ${task.cpus} -r $genome_file -i $bam_neigh $shard -q $bqsr  --algo DNAscope --emit_mode GVCF ${shard_name}_${id}.gvcf.gz`

The GVCF shards are merged per sample, and then merged per family into a multisample vcf ready for normalization annotation steps.

### Mutect2, mitochondria calling

Due to limitations in heteroplasmid calls in dnascope we use an additional caller only for mitochondrial calling. GATK mutect2.
* input - mitochondria intersected deduped bam
* references - genome fasta for mitochondria (subset of full used for alignment)
* output - multisample/single sample vcf

This is merged, using picard, with the dnascope vcf after some specific mitochondrial [annotations](annotation.md)

![sentieonhomepage](https://support.sentieon.com/appnotes/_images/distributed_mode-fig3-2.png)