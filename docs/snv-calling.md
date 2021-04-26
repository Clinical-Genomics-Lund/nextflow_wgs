## Distributed Locus Collector, Deduplication, Base Recalibration Scoring and Variant calling with dnascope

To speed up bam post alignment quality processing we use Sentieon distributed mode. This again splits the genome into smaller parts. We split the genome into 5 equally sized chunks. The chunks are defined in nextflow.config under params as two parameters; genomic_shards_file = "$baseDir/shards_5_38.csv" and genomic_shards_num = 5. The first is a file located in the main repository, and has to be in the same folder as main.nf. The basic process consists of 4 steps:


### LocusCollector

An algorithm that collects read information used later by the deduplication algorithm.
* input - BAM
* output - duplication likelihoods

`sentieon driver

--temp_dir /local/scratch/ \

-t ${task.cpus} \

-i $bam $shard \

--algo LocusCollector \

--fun score_info ${shard_name}_${id}.score`


![sentieonhomepage](https://support.sentieon.com/appnotes/_images/distributed_mode-fig3-2.png)