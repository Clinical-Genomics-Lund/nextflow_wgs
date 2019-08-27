#!/usr/bin/env nextflow

genome_file = params.refpath+"/hg19/fasta/human_g1k_v37_decoy.fasta"
name        = params.name
K_size      = 100000000
KNOWN1 = params.refpath+"/hg19/annotation_dbs/1000G_phase1.indels.b37.vcf.gz"
KNOWN2 = params.refpath+"/hg19/annotation_dbs/Mills_and_1000G_gold_standard.indels.b37.vcf.gz"
sentieon_model = params.refpath+"/sw/sentieon/SentieonDNAscopeModelBeta0.4a-201808.05.model";
bwa_num_shards = params.bwa_shards
bwa_shards = Channel.from( 0..bwa_num_shards-1 )

genomic_num_shards = params.genomic_shards_num

// Check if paired or unpaired analysis
mode = "paired"
fastq = Channel.create()
if(params.fastq_N_R1 && params.fastq_N_R2) {
    fastq = [['T', file(params.fastq_T_R1), file(params.fastq_T_R2)],
	     ['N', file(params.fastq_N_R1), file(params.fastq_N_R2)]]
}
else {
    fastq = [['T', file(params.fastq_T_R1), file(params.fastq_T_R2)]]
    mode = "unpaired"
}


if(genome_file ){
    bwaId = Channel
            .fromPath("${genome_file}.bwt")
            .ifEmpty { exit 1, "BWA index not found: ${genome_file}.bwt" }
}


Channel
    .fromPath(params.genomic_shards_file)
    .splitCsv(header:false)
    .into { shards1; shards2; shards3; shards4; shards5; }

// A channel to pair neighbouring bams and vcfs. 0 and top value removed later
// Needs to be 0..n+1 where n is number of shards in shards.csv
Channel
    .from( 0..genomic_num_shards )
    .collate( 3,1, false )
    .into{ shardie1; shardie2 }


bwa_in = bwa_shards.combine(fastq)

process bwa_align {
    cpus 56
    memory '64 GB'
    time '24:00:00'

    input:
	set val(shard), val(type), file(r1), file(r2) from bwa_in

    output:
	set val(type), file("shard_${shard}.bwa.sort.bam"), file("shard_${shard}.bwa.sort.bam.bai") into bwa_shards_ch

    """
    sentieon bwa mem -M -R '@RG\\tID:${name}_${type}\\tSM:${name}_${type}\\tPL:illumina' -K $K_size -t ${task.cpus} -p $genome_file '<sentieon fqidx extract -F $shard/$bwa_num_shards -K $K_size $r1 $r2' | sentieon util sort -r $genome_file -o shard_${shard}.bwa.sort.bam -t ${task.cpus} --sam2bam -i -
    """
}


process bwa_merge_shards {
    input:
	set val(type), file(shard), file(shard_bai) from bwa_shards_ch.groupTuple()

    output:
    file("merged.bam") into merged_bam
    file("merged.bam.bai") into merged_bai
    script:
	bams = shard.sort(false).join(' ')
    """
    sentieon util merge -o merged.bam ${bams}
    """
}


merged_bam.into{merged_bam1; merged_bam2}
merged_bai.into{merged_bai1; merged_bai2}

process locus_collector {
    cpus 16

    input:
    set val(shard_name), val(shard) from shards1
    each file(bam) from merged_bam1
    each file(bai) from merged_bai1

    output:
    set val(group), file("${shard_name}.score.gz"), file("${shard_name}.score.gz.tbi") into locus_collector_scores
    
    script:
    group = "scores"
    """
    sentieon driver -t ${task.cpus} -i $bam $shard --algo LocusCollector --fun score_info ${shard_name}.score.gz
    """
}


locus_collector_scores
   .groupTuple()
   .set{ all_scores }


process dedup {
    cpus 16

    input:
    set val(group), file(score), file(tbi), val(shard_name), val(shard) from all_scores.combine(shards2)
    each file(bam) from merged_bam2
    each file(bai) from merged_bai2

    output:
    set val(bam_group), file("${shard_name}.bam"), file("${shard_name}.bam.bai") into shard_dedup_bam

    script:
    scores = score.join(' --score_info ')
    bam_group = "bams"
    """
    sentieon driver -t ${task.cpus} -i $bam $shard --algo Dedup --score_info $scores --rmdup ${shard_name}.bam
    """
}

shard_dedup_bam
    .groupTuple()
    .into{ all_dedup_bams1; all_dedup_bams2;  }



process bqsr {
    cpus 16

    input:
    set val(shard_name), val(shard), val(group), file(bams), file(bai) from shards3.combine(all_dedup_bams1)
    val(combo) from shardie1
    output:
    file("${shard_name}.bqsr.table") into bqsr_table
    script:
    combo = (combo - 0) //first dummy value
    combo = (combo - 6) //last dummy value
    commons = (combo.collect{ "${it}.bam" } - bams)   //add .bam to each shardie, remove all other bams
    bam_neigh = commons.join(' -i ')
    """
    sentieon driver -t ${task.cpus} -r $genome_file -i $bam_neigh $shard --algo QualCal -k $KNOWN1 -k $KNOWN2 ${shard_name}.bqsr.table
    """
}


process merge_bqsr {
    
    input:
    file(tables) from bqsr_table.collect()
    output:
    file("merged.bqsr.table") into bqsr_merged
    """
    sentieon driver --passthru --algo QualCal --merge merged.bqsr.table $tables
    """
}

process dnascope {
    cpus 16

    input:
    set val(shard_name), val(shard), val(group), file(bams), file(bai) from shards4.combine(all_dedup_bams2)
    val(combo) from shardie2
    each file(bqsr) from bqsr_merged
    output:
    file("${shard_name}.dnascope.vcf") into vcf_shard
    file("${shard_name}.dnascope.vcf.idx") into vcf_idx
    script:
    combo = (combo - 0) //first dummy value
    combo = (combo - 6) //last dummy value
    commons = (combo.collect{ "${it}.bam" } - bams)   //add .bam to each shardie, remove all other bams
    bam_neigh = commons.join(' -i ')
    """
    sentieon driver -t ${task.cpus} -r $genome_file -i $bam_neigh $shard --algo DNAscope --model $sentieon_model ${shard_name}.vcf.tmp
    sentieon driver -t ${task.cpus} -r $genome_file $shard --algo DNAModelApply --model $sentieon_model -v ${shard_name}.vcf.tmp ${shard_name}.dnascope.vcf
    """
}

process merge_vcf {

    input:
    file(vcfs) from vcf_shard.toSortedList()
    file(idx) from vcf_idx.collect()
    output:
    file("${name}.dnascope.vcf") into complete_vcf
    """
    sentieon driver --passthru --algo DNAscope --merge ${name}.dnascope.vcf $vcfs
    """
}
