#!/usr/bin/env nextflow

genome_file = params.genome_file
K_size      = 100000000
KNOWN1 = params.refpath+"/annotation_dbs/1000G_phase1.indels.b37.vcf.gz"
KNOWN2 = params.refpath+"/annotation_dbs/Mills_and_1000G_gold_standard.indels.b37.vcf.gz"
sentieon_model = params.sentieon_model
bwa_num_shards = params.bwa_shards
bwa_shards = Channel.from( 0..bwa_num_shards-1 )

genomic_num_shards = params.genomic_shards_num

csv = file(params.csv)
mode = csv.countLines() > 2 ? "family" : "single"
println(mode)

Channel
    .fromPath(params.csv)
    .splitCsv(header:true)
    .map{ row-> tuple(row.group, row.id, row.read1, row.read2) }
    .set { fastq }

Channel
    .fromPath(params.csv)
    .splitCsv(header:true)
    .map{ row-> tuple(row.group, row.id, row.sex, row.mother, row.father, row.phenotype, row.diagnosis) }
    .into { ped; all_ids; yml_diag }



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
    .from( 0..(genomic_num_shards+1) )
    .collate( 3,1, false )
    .into{ shardie1; shardie2 }






process bwa_align {
    cpus 56
    memory '64 GB'

    input:
	set val(shard), val(group), val(id), r1, r2 from bwa_shards.combine(fastq)

    output:
	set val(id), val(group), file("${id}_${shard}.bwa.sort.bam"), file("${id}_${shard}.bwa.sort.bam.bai") into bwa_shards_ch

    """
    sentieon bwa mem -M -R '@RG\\tID:${group}_${id}\\tSM:${group}_${id}\\tPL:illumina' -K $K_size -t ${task.cpus} -p $genome_file '<sentieon fqidx extract -F $shard/$bwa_num_shards -K $K_size $r1 $r2' | sentieon util sort -r $genome_file -o ${id}_${shard}.bwa.sort.bam -t ${task.cpus} --sam2bam -i -
    """
}

process bwa_merge_shards {
    cpus 56
    publishDir '/fs1/results/bam/wgs'

    input:
        set val(id), val(group), file(shard), file(shard_bai) from bwa_shards_ch.groupTuple()

    output:
        set file("${id}_merged.bam"), file("${id}_merged.bam.bai") into merged_bam
        
    script:
        bams = shard.sort(false) { a, b -> a.getBaseName() <=> b.getBaseName() } .join(' ')

    """
    sentieon util merge -o ${id}_merged.bam ${bams}
    """
}


process locus_collector {
    cpus 16
    input:
    set file(bam), file(bai), val(shard_name), val(shard) from merged_bam.combine(shards1)
    //each file(bam) from merged_sort_bam
    
    output:
    set val(id), file("${shard_name}.score"), file("${shard_name}.score.idx") into locus_collector_scores
    set val(id), file(bam), file(bai) into tester
    script:
    // merged_bam1 recieves both bam and bai, remove the bai from input
    id = bam.getBaseName().tokenize("_")[0]
    group = "scores"
    """
    sentieon driver -t ${task.cpus} -i $bam $shard --algo LocusCollector --fun score_info ${shard_name}.score
    """
}


locus_collector_scores
    .groupTuple()
    .join(tester)
    .combine(shards2)
    .set{ all_scores }


process dedup {
    cpus 16

    input:
    set val(id), file(score), file(idx), file(bam), file(bai), val(shard_name), val(shard) from all_scores

    output:
    set val(id), file("${shard_name}_${id}.bam"), file("${shard_name}_${id}.bam.bai") into shard_dedup_bam

    script:
    scores = score.sort(false) { a, b -> a.getBaseName().tokenize(".")[0] as Integer <=> b.getBaseName().tokenize(".")[0] as Integer } .join(' --score_info')
    //scores = score.join(' --score_info ')
    bam_group = "bams"
    """
    sentieon driver -t ${task.cpus} -i $bam $shard --algo Dedup --score_info $scores --rmdup ${shard_name}_${id}.bam
    """
    
}

shard_dedup_bam
    .groupTuple()
    .into{ all_dedup_bams1; all_dedup_bams2;  }
//merge shards with shard combinations
shards3
    .merge(tuple(shardie1))
    .into{ shard_shard; shard_shard2 }

process bqsr {
    cpus 16

    input:
    set val(id), file(bams), file(bai), val(shard_name), val(shard), val(one), val(two), val(three) from all_dedup_bams1.combine(shard_shard)
    output:
    set val(id), file("${shard_name}_${id}.bqsr.table") into bqsr_table
    script:
    combo = [one, two, three]
    combo = (combo - 0) //first dummy value
    combo = (combo - (genomic_num_shards+1)) //last dummy value
    commons = (combo.collect{ "${it}_${id}.bam" } - bams)   //add .bam to each shardie, remove all other bams
    bam_neigh = commons.join(' -i ')
    """
    sentieon driver -t ${task.cpus} -r $genome_file -i $bam_neigh $shard --algo QualCal -k $KNOWN1 -k $KNOWN2 ${shard_name}_${id}.bqsr.table
    """
}


process merge_bqsr {
    publishDir '/fs1/results/bam/wgs/bqsr_tables'
    input:
    set id, file(tables) from bqsr_table.groupTuple()
    output:
    set val(id), file("${id}_merged.bqsr.table") into bqsr_merged
    """
    sentieon driver --passthru --algo QualCal --merge ${id}_merged.bqsr.table $tables
    """
}

all_dedup_bams2
    .join(bqsr_merged)
    .set{ all_dedup_bams3 }


all_dedup_bams3
    .combine(shard_shard2)
    .set{ bam_shard_shard }

process dnascope {
    cpus 16

    input:
    set id, file(bams), file(bai), file(bqsr), val(shard_name), val(shard), val(one), val(two), val(three) from bam_shard_shard
    output:
    set id, file("${shard_name}.vcf"), file("${shard_name}.vcf.idx") into vcf_shard
    script:
    combo = [one, two, three]
    combo = (combo - 0) //first dummy value
    combo = (combo - (genomic_num_shards+1)) //last dummy value
    commons = (combo.collect{ "${it}_${id}.bam" })   //add .bam to each shardie, remove all other bams
    bam_neigh = commons.join(' -i ')
    type = mode == "family" ? "--emit_mode GVCF" : ""
    """
    sentieon driver -t ${task.cpus} -r $genome_file -i $bam_neigh $shard -q $bqsr --algo DNAscope $type ${shard_name}.vcf
    """
}

process merge_vcf {

    input:
	set id, file(vcfs), file(idx) from vcf_shard.groupTuple()
    output:
	file("${id}.dnascope.vcf") into complete_vcf

    script:
    vcfs_sorted = vcfs.sort(false) { a, b -> a.getBaseName().tokenize(".")[0] as Integer <=> b.getBaseName().tokenize(".")[0] as Integer } .join(' ')
    
    """
    sentieon driver --passthru --algo DNAscope --merge ${id}.dnascope.vcf $vcfs_sorted
    """
}

