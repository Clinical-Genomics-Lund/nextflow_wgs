#!/usr/bin/env nextflow

genome_file = params.genome_file
K_size      = 100000000
KNOWN1 = params.refpath+"/annotation_dbs/1000G_phase1.indels.b37.vcf.gz"
KNOWN2 = params.refpath+"/annotation_dbs/Mills_and_1000G_gold_standard.indels.b37.vcf.gz"
sentieon_model = params.refpath+"/sw/sentieon/SentieonDNAscopeModelBeta0.4a-201808.05.model";
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

//bwa_in = bwa_shards.combine(fastq)




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


merged_bam.into{merged_bam1; merged_bam2}

process locus_collector {
    cpus 16
    //publishDir '/fs1/results/bam/wgs/scores'
    input:
    set val(shard_name), val(shard) from shards1
    each file(bam) from merged_bam1
    

    output:
    set val(id), file("${id}_${shard_name}.score"), file("${id}_${shard_name}.score.idx") into locus_collector_scores
    set val(id), file(bams), file(bai) into tester
    script:
    // merged_bam1 recieves both bam and bai, remove the bai from input
    bams = bam[0]
    bai = bam[1]
    id = bams.getBaseName().tokenize("_")[0]
    group = "scores"
    """
    sentieon driver -t ${task.cpus} -i $bams $shard --algo LocusCollector --fun score_info ${id}_${shard_name}.score
    """
}


 locus_collector_scores
     .groupTuple()
     .set{ all_scores }


process dedup {
    cpus 16

    input:
    set val(id), file(score), file(idx), val(shard_name), val(shard), file(bam), file(bai) from all_scores.combine(shards2).join(tester)

    output:
    set val(id), file("${shard_name}_${id}.bam"), file("${shard_name}_${id}.bam.bai") into shard_dedup_bam

    script:
    scores = score.join(' --score_info ')
    bam_group = "bams"
    """
    sentieon driver -t ${task.cpus} -i $bam $shard --algo Dedup --score_info $scores --rmdup ${shard_name}_${id}.bam
    """
    
}

shard_dedup_bam
    .groupTuple()
    .into{ all_dedup_bams1; all_dedup_bams2;  }

shards3
    .merge(tuple(shardie1))
    .into{ shard_shard }

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
    //file(tables) from bqsr_table.collect()
    set id, file(tables) from bqsr_table.groupTuple()
    output:
    file("${id}_merged.bqsr.table") into bqsr_merged
    """
    sentieon driver --passthru --algo QualCal --merge ${id}_merged.bqsr.table $tables
    """
}

// process dnascope {
//     cpus 16

//     input:
//     set val(shard_name), val(shard), val(group), file(bams), file(bai) from shards4.combine(all_dedup_bams2)
//     val(combo) from shardie2
//     each file(bqsr) from bqsr_merged
//     output:
//     file("${shard_name}.dnascope.vcf") into vcf_shard
//     file("${shard_name}.dnascope.vcf.idx") into vcf_idx
//     script:
//     combo = (combo - 0) //first dummy value
//     combo = (combo - 6) //last dummy value
//     commons = (combo.collect{ "${it}.bam" } - bams)   //add .bam to each shardie, remove all other bams
//     bam_neigh = commons.join(' -i ')
//     """
//     sentieon driver -t ${task.cpus} -r $genome_file -i $bam_neigh $shard --algo DNAscope --model $sentieon_model ${shard_name}.vcf.tmp
//     sentieon driver -t ${task.cpus} -r $genome_file $shard --algo DNAModelApply --model $sentieon_model -v ${shard_name}.vcf.tmp ${shard_name}.dnascope.vcf
//     """
// }

// process merge_vcf {

//     input:
// 	file(vcfs) from vcf_shard.collect()
//         file(idx) from vcf_idx.collect()
//     output:
// 	file("${name}.dnascope.vcf") into complete_vcf

//     script:
//         vcfs_sorted = vcfs.sort(false) { a, b -> a.getBaseName().tokenize(".")[0] as Integer <=> b.getBaseName().tokenize(".")[0] as Integer } .join(' ')
    
//     """
//     sentieon driver --passthru --algo DNAscope --merge ${name}.dnascope.vcf $vcfs_sorted
//     """
// }
