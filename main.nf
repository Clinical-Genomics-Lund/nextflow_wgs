#!/usr/bin/env nextflow

genome_file = file(params.fasta)
regions_bed = file(params.bed)
name        = params.name
K_size      = 100000000

bwa_num_shards = 2
bwa_shards = Channel.from( 0..bwa_num_shards-1 )

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


if(params.fasta ){
    bwaId = Channel
            .fromPath("${params.fasta}.bwt")
            .ifEmpty { exit 1, "BWA index not found: ${params.fasta}.bwt" }
}


// dedup_shards = Channel.create()
// dedup_shards = [["1", "--shard 1:1-249250621 --shard 2:1-243199373 --shard 3:1-198022430 --shard 4:1-191154276 --shard 5:1-118373300"],
// 	        ["2", "--shard 5:118373301-180915260 --shard 6:1-171115067 --shard 7:1-159138663 --shard 8:1-146364022 --shard 9:1-141213431 --shard 10:1-135534747 --shard 11:1-135006516 --shard 12:1-49085594"],
// 		["3", "--shard 12:49085595-133851895 --shard 13:1-115169878 --shard 14:1-107349540 --shard 15:1-102531392 --shard 16:1-90354753 --shard 17:1-81195210 --shard 18:1-78077248 --shard 19:1-59128983 --shard 20:1-63025520 --shard 21:1-48129895 --shard 22:1-51304566 --shard X:1-118966714"],
// 		["4", "--shard X:118966715-155270560 --shard Y:1-59373566 --shard MT:1-16569 --shard GL000207.1:1-4262 --shard GL000226.1:1-15008 --shard GL000229.1:1-19913 --shard GL000231.1:1-27386 --shard GL000210.1:1-27682 --shard GL000239.1:1-33824 --shard GL000235.1:1-34474 --shard GL000201.1:1-36148 --shard GL000247.1:1-36422 --shard GL000245.1:1-36651 --shard GL000197.1:1-37175 --shard GL000203.1:1-37498 --shard GL000246.1:1-38154 --shard GL000249.1:1-38502 --shard GL000196.1:1-38914 --shard GL000248.1:1-39786 --shard GL000244.1:1-39929 --shard GL000238.1:1-39939 --shard GL000202.1:1-40103 --shard GL000234.1:1-40531 --shard GL000232.1:1-40652 --shard GL000206.1:1-41001 --shard GL000240.1:1-41933 --shard GL000236.1:1-41934 --shard GL000241.1:1-42152 --shard GL000243.1:1-43341 --shard GL000242.1:1-43523 --shard GL000230.1:1-43691 --shard GL000237.1:1-45867 --shard GL000233.1:1-45941 --shard GL000204.1:1-81310 --shard GL000198.1:1-90085 --shard GL000208.1:1-92689 --shard GL000191.1:1-106433 --shard GL000227.1:1-128374 --shard GL000228.1:1-129120 --shard GL000214.1:1-137718 --shard GL000221.1:1-155397 --shard GL000209.1:1-159169 --shard GL000218.1:1-161147 --shard GL000220.1:1-161802 --shard GL000213.1:1-164239 --shard GL000211.1:1-166566 --shard GL000199.1:1-169874 --shard GL000217.1:1-172149 --shard GL000216.1:1-172294 --shard GL000215.1:1-172545 --shard GL000205.1:1-174588 --shard GL000219.1:1-179198 --shard GL000224.1:1-179693 --shard GL000223.1:1-180455 --shard GL000195.1:1-182896 --shard GL000212.1:1-186858 --shard GL000222.1:1-186861 --shard GL000200.1:1-187035 --shard GL000193.1:1-189789 --shard GL000194.1:1-191469 --shard GL000225.1:1-211173 --shard GL000192.1:1-547496 --shard NC_007605:1-171823 --shard hs37d5:1-35477943"],
// 		["5", "--shard NO_COOR"]]
Channel
    .fromPath(params.csv)
    .splitCsv(header:false)
    .into { shards1; shards2; shards3; }



bwa_in = bwa_shards.combine(fastq)

process bwa_align {
    cpus 16
    //memory '25 GB'

    input:
	set val(shard), val(type), file(r1), file(r2) from bwa_in

    output:
	set val(type), file("shard_${shard}.bwa.sort.bam"), file("shard_${shard}.bwa.sort.bam.bai") into bwa_shards_ch

    """
    sentieon bwa mem -M -R '@RG\\tID:${name}_${type}\\tSM:${name}_${type}\\tPL:illumina' -K $K_size -t ${task.cpus} -p $genome_file "<sentieon fqidx extract -F $shard/$bwa_num_shards -K $K_size $r1 $r2" | sentieon util sort -r $genome_file -o shard_${shard}.bwa.sort.bam -t ${task.cpus} --sam2bam -i -
    """
}


process bwa_merge_shards {
    input:
	set val(type), file(shard), file(shard_bai) from bwa_shards_ch.groupTuple()

    output:
    file("merged.bam") into merged_bam
    file("merged.bam.bai") into merged_bai
    script:
    bams = shard.join(' ')
    """
    sentieon util merge -o merged.bam ${bams}
    """
}

//(locus_collector_bam, dedup_bam) = merged_bam.into(2)
//(locus_collector_bam, dedup_bam) = merged_bam.into(2)


merged_bam.into{merged_bam1; merged_bam2}
merged_bai.into{merged_bai1; merged_bai2}
//(locus_collector_in, dedup_in) = merged_bam1.combine(shards1).into(2)

process locus_collector {
    cpus 16

    input:
	set val(shard_name), val(shard) from shards1
    each file(bam) from merged_bam1
    each file(bai) from merged_bai1

    output:
    set val(group), file("${shard_name}.score.gz") into locus_collector_scores
   // file("${shard_name}.score.gz.tbi") into locus_collector_tbi
    script:
    group = "scores"
    """
    sentieon driver -t ${task.cpus} -i $bam $shard --algo LocusCollector --fun score_info ${shard_name}.score.gz
    """
    //sentieon driver -t ${task.cpus} -i $bam $shard --algo LocusCollector --fun score_info ${type}.${shard_name}.score.gz
}


locus_collector_scores
   .groupTuple()
   .set{ all_scores }
// //all_scores.printsln()

process dedup {
    cpus 16

    input:
	set val(group), file(score), val(shard_name), val(shard) from all_scores.combine(shards2)
    each file(bam) from merged_bam2
    each file(bai) from merged_bai2

    output:
    //file("dedup.${shard_name}.bam") into merge_deduped_bams
    set val(bam_group), file("${shard_name}.bam") into shard_dedup_bam
    
    //sentieon driver -t 16 -i $bam $shard --algo Dedup --score_info ${score.join(" --score_info " )} --rmdup dedup.${type}.${shard_name}.bam
    //sentieon driver -t 16 -i $bam $shard --algo Dedup --score_info ${all_scores.join(" --score_info ")} --rmdup dedup.${type}.${shard_name}.bam
    script:
    scores = score.join(' --score_info ')
    bam_group = "bams"
    """
    echo " $shard --score_info $scores" > ${shard_name}.bam
    """

}
shard_dedup_bam
    .groupTuple()
    .set{ all_dedup_bams }

shard_combo = Channel.create()
shard_combo = [["1.bam,2.bam"],["1.bam,2.bam,3.bam"],["2.bam,3.bam,4.bam"],["3.bam,4.bam,5.bam"],["4.bam,5.bam"]]



process bqsr {
    cpus 16

    input:
    set val(group), file(bams), val(shard_combo) from all_dedup_bams.combine(shard_combo)
    output:
    file("test")
    script:
    //bam = bams.join('')
    combos = shard_combo.split(',')
    commons = (combos - bams) 
    """
    echo $commons > test
    """
}

