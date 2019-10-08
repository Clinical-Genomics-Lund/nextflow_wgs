#!/usr/bin/env nextflow

genome_file = params.genome_file
K_size      = 100000000
KNOWN1 = params.refpath+"/annotation_dbs/1000G_phase1.indels.b37.vcf.gz"
KNOWN2 = params.refpath+"/annotation_dbs/Mills_and_1000G_gold_standard.indels.b37.vcf.gz"
sentieon_model = params.sentieon_model
bwa_num_shards = params.bwa_shards
bwa_shards = Channel.from( 0..bwa_num_shards-1 )
OUTDIR = params.outdir

// temporary, resume does not work for dedup, skipping alignment
name = "4014-14"
vcfdone = file("/fs1/results/vcf/wgs/4014-14.combined.gvcf")
Channel
    .from(name, vcfdone)
    .toList()
    .set{ vcf_temp }

genomic_num_shards = params.genomic_shards_num


CADD = params.CADD
VEP_FASTA = params.VEP_FASTA
MAXENTSCAN = params.MAXENTSCAN
VEP_CACHE = params.VEP_CACHE
GNOMAD = params.GNOMAD
GERP = params.GERP
PHYLOP =  params.PHYLOP
PHASTCONS = params.PHASTCONS

SNPSIFT = params.SNPSIFT
CLINVAR = params.CLINVAR
SWEGEN = params.SWEGEN
SPIDEX = params.SPIDEX
inter_bed = params.intersect_bed
rank_model = params.rank_model
rank_model_s = params.rank_model_s

scoutbed = params.scoutbed



csv = file(params.csv)
mode = csv.countLines() > 2 ? "family" : "single"
println(mode)

Channel
    .fromPath(params.csv)
    .splitCsv(header:true)
    .map{ row-> tuple(row.group, row.id, row.read1, row.read2) }
    .into { fastq; vcf_info }

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




//_${shard}

// Align fractions of fastq files with BWA
process bwa_align {
    cpus 56
    memory '64 GB'

    input:
	set val(shard), val(group), val(id), r1, r2 from bwa_shards.combine(fastq)
    
    output:
	set val(id), file("${id}_${shard}.bwa.sort.bam"), file("${id}_${shard}.bwa.sort.bam.bai") into bwa_shards_ch
    
    when:
    params.align

    """
    sentieon bwa mem -M -R '@RG\\tID:${id}\\tSM:${id}\\tPL:illumina' -K $K_size -t ${task.cpus} -p $genome_file '<sentieon fqidx extract -F $shard/$bwa_num_shards -K $K_size $r1 $r2' | sentieon util sort -r $genome_file -o ${id}_${shard}.bwa.sort.bam -t ${task.cpus} --sam2bam -i -
    """
}

// Merge the fractioned bam files
process bwa_merge_shards {
    cpus 56

    input:
    set val(id), file(shard), file(shard_bai) from bwa_shards_ch.groupTuple()

    output:
    set id, file("${id}_merged.bam"), file("${id}_merged.bam.bai") into merged_bam
        
    script:
    bams = shard.sort(false) { a, b -> a.getBaseName() <=> b.getBaseName() } .join(' ')

    """
    sentieon util merge -o ${id}_merged.bam ${bams}
    """
}
merged_bam.into{merged_bam2; qc_bam;}


// Collect information that will be used by to remove duplicate reads.
// The output of this step needs to be uncompressed (Sentieon manual uses .gz)
// or the command will occasionally crash in Sentieon 201808.07 (works in earlier)
process locus_collector {
    cpus 16

    input:
    set id, file(bam), file(bai), val(shard_name), val(shard) from merged_bam2.combine(shards1)

    output:
    set val(id), file("${shard_name}_${id}.score"), file("${shard_name}_${id}.score.idx") into locus_collector_scores
    set val(id), file(bam), file(bai) into merged_bam_id

    """
    sentieon driver -t ${task.cpus} -i $bam $shard --algo LocusCollector --fun score_info ${shard_name}_${id}.score
    """
}



locus_collector_scores
    .groupTuple()
    .join(merged_bam_id)
    .combine(shards2)
    .set{ all_scores }


// Remove duplicate reads
process dedup {
    cpus 16
    cache 'deep'

    input:
    set val(id), file(score), file(idx), file(bam), file(bai), val(shard_name), val(shard) from all_scores

    output:
    set val(id), file("${shard_name}_${id}.bam"), file("${shard_name}_${id}.bam.bai") into shard_dedup_bam
    set id, file("${shard_name}_${id}_dedup_metrics.txt") into dedup_metrics

    script:
    scores = score.sort(false) { a, b -> a.getBaseName().tokenize("_")[0] as Integer <=> b.getBaseName().tokenize("_")[0] as Integer } .join(' --score_info ')

    """
    sentieon driver -t ${task.cpus} -i $bam $shard --algo Dedup --score_info $scores --metrics ${shard_name}_${id}_dedup_metrics.txt --rmdup ${shard_name}_${id}.bam
    """
    
}

shard_dedup_bam
    .groupTuple()
    .into{ all_dedup_bams1; all_dedup_bams2; all_dedup_bams4 }
//merge shards with shard combinations
shards3
    .merge(tuple(shardie1))
    .into{ shard_shard; shard_shard2 }

process dedup_metrics_merge {

    input:
    set id, file(dedup) from dedup_metrics.groupTuple()

    output:
    set id, file("dedup_metrics.txt") into merged_dedup_metrics

    """
    sentieon driver --passthru --algo Dedup --merge dedup_metrics.txt $dedup
    """
}

//Collect various QC data: TODO MOVE qc_sentieon to container!
process sentieon_qc {
    cpus 56
    memory '64 GB'
    publishDir "${OUTDIR}/postmap/wgs", mode: 'copy' , overwrite: 'true'
    publishDir "${OUTDIR}/cron/qc", mode: 'copy' , overwrite: 'true'

    input:
	set id, file(bam), file(bai), file(dedup) from qc_bam.join(merged_dedup_metrics)

    output:
    file("${id}.QC") into qc_done

    """
    sentieon driver \\
        -r $genome_file -t ${task.cpus} \\
        -i ${bam} \\
        --algo MeanQualityByCycle mq_metrics.txt \\
        --algo QualDistribution qd_metrics.txt \\
        --algo GCBias --summary gc_summary.txt gc_metrics.txt \\
        --algo AlignmentStat aln_metrics.txt \\
        --algo InsertSizeMetricAlgo is_metrics.txt \\
        --algo WgsMetricsAlgo wgs_metrics.txt
    /fs1/pipelines/wgs_germline/annotation/qc_sentieon.pl $id wgs > ${id}.QC
    """

}

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
    commons = combo.collect{ "${it}_${id}.bam" }   //add .bam to each shardie, remove all other bams
    bam_neigh = commons.join(' -i ')

    """
    sentieon driver \\
        -t ${task.cpus} \\
        -r $genome_file \\
        -i $bam_neigh $shard \\
        --algo QualCal -k $KNOWN1 -k $KNOWN2 ${shard_name}_${id}.bqsr.table
    """
}

// Merge the bqrs shards
process merge_bqsr {
    publishDir "${OUTDIR}/bam/wgs/bqsr_tables"

    input:
    set id, file(tables) from bqsr_table.groupTuple()

    output:
    set val(id), file("${id}_merged.bqsr.table") into bqsr_merged

    """
    sentieon driver \\
        --passthru \\
        --algo QualCal \\
        --merge ${id}_merged.bqsr.table $tables
    """
}
bqsr_merged
    .groupTuple()
    .into{ bqsr_merged1; bqsr_merged2;}

all_dedup_bams2
    .join(bqsr_merged1)
    .set{ all_dedup_bams3 }


all_dedup_bams3
    .combine(shard_shard2)
    .set{ bam_shard_shard }


process merge_dedup_bam {
    cpus 16

    publishDir "${OUTDIR}/bam/wgs/", mode: 'copy', overwrite: 'true'

    input:
    set val(id), file(bams), file(bais) from all_dedup_bams4

    output:
    set id, file("${id}_merged_dedup.bam"), file("${id}_merged_dedup.bam.bai") into merged_dedup_bam

    script:
    bams_sorted_str = bams.sort(false) { a, b -> a.getBaseName().tokenize("_")[0] as Integer <=> b.getBaseName().tokenize("_")[0] as Integer } .join(' -i ')
    """
    sentieon util merge -i ${bams_sorted_str} -o ${id}_merged_dedup.bam --mergemode 10
    """
}


process bam_recal {
    cpus 56
    publishDir "${OUTDIR}/bam/wgs/", mode: 'copy', overwrite: 'true'

    input:
    set id, file(bam), file(bai), file(table) from merged_dedup_bam.join(bqsr_merged2)

    output:
    set group, id, file("${id}_recal.bam"), file("${id}_recal.bam.bai"), file("${id}_recal_post*") into merged_recal_dedup_bam

    script:
    group = "bams"

    """
    sentieon driver \\
        -t ${task.cpus} \\
        -r $genome_file \\
        -i $bam \\
        -q $table \\
        --algo QualCal -k $KNOWN1 -k $KNOWN2 ${id}_recal_post \\
        --algo ReadWriter ${id}_recal.bam
    """
}

merged_recal_dedup_bam.into{ mrdb1; mrdb2; mrdb3; }

// process sambamba {
//     cpus 16
//     input:
//     set id, file(bam), file(bai), file(recalval) from mrdb1
//     output:
//     file("${id}_.bwa.chanjo.cov")
//     """
//     sambamba depth region -t ${task.cpus} -L $scoutbed -T 10 -T 15 -T 20 -T 50 -T 100 $bam > ${id}_.bwa.chanjo.cov
//     """
// }


// Do variant calling using DNAscope, sharded
process dnascope {
    cpus 16

    input:
    set id, file(bams), file(bai), file(bqsr), val(shard_name), val(shard), val(one), val(two), val(three) from bam_shard_shard

    output:
    set id, file("${shard_name}_${id}.vcf"), file("${shard_name}_${id}.vcf.idx") into vcf_shard

    script:
    combo = [one, two, three]
    combo = (combo - 0) //first dummy value
    combo = (combo - (genomic_num_shards+1)) //last dummy value
    commons = (combo.collect{ "${it}_${id}.bam" })   //add .bam to each shardie, remove all other bams
    bam_neigh = commons.join(' -i ')
    type = mode == "family" ? "--emit_mode GVCF" : ""

    """
    /opt/sentieon-genomics-201711.05/bin/sentieon driver \\
        -t ${task.cpus} \\
        -r $genome_file \\
        -i $bam_neigh $shard \\
        -q $bqsr \\
        --algo DNAscope $type ${shard_name}_${id}.vcf
    """
}

// Merge vcf shards
process merge_vcf {
    cpus 16

    input:
	set id, file(vcfs), file(idx) from vcf_shard.groupTuple()

    output:
	set group, file("${id}.dnascope.vcf"), file("${id}.dnascope.vcf.idx") into complete_vcf

    script:
    group = "vcfs"
    vcfs_sorted = vcfs.sort(false) { a, b -> a.getBaseName().tokenize("_")[0] as Integer <=> b.getBaseName().tokenize("_")[0] as Integer } .join(' ')
    """
    /opt/sentieon-genomics-201711.05/bin/sentieon driver \\
        -t ${task.cpus} \\
        --passthru \\
        --algo DNAscope \\
        --merge ${id}.dnascope.vcf $vcfs_sorted
    """
}

complete_vcf
    .groupTuple()
    .set{ gvcfs }

process gvcf_combine {
    cpus 16

    input:
    set id, file(vcf), file(idx) from gvcfs
    set val(group), val(id), r1, r2 from vcf_info

    output:
    set group, file("${group}.combined.gvcf"), file("${group}.combined.gvcf.idx") into g_gvcf
    when:

    script:
    // Om fler än en vcf, GVCF combine annars döp om och skickade vidare
    if (mode == "family" ) {
        ggvcfs = vcf.join(' -v ')

        """
        sentieon driver \\
            -t ${task.cpus} \\
            -r $genome_file \\
            --algo GVCFtyper \\
            -v $ggvcfs ${group}.combined.gvcf
        """
    }
    // annars ensam vcf, skicka vidare
    else {
        ggvcf = vcf.join('')
        gidx = idx.join('')

        """
        mv ${ggvcf} ${group}.combined.gvcf
        mv ${gidx} ${group}.combined.gvcf.idx
        """
    }
}

// skapa en pedfil, ändra input istället för sök ersätt?
process create_ped {
    input:
    set group, id, sex, mother, father, phenotype, diagnosis from ped

    output:
    file("${group}.ped") into ped_ch

    script:
    if ( sex =~ /F/) {
        sex = "2"
    }
    else {
        sex = "1"
    }
    if ( phenotype =~ /unaffected/ ) {
        phenotype = "1"
    }
    else {
        phenotype = "2"
    }
    if ( father == "" ) {
        father = "0"
    }
    if ( mother == "" ) {
        mother = "0"
    }

    """
    echo "${group}\t${id}\t${father}\t${mother}\t${sex}\t${phenotype}" > ${group}.ped
    """
}

// collects each individual's ped-line and creates one ped-file
ped_ch
    .collectFile(sort: true, storeDir: "${OUTDIR}/ped/wgs")
    .into{ ped_mad; ped_peddy; ped_inher; ped_scout }


//madeline ped om familj
process madeline {
    publishDir "${OUTDIR}/ped/wgs", mode: 'copy' , overwrite: 'true'

    input:
    file(ped) from ped_mad

    output:
    file("${ped}.madeline.xml") into madeline_ped

    when:
    mode == "family"

    """
    ped_parser \\
        -t ped $ped \\
        --to_madeline \\
        -o ${ped}.madeline
    madeline2 \\
        -L "IndividualId" ${ped}.madeline \\
        -o ${ped}.madeline \\
        -x xml
    """
}

// Intersect VCF, exome/clinvar introns

process intersect {

    input:
    set group, file(vcf), file(idx) from g_gvcf

    output:
    set group, file("${group}.intersected.vcf") into intersected_vcf

    when:
    params.annotate

    """
    bedtools intersect -a $vcf -b $inter_bed -header > ${group}.intersected.vcf
    """

}



// Splitting & normalizing variants:

process split_normalize {
    cpus 16

    input:
    //set group, file(vcf) from vcf_temp
    set group, file(vcf) from intersected_vcf

    output:
    set group, file("${group}.norm.DPAF.vcf") into split
    
    """
    vcfbreakmulti ${vcf} > ${group}.multibreak.vcf
    bcftools norm -m-both -c w -O v -f $genome_file -o ${group}.norm.vcf ${group}.multibreak.vcf
    exome_DPAF_filter.pl ${group}.norm.vcf > ${group}.norm.DPAF.vcf
    """

}


process annotate_vep {
    container = '/fs1/resources/containers/container_VEP.sif'
    cpus 56

    input:
    set group, file(vcf) from split

    output:
    set group, file("${group}.vep.vcf") into vep

    """
    vep \\
        -i ${vcf} \\
        -o ${group}.vep.vcf \\
        --offline \\
        --merged \\
        --everything \\
        --vcf \\
        --no_stats \\
        --fork ${task.cpus} \\
        --force_overwrite \\
        --plugin CADD,$CADD \\
        --plugin LoFtool \\
        --plugin MaxEntScan,$MAXENTSCAN,SWA,NCSS \\
        --fasta $VEP_FASTA \\
        --dir_cache $VEP_CACHE \\
        --dir_plugins $VEP_CACHE/Plugins \\
        --distance 200 \\
        -cache \\
        -custom $GNOMAD \\
        -custom $GERP \\
        -custom $PHYLOP \\
        -custom $PHASTCONS
    """
}

// Annotating variants with SnpSift 2

process snp_sift {
    cpus 16

    input:
    set group, file(vcf) from vep

    output:
    set group, file("${group}.clinvar.vcf") into snpsift

    """
    SnpSift -Xmx60g annotate $CLINVAR \\
        -info CLNSIG,CLNACC,CLNREVSTAT $vcf > ${group}.clinvar.vcf
    """

}

// // Adding SweGen allele frequencies
process swegen_all {
    cpus 16

    input:
    set group, file(vcf) from snpsift

    output:
    set group, file("${group}.swegen.vcf") into sweall

    """
    SnpSift -Xmx60g annotate $SWEGEN \\
        -name swegen \\
        -info AF $vcf > ${group}.swegen.vcf
    """
}
// Annotating variants with Genmod
process annotate_genmod {
    cpus 16

    input:
    set group, file(vcf) from sweall

    output:
    set group, file("${group}.genmod.vcf") into genmod
    
    """
    genmod annotate --spidex $SPIDEX --annotate_regions $vcf -o ${group}.genmod.vcf
    """
}

// # Annotating variant inheritance models:
process inher_models {
    cpus 8
    memory '64 GB'

    input:
    set group, file(vcf) from genmod
    file(ped) from ped_inher

    output:
    set group, file("${group}.models.vcf") into inhermod

    """
    genmod models $vcf -p ${task.cpus} -f $ped > ${group}.models.vcf
    """
}


// Extracting most severe consequence: 
// Modifying annotations by VEP-plugins, and adding to info-field: 
// Modifying CLNSIG field to allow it to be used by genmod score properly:
process modify_vcf {
    cpus 16

    input:
    set group, file(vcf) from inhermod

    output:
    set group, file("${group}.mod.vcf") into mod_vcf

    """
    /opt/bin/modify_vcf_nexomeflow.pl $vcf > ${group}.mod.vcf
    """
} 


// Adding loqusdb allele frequency to info-field: 
// ssh needs to work from anywhere, filesystems mounted on cmdscout
process loqdb {
    cpus 16
    queue 'bigmem'

    input:
    set group, file(vcf) from mod_vcf

    output:
    set group, file("${group}.loqdb.vcf") into loqdb_vcf

    """
    export PORT_CMDSCOUT1_MONGODB=33001 #TA BORT VÄLDIGT FULT
    /opt/bin/loqus_db_filter.pl $vcf PORT_CMDSCOUT1_MONGODB > ${group}.loqdb.vcf
    """
}
// Marking splice INDELs: 
// Annotating delins with cadd: 
process mark_splice {
    cpus 16

    input:
    set group, file(vcf) from loqdb_vcf

    output:
    set group, file("${group}.marksplice.vcf") into splice_marked

    """
    /opt/bin/mark_spliceindels.pl $vcf > ${group}.marksplice.vcf
    """
}

process add_cadd {
    cpus 16
    container = '/fs1/resources/containers/container_cadd_v1.5.sif'
    containerOptions '--bind /tmp/ --bind /local/'

    input:
    set group, file(vcf) from splice_marked

    output:
    set group, file("${group}.marksplice.cadd.vcf") into splice_cadd

    """
    echo test
    /opt/add_missing_CADDs_1.4.sh -i $vcf -o ${group}.marksplice.cadd.vcf -t ./tmp
    """
}

// Scoring variants: 
// Adjusting compound scores: 
// Sorting VCF according to score: 

process genmodscore {
    cpus 16

    input:
    set group, file(vcf) from splice_cadd

    output:
    set group, file("${group}.scored.vcf") into scored_vcf

    script:
    if (mode == "family") {
        """
        genmod score -i $group -c $rank_model -r $vcf -o ${group}.score1.vcf
        genmod compound ${group}.score1.vcf > ${group}.score2.vcf
        genmod sort -p -f $group ${group}.score2.vcf -o ${group}.scored.vcf
        """
    }
    else {
        """
        genmod score -i $group -c $rank_model_s -r $vcf -o ${group}.score1.vcf
        genmod sort -p -f $group ${group}.score1.vcf -o ${group}.scored.vcf
        """
    }

}

// Bgzipping and indexing VCF: 
process vcf_completion {
    cpus 16
    publishDir "${OUTDIR}/vcf/wgs/", mode: 'copy', overwrite: 'true'

    input:
    set group, file(vcf) from scored_vcf

    output:
    set group, file("${group}.scored.vcf.gz"), file("${group}.scored.vcf.gz.tbi") into vcf_done

    """
    bgzip -@ ${task.cpus} $vcf -f
    tabix ${vcf}.gz -f
    """
}

vcf_done.into {
    vcf_done1
    vcf_done2
    vcf_done3
}

// Running PEDDY: 
process peddy {
    publishDir "${OUTDIR}/ped/wgs", mode: 'copy' , overwrite: 'true'
    cpus 6

    input:
    file(ped) from ped_peddy
    set group, file(vcf), file(idx) from vcf_done1

    output:
    set file("${group}.ped_check.csv"),file("${group}.background_pca.json"),file("${group}.peddy.ped"),file("${group}.html"), file("${group}.het_check.csv"), file("${group}.sex_check.csv"), file("${group}.vs.html") into peddy_files
    
    """
    source activate peddy
    python -m peddy -p ${task.cpus} $vcf $ped --prefix $group
    """
}


// Running gSNP:
// process gnsp {
//     //container = 'container_mongodb.sif'
//     publishDir "${OUTDIR}/tmp/wgs/gSNP", mode: 'copy' , overwrite: 'true'
//     input:
//     set group, file(vcf), file(idx) from vcf_done3
//     set group, id, sex, mother, father, phenotype, diagnosis from all_ids.groupTuple()
//     output:
//     set file("${group}_gSNP.tsv"), file("${group}_gSNP.png")
//     when:
//     mode == "family"
//     script:
//     ids = id.join(' ')
//     """
//     perl /opt/bin/gSNP.pl $vcf ${group}_gSNP $ids
//     """
// }

// mrdb2
//     .groupTuple()
//     .set{ yaml_bam }
// Uploading case to scout:
process create_yaml {
    queue 'bigmem'
    publishDir "${OUTDIR}/json/wgs", mode: 'copy' , overwrite: 'true'
    publishDir "${OUTDIR}/cron/scout", mode: 'copy' , overwrite: 'true'

    input:
    set group, id, file(bam), file(bai), file(crap) from mrdb2.groupTuple()
    set group, file(vcf), file(idx) from vcf_done2
    set file(ped_check),file(json),file(peddy_ped),file(html), file(hetcheck_csv), file(sexcheck), file(vs_html) from peddy_files
    file(ped) from ped_scout
    file(xml) from madeline_ped.ifEmpty('single')
    set group, id, sex, mother, father, phenotype, diagnosis from yml_diag

    output:
    set group, file("${group}.yaml") into yaml

    script:
    bams = bam.join(',')
    madde = xml.name != 'single' ? "$xml" : "single"

    """
    export PORT_CMDSCOUT1_MONGODB=33001 #TA BORT VÄLDIGT FULT
    /fs1/pipelines/wgs_germline/annotation/create_yml.pl $bams $ped $group $vcf $madde $peddy_ped $ped_check $sexcheck $OUTDIR $diagnosis PORT_CMDSCOUT1_MONGODB > ${group}.yaml
    """
    ///opt/bin/create_yml.pl
}