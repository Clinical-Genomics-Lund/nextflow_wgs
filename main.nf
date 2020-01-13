#!/usr/bin/env nextflow

// GENERAL PATHS //
OUTDIR = params.outdir+'/'+params.subdir
CRONDIR = params.crondir

// SENTIEON CONFIGS //
K_size      = 100000000
sentieon_model = params.sentieon_model
bwa_num_shards = params.bwa_shards
bwa_shards = Channel.from( 0..bwa_num_shards-1 )
genomic_num_shards = params.genomic_shards_num

// FASTA //
genome_file = params.genome_file

// VEP REFERENCES AND ANNOTATION DBS //
CADD = params.CADD
VEP_FASTA = params.VEP_FASTA
MAXENTSCAN = params.MAXENTSCAN
VEP_CACHE = params.VEP_CACHE
GNOMAD = params.GNOMAD
PHYLOP =  params.PHYLOP
PHASTCONS = params.PHASTCONS

// ANNOTATION DBS GENERAL //
CLINVAR = params.CLINVAR
KNOWN2 = params.KNOWN2

// RANK MODELS //
rank_model = params.rank_model
rank_model_s = params.rank_model_s

// BED FILES //
inter_bed = params.intersect_bed
scoutbed = params.scoutbed

PON = [F: params.GATK_PON_FEMALE, M: params.GATK_PON_MALE]


csv = file(params.csv)
mode = csv.countLines() > 2 ? "family" : "single"
println(mode)

Channel
    .fromPath(params.csv)
    .splitCsv(header:true)
    .map{ row-> tuple(row.group, row.id, file(row.read1), file(row.read2)) }
    .into { fastq_sharded; fastq; vcf_info }

Channel
    .fromPath(params.csv)
    .splitCsv(header:true)
    .map{ row-> tuple(row.id, row.diagnosis, row.read1, row.read2) }
    .set{ qc_extra }

Channel
    .fromPath(params.csv)
    .splitCsv(header:true)
    .map{ row-> tuple(row.group, row.id, row.sex, row.mother, row.father, row.phenotype, row.diagnosis) }
    .into { ped; all_ids; yml_diag; meta_upd }


Channel
    .fromPath(params.csv)
    .splitCsv(header:true)
    .map{ row-> tuple(row.group, row.id, row.sex, row.type) }
    .set { meta_gatkcov }



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


// Align fractions of fastq files with BWA
process bwa_align_sharded {
	cpus 50
	memory '64 GB'

	input:
		set val(shard), val(group), val(id), r1, r2 from bwa_shards.combine(fastq_sharded)

	output:
		set val(id), file("${id}_${shard}.bwa.sort.bam"), file("${id}_${shard}.bwa.sort.bam.bai") into bwa_shards_ch

	when:
		params.align && params.shardbwa

	"""
	sentieon bwa mem -M \\
		-R '@RG\\tID:${id}\\tSM:${id}\\tPL:illumina' \\
		-K $K_size \\
		-t ${task.cpus} \\
		-p $genome_file '<sentieon fqidx extract -F $shard/$bwa_num_shards -K $K_size $r1 $r2' | sentieon util sort \\
		-r $genome_file \\
		-o ${id}_${shard}.bwa.sort.bam \\
		-t ${task.cpus} --sam2bam -i -
	"""
}

// Merge the fractioned bam files
process bwa_merge_shards {
	cpus 50

	input:
		set val(id), file(shard), file(shard_bai) from bwa_shards_ch.groupTuple()

	output:
		set id, file("${id}_merged.bam"), file("${id}_merged.bam.bai") into merged_bam, qc_merged_bam

	when:
		params.shardbwa
    
	script:
		bams = shard.sort(false) { a, b -> a.getBaseName() <=> b.getBaseName() } .join(' ')

	"""
	sentieon util merge -o ${id}_merged.bam ${bams}
	"""
}

// ALTERNATIVE PATH: Unsharded BWA, utilize local scratch space.
process bwa_align {
	cpus 50
	memory '64 GB'
	scratch true
	stageInMode 'copy'
	stageOutMode 'copy'

	input:
		set val(group), val(id), file(r1), file(r2) from fastq

	output:
		set id, file("${id}_merged.bam"), file("${id}_merged.bam.bai") into bam, qc_bam

	when:
		params.align && !params.shardbwa

	"""
	sentieon bwa mem \\
		-M \\
		-R '@RG\\tID:${id}\\tSM:${id}\\tPL:illumina' \\
		-t ${task.cpus} \\
		$genome_file $r1 $r2 \\
		| sentieon util sort \\
		-r $genome_file \\
		-o ${id}_merged.bam \\
		-t ${task.cpus} --sam2bam -i -
	"""
}



// Collect information that will be used by to remove duplicate reads.
// The output of this step needs to be uncompressed (Sentieon manual uses .gz)
// or the command will occasionally crash in Sentieon 201808.07 (works in earlier)
process locus_collector {
	cpus 16
	errorStrategy 'retry'
	maxErrors 5

	input:
		set id, file(bam), file(bai), val(shard_name), val(shard) from bam.mix(merged_bam).combine(shards1)

	output:
		set val(id), file("${shard_name}_${id}.score"), file("${shard_name}_${id}.score.idx") into locus_collector_scores
		set val(id), file(bam), file(bai) into merged_bam_id

	"""
	sentieon driver \\
		-t ${task.cpus} \\
		-i $bam $shard \\
		--algo LocusCollector \\
		--fun score_info ${shard_name}_${id}.score
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
	sentieon driver \\
		-t ${task.cpus} \\
		-i $bam $shard \\
		--algo Dedup --score_info $scores \\
		--metrics ${shard_name}_${id}_dedup_metrics.txt \\
		--rmdup ${shard_name}_${id}.bam
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
	cpus 54
	memory '64 GB'
	publishDir "${OUTDIR}/qc", mode: 'copy' , overwrite: 'true'

	input:
		set id, file(bam), file(bai), file(dedup) from qc_bam.mix(qc_merged_bam).join(merged_dedup_metrics)

	output:
		set id, file("${id}.QC") into qc_cdm

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
	qc_sentieon.pl $id wgs > ${id}.QC
	"""
}


// Load QC data into CDM (via middleman)
process qc_to_cdm {
        cpus 1
	errorStrategy 'retry'
	maxErrors 5
    	publishDir "${CRONDIR}/qc", mode: 'copy' , overwrite: 'true'
	
	input:
		set id, file(qc) from qc_cdm
		set id, diagnosis, r1, r2 from qc_extra

	output:
		file("${id}.cdm") into cdm_done

	script:
		parts = r1.split('/')
		idx =  parts.findIndexOf {it ==~ /......_......_...._........../}
		rundir = parts[0..idx].join("/")

	"""
	echo "--run-folder $rundir --sample-id $id --subassay $diagnosis --assay wgs --qc ${OUTDIR}/qc/${id}.QC" > ${id}.cdm
	"""
}



process bqsr {
	cpus 16
	errorStrategy 'retry'
	maxErrors 5

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
		--algo QualCal -k $KNOWN2 ${shard_name}_${id}.bqsr.table
	"""
}

// Merge the bqrs shards
process merge_bqsr {
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
	cpus 1
	publishDir "${OUTDIR}/bam", mode: 'copy', overwrite: 'true'

	input:
		set val(id), file(bams), file(bais) from all_dedup_bams4

	output:
		set group, id, file("${id}_merged_dedup.bam"), file("${id}_merged_dedup.bam.bai") into chanjo_bam, expansionhunter_bam, yaml_bam, cov_bam

	script:
		bams_sorted_str = bams.sort(false) { a, b -> a.getBaseName().tokenize("_")[0] as Integer <=> b.getBaseName().tokenize("_")[0] as Integer } .join(' -i ')
		group = "bams"

	"""
	sentieon util merge -i ${bams_sorted_str} -o ${id}_merged_dedup.bam --mergemode 10
	"""
}


// Calculate coverage for chanjo
process chanjo_sambamba {
	cpus 16
	memory '64 GB'
	publishDir "${OUTDIR}/cov"

	input:	
		set group, id, file(bam), file(bai) from chanjo_bam

	output:
		file("${id}.bwa.chanjo.cov") into chanjocov

	"""
	sambamba depth region -t ${task.cpus} -L $scoutbed -T 10 -T 15 -T 20 -T 50 -T 100 $bam > ${id}.bwa.chanjo.cov
	"""
}




////////////////////////////////////////////////////////////////////////
////////////////////////// EXPANSION HUNTER ////////////////////////////
////////////////////////////////////////////////////////////////////////

// call STRs using ExpansionHunter
process expansionhunter {
	input:
		set group, id, file(bam), file(bai) from expansionhunter_bam

	output:
		set group, id, file("${id}.eh.vcf") into expansionhunter_vcf

	"""
	ExpansionHunter \
		--reads $bam \
		--reference $genome_file \
		--variant-catalog $params.expansionhunter_catalog \
		--output-prefix ${id}.eh
	"""
}

// annotate expansionhunter vcf
process stranger {
	input:
		set group, id, file(eh_vcf) from expansionhunter_vcf

	output:
		set group, id, file("${id}.eh.stranger.vcf") into expansionhunter_vcf_anno

	"""
	source activate py3-env
	stranger ${eh_vcf} > ${id}.eh.stranger.vcf
	"""
}

// split multiallelic sites in expansionhunter vcf
// FIXME: Use env variable for picard path...
process vcfbreakmulti_expansionhunter {
	publishDir "${OUTDIR}/vcf", mode: 'copy' , overwrite: 'true'

	input:
		set group, id, file(eh_vcf_anno) from expansionhunter_vcf_anno

	output:
		file "${id}.expansionhunter.vcf.gz" into expansionhunter_scout

	"""
	java -jar /opt/conda/envs/CMD-WGS/share/picard-2.21.2-1/picard.jar RenameSampleInVcf INPUT=${eh_vcf_anno} OUTPUT=${eh_vcf_anno}.rename.vcf NEW_SAMPLE_NAME=${id}
	vcfbreakmulti ${eh_vcf_anno}.rename.vcf > ${id}.expansionhunter.vcf
	bgzip ${id}.expansionhunter.vcf
	tabix ${id}.expansionhunter.vcf.gz
	"""
}
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////






// Do variant calling using DNAscope, sharded
process dnascope {
	cpus 16

	input:
		set id, file(bams), file(bai), file(bqsr), val(shard_name), val(shard), val(one), val(two), val(three) from bam_shard_shard

	output:
		set id, file("${shard_name}_${id}.gvcf"), file("${shard_name}_${id}.gvcf.idx") into vcf_shard

	script:
		combo = [one, two, three] // one two three take on values 0 1 2, 1 2 3...30 31 32
		combo = (combo - 0) //first dummy value removed (0)
		combo = (combo - (genomic_num_shards+1)) //last dummy value removed (32)
		commons = (combo.collect{ "${it}_${id}.bam" })   //add .bam to each combo to match bam files from input channel
		bam_neigh = commons.join(' -i ') 

	"""
	/opt/sentieon-genomics-201711.05/bin/sentieon driver \\
		-t ${task.cpus} \\
		-r $genome_file \\
		-i $bam_neigh $shard \\
		-q $bqsr \\
		--algo DNAscope --emit_mode GVCF ${shard_name}_${id}.gvcf
	"""
}

// Merge gvcf shards
process merge_gvcf {
    cpus 16
	publishDir "${OUTDIR}/gvcf", mode: 'copy' , overwrite: 'true'

    input:
		set id, file(vcfs), file(idx) from vcf_shard.groupTuple()

    output:
		set group, file("${id}.dnascope.gvcf"), file("${id}.dnascope.gvcf.idx") into complete_vcf

    script:
		group = "vcfs"
		vcfs_sorted = vcfs.sort(false) { a, b -> a.getBaseName().tokenize("_")[0] as Integer <=> b.getBaseName().tokenize("_")[0] as Integer } .join(' ')

    """
    /opt/sentieon-genomics-201711.05/bin/sentieon driver \\
        -t ${task.cpus} \\
        --passthru \\
        --algo DNAscope \\
        --merge ${id}.dnascope.gvcf $vcfs_sorted
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
	set group, file("${group}.combined.vcf"), file("${group}.combined.vcf.idx") into combined_vcf

    script:
		all_gvcfs = vcf.join(' -v ')

	"""
	sentieon driver \\
		-t ${task.cpus} \\
		-r $genome_file \\
		--algo GVCFtyper \\
		-v $all_gvcfs ${group}.combined.vcf
	"""
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
    .collectFile(sort: true, storeDir: "${OUTDIR}/ped/")
    .into{ ped_mad; ped_peddy; ped_inher; ped_scout; ped_loqus }


//madeline ped om familj
process madeline {
	publishDir "${OUTDIR}/ped", mode: 'copy' , overwrite: 'true'

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

// Splitting & normalizing variants:
process split_normalize {
	cpus 1
	publishDir "${OUTDIR}/vcf", mode: 'copy', overwrite: 'true'

	input:
		set group, file(vcf), file(idx) from combined_vcf

	output:
		set group, file("${group}.norm.uniq.DPAF.vcf") into split_norm, vcf_gnomad

	"""
	vcfbreakmulti ${vcf} > ${group}.multibreak.vcf
	bcftools norm -m-both -c w -O v -f $genome_file -o ${group}.norm.vcf ${group}.multibreak.vcf
	vcfstreamsort ${group}.norm.vcf | vcfuniq > ${group}.norm.uniq.vcf
	wgs_DPAF_filter.pl ${group}.norm.uniq.vcf > ${group}.norm.uniq.DPAF.vcf
	"""

}

// Intersect VCF, exome/clinvar introns
process intersect {

	input:
		set group, file(vcf) from split_norm

	output:
		set group, file("${group}.intersected.vcf") into split_vep, split_cadd, vcf_loqus

	when:
		params.annotate

	"""
	bedtools intersect -a $vcf -b $inter_bed -u -header > ${group}.intersected.vcf
	"""

}

process add_to_loqusdb {
	cpus 1
	publishDir "${CRONDIR}/loqus", mode: 'copy' , overwrite: 'true'

	input:
		set group, file(vcf) from vcf_loqus
		file(ped) from ped_loqus

	output:
		file("${group}.loqus") into loqusdb_done

	"""
	echo "loqusdb load -f ${ped.toRealPath()} --variant-file ${vcf.toRealPath()}" > ${group}.loqus
	"""
}

process annotate_vep {
	container = '/fs1/resources/containers/ensembl-vep_latest.sif'
	cpus 54

	input:
		set group, file(vcf) from split_vep

	output:
		set group, file("${group}.vep.vcf") into vep

	"""
	vep \\
		-i ${vcf} \\
		-o ${group}.vep.vcf \\
		--offline \\
		--everything \\
		--merged \\
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
		-custom $PHYLOP \\
		-custom $PHASTCONS
	"""
}

// Annotating variants with clinvar

process annotate_clinvar {
        cpus 1
        memory '32GB'

	input:
		set group, file(vcf) from vep

	output:
		set group, file("${group}.clinvar.vcf") into snpsift

	"""
	SnpSift -Xmx60g annotate $CLINVAR \\
		-info CLNSIG,CLNACC,CLNREVSTAT $vcf > ${group}.clinvar.vcf
	"""

}

// Annotating variants with Genmod
process annotate_genmod {
	cpus 2

	input:
		set group, file(vcf) from snpsift

	output:
		set group, file("${group}.genmod.vcf") into genmod

	"""
	genmod annotate --genome-build 38 --spidex $SPIDEX --annotate_regions $vcf -o ${group}.genmod.vcf
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
	cpus 1

	input:
		set group, file(vcf) from inhermod

	output:
		set group, file("${group}.mod.vcf") into mod_vcf

	"""
	modify_vcf_scout.pl $vcf > ${group}.mod.vcf
	"""
} 


// Adding loqusdb allele frequency to info-field: 
// ssh needs to work from anywhere, filesystems mounted on cmdscout
process loqdb {
	cpus 1
	queue 'bigmem'
	errorStrategy 'retry'
	maxErrors 5

	input:
		set group, file(vcf) from mod_vcf

	output:
		set group, file("${group}.loqdb.vcf") into loqdb_vcf

	"""
	export PORT_CMDSCOUT2_MONGODB=33002 #TA BORT VÄLDIGT FULT
	/opt/bin/loqus_db_filter.pl $vcf PORT_CMDSCOUT2_MONGODB > ${group}.loqdb.vcf
	"""
}
// Marking splice INDELs: 
process mark_splice {
	cpus 1

	input:
		set group, file(vcf) from loqdb_vcf

	output:
		set group, file("${group}.marksplice.vcf") into splice_marked

	"""
	/opt/bin/mark_spliceindels.pl $vcf > ${group}.marksplice.vcf
	"""
}

// Extract all INDELs from VCF for CADD annotation
process extract_indels_for_cadd {
	cpus 1

	input:
		set group, file(vcf) from split_cadd
	
	output:
		set group, file("${group}.only_indels.vcf") into indel_cadd_vcf

	"""
	bcftools view $vcf -V snps -o ${group}.only_indels.vcf 
	"""    
}

// Calculate CADD scores for all indels
process calculate_indel_cadd {
	cpus 1
	container = '/fs1/resources/containers/container_cadd_v1.5.sif'
	containerOptions '--bind /tmp/ --bind /local/'

	input:
		set group, file(vcf) from indel_cadd_vcf

	output:
		set group, file("${group}.indel_cadd.gz") into indel_cadd

	"""
        source activate cadd-env
        /opt/cadd/CADD.sh -g GRCh38 -o ${group}.indel_cadd.gz $vcf
	"""
}

// Add the calculated indel CADDs to the vcf
process add_cadd_scores_to_vcf {
	cpus 4

	input: 
		set group, file(vcf) from splice_marked
		set group, file(cadd_scores) from indel_cadd

	output:
		set group, file("${group}.cadd.vcf") into indel_cadd_added

	"""
	gunzip -c $cadd_scores > cadd
	bgzip -@ ${task.cpus} cadd
	tabix -p vcf cadd.gz
	genmod annotate --cadd-file cadd.gz $vcf > ${group}.cadd.vcf
	"""
}

// Scoring variants: 
// Adjusting compound scores: 
// Sorting VCF according to score: 
process genmodscore {
	cpus 2

	input:
		set group, file(vcf) from indel_cadd_added

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
	publishDir "${OUTDIR}/vcf", mode: 'copy', overwrite: 'true'

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
	publishDir "${OUTDIR}/ped", mode: 'copy' , overwrite: 'true'
	cpus 6

	input:
		file(ped) from ped_peddy
		set group, file(vcf), file(idx) from vcf_done1

	output:
		set file("${group}.ped_check.csv"),file("${group}.background_pca.json"),file("${group}.peddy.ped"),file("${group}.html"), file("${group}.het_check.csv"), file("${group}.sex_check.csv"), file("${group}.vs.html") into peddy_files

	"""
	source activate py3-env
	python -m peddy -p ${task.cpus} $vcf $ped --prefix $group
	"""
}




// Extract all variants (from whole genome) with a gnomAD af > x%
process fastgnomad {
	cpus 2
	memory '16 GB'

	publishDir "${OUTDIR}/vcf", mode: 'copy', overwrite: 'true'

    input:
		set group, file(vcf) from vcf_gnomad

	output:
		set group, file("${group}.SNPs.vcf") into vcf_upd, vcf_roh

	"""
	gzip -c $vcf > ${vcf}.gz
	annotate -g $params.FASTGNOMAD_REF -i ${vcf}.gz > ${group}.SNPs.vcf
	"""
	
}


// Call UPD regions from SNP vcf
process upd {
	input:
		set gr, file(vcf) from vcf_upd
		set group, id, sex, mother, father, phenotype, diagnosis from meta_upd

	output:
		file("upd.bed") into upd_plot

	when:
		mode == "family"

	"""
	upd --vcf $vcf --proband $id --mother $mother --father $father --af-tag GNOMADAF regions > upd.bed
	"""
}


// Call ROH regions from SNP vcf
process roh {
	input:
		set gr, file(vcf) from vcf_roh

	output:
		set gr, file("roh.txt") into roh_plot

	"""
    bcftools roh --rec-rate 1e-9 --AF-tag GNOMADAF ${vcf} -o roh.txt
	"""
}

// Create coverage profile using GATK
process gatkcov {
	publishDir "${OUTDIR}/cov", mode: 'copy' , overwrite: 'true'    
    
	cpus 2
	memory '16 GB'

	input:
		set id, group, file(bam), file(bai), gr, sex, type from cov_bam.join(meta_gatkcov, by:1)

	output:
		set group, id, type, file("${id}.standardizedCR.tsv"), file("${id}.denoisedCR.tsv") into cov_plot

	when:
	params.gatkcov

	"""
	source activate gatk4-env

	gatk CollectReadCounts \
		-I $bam -L $params.COV_INTERVAL_LIST \
		--interval-merging-rule OVERLAPPING_ONLY -O ${bam}.hdf5

	gatk --java-options "-Xmx12g" DenoiseReadCounts \
		-I ${bam}.hdf5 --count-panel-of-normals ${PON[sex]} \
		--standardized-copy-ratios ${id}.standardizedCR.tsv \
		--denoised-copy-ratios ${id}.denoisedCR.tsv

	gatk PlotDenoisedCopyRatios \
		--standardized-copy-ratios ${id}.standardizedCR.tsv \
		--denoised-copy-ratios ${id}.denoisedCR.tsv \
		--sequence-dictionary $params.GENOMEDICT \
		--minimum-contig-length 46709983 --output . --output-prefix $id
	"""
}


// Plot ROH, UPD and coverage in a genomic overview plot
process overview_plot {
	publishDir "${OUTDIR}/plots", mode: 'copy' , overwrite: 'true'

	input:
		file(upd) from upd_plot
		set gr, file(roh) from roh_plot
		set group, id, type, file(cov_stand), file(cov_denoised) from cov_plot.groupTuple().view()


	output:
		file("${id[proband_idx]}.genomic_overview.png")

	script:
		proband_idx = type.findIndexOf{ it == "proband" }

	"""
	genome_plotter.pl  --dict $params.GENOMEDICT \\
		 --sample ${id[proband_idx]} \\
		 --upd $upd \\
		 --roh $roh \\
		 --cov ${cov_denoised[proband_idx]} \\
		 --out ${id[proband_idx]}.genomic_overview.png
	"""
}



process create_yaml {
	queue 'bigmem'
	publishDir "${OUTDIR}/yaml", mode: 'copy' , overwrite: 'true'
	publishDir "${CRONDIR}/scout", mode: 'copy' , overwrite: 'true'
	errorStrategy 'retry'
	maxErrors 5

	input:
		set group, id, file(bam), file(bai) from yaml_bam.groupTuple()
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
	export PORT_CMDSCOUT2_MONGODB=33002 #TA BORT VÄLDIGT FULT
	which create_yml.pl
	create_yml.pl \\
		$bams \\
		$group \\
		$OUTDIR \\
		$diagnosis \\
		PORT_CMDSCOUT2_MONGODB \\
		> ${group}.yaml
	"""
}
