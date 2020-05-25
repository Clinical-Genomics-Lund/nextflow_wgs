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

PON = [F: params.GATK_PON_FEMALE, M: params.GATK_PON_MALE]

// Count lines of input csv, if more than 2(header + 1 ind) then mode is set to family //
csv = file(params.csv)
mode = csv.countLines() > 2 ? "family" : "single"
trio = csv.countLines() > 3 ? true : false
println(csv)
println(mode)
println(trio)

workflow.onComplete {

	def msg = """\
		Pipeline execution summary
		---------------------------
		Completed at: ${workflow.complete}
		Duration    : ${workflow.duration}
		Success     : ${workflow.success}
		scriptFile  : ${workflow.scriptFile}
		workDir     : ${workflow.workDir}
		exit status : ${workflow.exitStatus}
		errorMessage: ${workflow.errorMessage}
		errorReport :
		"""
		.stripIndent()
	def error = """\
		${workflow.errorReport}
		"""
		.stripIndent()

	base = csv.getBaseName()
	logFile = file("/fs1/results/cron/logs/" + base + ".complete")
	logFile.text = msg
	logFile.append(error)
}

// Input channels for alignment, variant calling and annotation //
Channel
	.fromPath(params.csv)
	.splitCsv(header:true)
	.map{ row-> tuple(row.group, row.id, file(row.read1), file(row.read2)) }
	.into { input_files; vcf_info }

fastq = Channel.create()
bam_choice = Channel.create()
vcf_choice = Channel.create()
fastq_sharded = Channel.create()
gvcf_choice = Channel.create()
fastq_umi = Channel.create()

// If input files are fastq -> normal path. Flags affecting; --shardbwa (sharded bwa) --align(req), --varcall(if variant calling is to be done) and --annotate(if --varcall)
// bam -> skips align and is variant called (if --varcall is present) and annotated (if --annotate is present)
// vcf skips align + varcall and is only annotated (if --annotate is present)

// If .bam -> value 1, else if .vcf -> value 2 else if .gvcf -> value 4 else if none(.fq.gz) if params.shardbwa true -> value 3 otherwise 0
input_files.view().choice(fastq, bam_choice, vcf_choice, fastq_sharded, gvcf_choice, fastq_umi) { it[2]  =~ /\.bam/ ? 1 : ( it[2] =~ /\.vcf/ ? 2 : ( it[2] =~ /\.gvcf/ ? 4 : (params.shardbwa ? 3 : (params.umi ? 5 : 0))))  }

bam_choice
	.into{ expansionhunter_bam_choice; dnascope_bam_choice; chanjo_bam_choice }

// Input channels for various meta information //
Channel
	.fromPath(params.csv)
	.splitCsv(header:true)
	.map{ row-> tuple(row.id, row.diagnosis, row.read1, row.read2) }
	.set{ qc_extra }

Channel
	.fromPath(params.csv)
	.splitCsv(header:true)
	.map{ row-> tuple(row.group, row.id, row.sex, row.mother, row.father, row.phenotype, row.diagnosis, row.type, row.assay, row.clarity_sample_id, (row.containsKey("ffpe") ? row.ffpe : false), (row.containsKey("analysis") ? row.analysis : false) ) }
	.into { ped; yml_diag; meta_upd; meta_str }


Channel
	.fromPath(params.csv)
	.splitCsv(header:true)
	.map{ row-> tuple(row.group, row.id, row.sex, row.type) }
	.into { meta_gatkcov; meta_exp}


// Check whether genome assembly is indexed //
if(genome_file ){
	bwaId = Channel
			.fromPath("${genome_file}.bwt")
			.ifEmpty { exit 1, "BWA index not found: ${genome_file}.bwt" }
}

// Create genomic shards //
Channel
	.fromPath(params.genomic_shards_file)
	.splitCsv(header:false)
	.into { locuscollector_shards; dedup_shards; genomicshards }

// A channel to pair neighbouring bams and vcfs. 0 and top value removed later
// Needs to be 0..n+1 where n is number of shards in shards.csv
Channel
	.from( 0..(genomic_num_shards+1) )
	.collate( 3,1, false )
	.set{ neighbour_shards }

//merge genomic shards with neighbouring shard combinations
genomicshards
	.merge(tuple(neighbour_shards))
	.into{ bqsr_shard_shard; varcall_shard_shard }

process fastp {
	cpus 10
	tag "$id"
	container = '/fs1/resources/containers/container_twist-brca.sif'
	containerOptions = '--bind /fs1/'

	when:
		params.umi

	input:
		set group, val(id), r1, r2 from fastq_umi

	output:
		set group, val(id), file("${id}_R1_a_q_u_trimmed.fq.gz"),file("${id}_R2_a_q_u_trimmed.fq.gz") into fastq_trimmed

	script:
		"""
		fastp -i $r1 -I $r2 --stdout \\
			-U --umi_loc=per_read --umi_len=3 \\
			-w ${task.cpus} \\
		| fastp --stdin --interleaved_in -f 2 -F 2 \\
			-o ${id}_R1_a_q_u_trimmed.fq.gz \\
			-O ${id}_R2_a_q_u_trimmed.fq.gz \\
			-l 30 \\
			-w ${task.cpus}
		"""
}

// Align fractions of fastq files with BWA
process bwa_align_sharded {
	cpus 50
	memory '64 GB'
	tag "$id $shard"

	input:
		set val(shard), val(group), val(id), r1, r2 from bwa_shards.combine(fastq_sharded)

	output:
		set val(id), group, file("${id}_${shard}.bwa.sort.bam"), file("${id}_${shard}.bwa.sort.bam.bai") into bwa_shards_ch

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
	tag "$id"

	input:
		set val(id), group, file(shard), file(shard_bai) from bwa_shards_ch.groupTuple(by: [0,1])

	output:
		set id, group, file("${id}_merged.bam"), file("${id}_merged.bam.bai") into merged_bam

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
	tag "$id"

	input:
		set val(group), val(id), file(r1), file(r2) from fastq.mix(fastq_trimmed)

	output:
		set id, group, file("${id}_merged.bam"), file("${id}_merged.bam.bai") into bam

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
	tag "$id ($shard_name)"

	input:
		set id, group, file(bam), file(bai), val(shard_name), val(shard) from bam.mix(merged_bam).combine(locuscollector_shards)

	output:
		set val(id), group, file("${shard_name}_${id}.score"), file("${shard_name}_${id}.score.idx") into locus_collector_scores
		set val(id), file(bam), file(bai) into merged_bam_id

	"""
	sentieon driver \\
		-t ${task.cpus} \\
		-i $bam $shard \\
		--algo LocusCollector \\
		--fun score_info ${shard_name}_${id}.score
	"""
}

// Remove duplicate reads
process dedup {
	cpus 16
	cache 'deep'
	tag "$id ($shard_name)"

	input:
		set val(id), group, file(score), file(idx), file(bam), file(bai), val(shard_name), val(shard) \
			from locus_collector_scores.groupTuple(by: [0,1]).join(merged_bam_id).combine(dedup_shards)

	output:
		set val(id), group, file("${shard_name}_${id}.bam"), file("${shard_name}_${id}.bam.bai") into all_dedup_bams_bqsr, all_dedup_bams_dnascope, all_dedup_bams_mergepublish
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

process dedup_metrics_merge {
	tag "$id"

	input:
		set id, file(dedup) from dedup_metrics.groupTuple()

	output:
		set id, file("dedup_metrics.txt") into merged_dedup_metrics

	"""
	sentieon driver --passthru --algo Dedup --merge dedup_metrics.txt $dedup
	"""
}

process bqsr {
	cpus 16
	errorStrategy 'retry'
	maxErrors 5
	tag "$id ($shard_name)"

	input:
		set val(id), group, file(bams), file(bai), val(shard_name), val(shard), val(one), val(two), val(three) from \
			all_dedup_bams_bqsr.groupTuple(by: [0,1]).combine(bqsr_shard_shard)

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
		--algo QualCal -k $params.KNOWN ${shard_name}_${id}.bqsr.table
	"""
}

// Merge the bqrs shards
process merge_bqsr {
	publishDir "${OUTDIR}/bqsr", mode: 'copy', overwrite: 'true'
	tag "$id"

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


process merge_dedup_bam {
	cpus 1
	publishDir "${OUTDIR}/bam", mode: 'copy', overwrite: 'true', pattern: '*.bam*'
	tag "$id"

	input:
		set val(id), group, file(bams), file(bais) from all_dedup_bams_mergepublish.groupTuple(by: [0,1])

	output:
		set group, id, file("${id}_merged_dedup.bam"), file("${id}_merged_dedup.bam.bai") into chanjo_bam, expansionhunter_bam, yaml_bam, cov_bam, bam_manta, bam_nator, bam_tiddit, bam_manta_panel, bam_delly_panel, bam_cnvkit_panel
		set id, group, file("${id}_merged_dedup.bam"), file("${id}_merged_dedup.bam.bai") into qc_bam, bam_melt
		file("${group}.INFO") into bam_INFO

	script:
		bams_sorted_str = bams.sort(false) { a, b -> a.getBaseName().tokenize("_")[0] as Integer <=> b.getBaseName().tokenize("_")[0] as Integer } .join(' -i ')
		bgroup = "bams"

	"""
	sentieon util merge -i ${bams_sorted_str} -o ${id}_merged_dedup.bam --mergemode 10
	echo "BAM	$id	${OUTDIR}/bam/${id}_merged_dedup.bam" > ${group}.INFO
	"""
}

//Collect various QC data: TODO MOVE qc_sentieon to container!
process sentieon_qc {
	cpus 54
	memory '64 GB'
	publishDir "${OUTDIR}/qc", mode: 'copy' , overwrite: 'true'
	tag "$id"
	cache 'deep'

	input:
		set id, group, file(bam), file(bai), file(dedup) from qc_bam.join(merged_dedup_metrics)

	output:
		set id, file("${id}.QC") into qc_cdm, qc_melt

	script:
		target = ""
		panel = ""
		cov = "WgsMetricsAlgo wgs_metrics.txt"
		assay = "wgs"
		if( params.onco ) {
			target = "--interval $params.intervals"
			panel = params.panelhs + "${bam}" + params.panelhs2 
			cov = "CoverageMetrics --cov_thresh 1 --cov_thresh 10 --cov_thresh 30 --cov_thresh 100 --cov_thresh 250 --cov_thresh 500 cov_metrics.txt"
			assay = "panel"
		}
	"""
	sentieon driver \\
		-r $genome_file $target \\
		-t ${task.cpus} \\
		-i ${bam} \\
		--algo MeanQualityByCycle mq_metrics.txt \\
		--algo QualDistribution qd_metrics.txt \\
		--algo GCBias --summary gc_summary.txt gc_metrics.txt \\
		--algo AlignmentStat aln_metrics.txt \\
		--algo InsertSizeMetricAlgo is_metrics.txt \\
		--algo $cov
	$panel
	qc_sentieon.pl $id $assay > ${id}.QC
	"""
}


// Load QC data into CDM (via middleman)
process qc_to_cdm {
	cpus 1
	errorStrategy 'retry'
	maxErrors 5
	publishDir "${CRONDIR}/qc", mode: 'copy' , overwrite: 'true'
	tag "$id"
	
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
	echo "--run-folder $rundir --sample-id $id --subassay $diagnosis --assay $params.assay --qc ${OUTDIR}/qc/${id}.QC" > ${id}.cdm
	"""
}

// Calculate coverage for chanjo
process chanjo_sambamba {
	cpus 16
	memory '64 GB'
	publishDir "${OUTDIR}/cov"
	tag "$id"

	when:
		params.varcall

	input:	
		set group, id, file(bam), file(bai) from chanjo_bam.mix(chanjo_bam_choice)

	output:
		file("${id}.bwa.chanjo.cov") into chanjocov

	"""
	sambamba depth region -t ${task.cpus} -L $params.scoutbed -T 10 -T 15 -T 20 -T 50 -T 100 ${bam.toRealPath()} > ${id}.bwa.chanjo.cov
	"""
}




////////////////////////////////////////////////////////////////////////
////////////////////////// EXPANSION HUNTER ////////////////////////////
////////////////////////////////////////////////////////////////////////

// call STRs using ExpansionHunter
process expansionhunter {
	tag "$id"
	cpus 2

	when:
		params.str
		
	input:
		set group, id, file(bam), file(bai), sex, type \
			from expansionhunter_bam.mix(expansionhunter_bam_choice).join(meta_exp, by: [0,1]).filter { item -> item[5] == 'proband' }

	output:
		set group, id, file("${id}.eh.vcf") into expansionhunter_vcf

	"""
	ExpansionHunter \
		--reads ${bam.toRealPath()} \
		--reference $genome_file \
		--variant-catalog $params.expansionhunter_catalog \
		--output-prefix ${id}.eh
	"""
}

// annotate expansionhunter vcf
process stranger {
	tag "$id"

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
	tag "$id"

	input:
		set group, id, file(eh_vcf_anno) from expansionhunter_vcf_anno
		set group, id, sex, mother, father, phenotype, diagnosis, type, assay, clarity_sample_id, ffpe, analysis from meta_str.filter{ item -> item[7] == 'proband' }

	output:
		file("${id}.expansionhunter.vcf.gz") into expansionhunter_scout
		file("${group}.INFO") into str_INFO

	script:
		if (mode == "family") {
			"""
			java -jar /opt/conda/envs/CMD-WGS/share/picard-2.21.2-1/picard.jar RenameSampleInVcf INPUT=${eh_vcf_anno} OUTPUT=${eh_vcf_anno}.rename.vcf NEW_SAMPLE_NAME=${id}
			vcfbreakmulti ${eh_vcf_anno}.rename.vcf > ${id}.expansionhunter.vcf.tmp
			familyfy_str.pl --vcf ${id}.expansionhunter.vcf.tmp --mother $mother --father $father --out ${id}.expansionhunter.vcf
			bgzip ${id}.expansionhunter.vcf
			tabix ${id}.expansionhunter.vcf.gz
			echo "STR	${OUTDIR}/vcf/${id}.expansionhunter.vcf.gz" > ${group}.INFO
			"""
		}
		else {
			"""
			java -jar /opt/conda/envs/CMD-WGS/share/picard-2.21.2-1/picard.jar RenameSampleInVcf INPUT=${eh_vcf_anno} OUTPUT=${eh_vcf_anno}.rename.vcf NEW_SAMPLE_NAME=${id}
			vcfbreakmulti ${eh_vcf_anno}.rename.vcf > ${id}.expansionhunter.vcf
			bgzip ${id}.expansionhunter.vcf
			tabix ${id}.expansionhunter.vcf.gz
			echo "STR	${OUTDIR}/vcf/${id}.expansionhunter.vcf.gz" > ${group}.INFO
			"""
		}
}
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
process melt_qc_val {
	tag "$id"

	when:
		params.onco

	input:
		set id, qc from qc_melt

	output:
		set id, val(INS_SIZE), val(MEAN_DEPTH), val(COV_DEV) into qc_melt_val, qc_cnvkit_val
	
	script:
		// Collect qc-data if possible from normal sample, if only tumor; tumor
		qc.readLines().each{
			if (it =~ /\"(ins_size_dev)\" : \"(\S+)\"/) {
				ins_dev = it =~ /\"(ins_size_dev)\" : \"(\S+)\"/
			}
			if (it =~ /\"(mean_coverage)\" : \"(\S+)\"/) {
				coverage = it =~ /\"(mean_coverage)\" : \"(\S+)\"/
			}
			if (it =~ /\"(ins_size)\" : \"(\S+)\"/) {
				ins_size = it =~ /\"(ins_size)\" : \"(\S+)\"/
			}
		}
		// might need to be defined for -resume to work "def INS_SIZE" and so on....
		INS_SIZE = ins_size[0][2]
		MEAN_DEPTH = coverage[0][2]
		COV_DEV = ins_dev[0][2]
		"""
		echo hej > hej
		"""
}

// MELT always give VCFs for each type of element defined in mei_list
// If none found -> 0 byte vcf. merge_melt.pl merges the three, if all empty
// it creates a vcf with only header from params.meltheader
// merge_melt.pl gives output ${id}.melt.merged.vcf
process melt {
	cpus 3
	errorStrategy 'retry'
	container = '/fs1/resources/containers/container_twist-brca.sif'
	tag "$id"

	input:
		set id, group, file(bam), file(bai), val(INS_SIZE), val(MEAN_DEPTH), val(COV_DEV) from bam_melt.join(qc_melt_val)

	when:
		params.onco

	output:
		set group, id, file("${id}.melt.merged.vcf") into melt_vcf

	"""
	java -jar  /opt/MELT.jar Single \\
		-bamfile $bam \\
		-r 150 \\
		-h $genome_file \\
		-n $params.bed_melt \\
		-z 50000 \\
		-d 50 -t $params.mei_list \\
		-w . \\
		-b 1/2/3/4/5/6/7/8/9/10/11/12/14/15/16/18/19/20/21/22 \\
		-c $MEAN_DEPTH \\
		-cov $COV_DEV \\
		-e $INS_SIZE
	merge_melt.pl $params.meltheader $id
	"""

}

// When rerunning sample from bam, dnascope has to be run unsharded. this is mixed together with all other vcfs in a trio //
process dnascope_bam_choice {
	cpus 54
	tag "$id"

	when:
		params.varcall

	input:
		set group, id, bam, bqsr from dnascope_bam_choice

	output:
		set group, ph, file("${id}.dnascope.gvcf.gz"), file("${id}.dnascope.gvcf.gz.tbi") into complete_vcf_choice

	script:
	vgroup = "vcfs"
	ph = "dnascope_choice"
	"""
	sentieon driver \\
		-t ${task.cpus} \\
		-r $genome_file \\
		-i ${bam.toRealPath()} \\
		-q $bqsr \\
		--algo DNAscope --emit_mode GVCF ${id}.dnascope.gvcf.gz
	"""
}

// Do variant calling using DNAscope, sharded
process dnascope {
	cpus 16
	tag "$id ($shard_name)"

	when:
		params.varcall

	input:
		set id, group, file(bams), file(bai), file(bqsr), val(shard_name), val(shard), val(one), val(two), val(three) \
			from all_dedup_bams_dnascope.groupTuple(by: [0,1]).join(bqsr_merged.groupTuple()).combine(varcall_shard_shard)
		
	output:
		set id, group, file("${shard_name}_${id}.gvcf.gz"), file("${shard_name}_${id}.gvcf.gz.tbi") into vcf_shard

	script:
		combo = [one, two, three] // one two three take on values 0 1 2, 1 2 3...30 31 32
		combo = (combo - 0) //first dummy value removed (0)
		combo = (combo - (genomic_num_shards+1)) //last dummy value removed (32)
		commons = (combo.collect{ "${it}_${id}.bam" })   //add .bam to each combo to match bam files from input channel
		bam_neigh = commons.join(' -i ') 

	"""
	sentieon driver \\
		-t ${task.cpus} \\
		-r $genome_file \\
		-i $bam_neigh $shard \\
		-q $bqsr \\
		--algo DNAscope --emit_mode GVCF ${shard_name}_${id}.gvcf.gz
	"""
}

// Merge gvcf shards
process merge_gvcf {
	cpus 16
	publishDir "${OUTDIR}/gvcf", mode: 'copy' , overwrite: 'true'
	tag "$id ($group)"

	input:
		set id, group, file(vcfs), file(idx) from vcf_shard.groupTuple(by: [0,1])

	output:
		set group, ph, file("${id}.dnascope.gvcf.gz"), file("${id}.dnascope.gvcf.gz.tbi") into complete_vcf
		set group, id, file("${id}.dnascope.gvcf.gz") into gvcf_gens

	script:
		vgroup = "vcfs"
		vcfs_sorted = vcfs.sort(false) { a, b -> a.getBaseName().tokenize("_")[0] as Integer <=> b.getBaseName().tokenize("_")[0] as Integer } .join(' ')
		ph = "normalpath"
	"""
	sentieon driver \\
		-t ${task.cpus} \\
		--passthru \\
		--algo DNAscope \\
		--merge ${id}.dnascope.gvcf.gz $vcfs_sorted
	"""
}

process gvcf_combine {
	cpus 16
	tag "$group"

	input:
		set vgroup, ph, file(vcf), file(idx) from complete_vcf.mix(complete_vcf_choice).mix(gvcf_choice).groupTuple()
		set val(group), val(id), r1, r2 from vcf_info

	output:
		set group, id, file("${group}.combined.vcf"), file("${group}.combined.vcf.idx") into combined_vcf

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

// Create ped from input variables //
process create_ped {
	tag "$group"

	input:
		set group, id, sex, mother, father, phenotype, diagnosis, type, assay, clarity_sample_id, ffpe, analysis from ped
		

	output:
		file("${group}.ped") into ped_ch
		set id, val(group) into madde_group
		file("${group}.INFO") into tissue_INFO

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
	echo "TISSUE $id $ffpe" > ${group}.INFO
	"""
}

// collects each individual's ped-line and creates one ped-file
ped_ch
	.collectFile(sort: true, storeDir: "${OUTDIR}/ped/")
	.into{ ped_mad; ped_peddy; ped_inher; ped_scout; ped_loqus; ped_prescore; ped_compound }


//madeline ped, run if family mode
process madeline {
	publishDir "${OUTDIR}/ped", mode: 'copy' , overwrite: 'true'

	input:
		file(ped) from ped_mad
		set id, val(group) from madde_group

	output:
		file("${ped}.madeline.xml") into madeline_ped
		file("${group}.INFO") into madde_INFO

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
	echo "MADDE ${ped}.madeline.xml" > ${group}.INFO
	"""
}

// Splitting & normalizing variants:
process split_normalize {
	cpus 1
	publishDir "${OUTDIR}/vcf", mode: 'copy', overwrite: 'true'
	tag "$group"
	cache 'deep'

	when:
		params.annotate

	input:
		set group, id, file(vcf), file(idx) from combined_vcf.mix(vcf_choice)

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
	tag "$group"

	input:
		set group, file(vcf) from split_norm

	output:
		set group, file("${group}.intersected.vcf") into split_vep, split_cadd, vcf_loqus, vcf_cnvkit

	script:

		"""
		bedtools intersect -a $vcf -b $params.intersect_bed -u -header > ${group}.intersected.vcf
		"""

}

process add_to_loqusdb {
	cpus 1
	publishDir "${CRONDIR}/loqus", mode: 'copy' , overwrite: 'true'
	tag "$group"

	input:
		set group, file(vcf) from vcf_loqus
		file(ped) from ped_loqus

	output:
		file("${group}.loqus") into loqusdb_done

	"""
	echo "loqusdb -db $params.loqusdb load -f ${ped.toRealPath()} --variant-file ${vcf.toRealPath()}" > ${group}.loqus
	"""
}

process annotate_vep {
	container = '/fs1/resources/containers/ensembl-vep_latest.sif'
	cpus 54
	tag "$group"

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
		--plugin CADD,$params.CADD \\
		--plugin LoFtool \\
		--plugin MaxEntScan,$params.MAXENTSCAN,SWA,NCSS \\
		--fasta $params.VEP_FASTA \\
		--dir_cache $params.VEP_CACHE \\
		--dir_plugins $params.VEP_CACHE/Plugins \\
		--distance 200 \\
		-cache \\
		-custom $params.GNOMAD_EXOMES,gnomADe,vcf,exact,0,AF_popmax,AF,popmax \\
		-custom $params.GNOMAD_GENOMES,gnomADg,vcf,exact,0,AF_popmax,AF,popmax \\
		-custom $params.PHYLOP \\
		-custom $params.PHASTCONS
	"""
}

// gene, clinvar, loqusdb, enigma(onco)
process vcfanno {
	cpus params.cpu_some
	errorStrategy 'retry'
	cache false
	memory '32GB'
	time '2h'

	input:
		set group, file(vcf) from vep

	output:
		set group, file("${group}.clinvar.loqusdb.gene.vcf") into vcfanno_vcf

	"""
	vcfanno_linux64 -lua /fs1/resources/ref/hg19/bed/scout/sv_tracks/silly.lua $params.vcfanno $vcf > ${group}.clinvar.loqusdb.gene.vcf
	"""
}

// Annotating variants with clinvar
// process annotate_clinvar {
// 	cpus 1
// 	memory '32GB'
// 	tag "$group"

// 	input:
// 		set group, file(vcf) from vep

// 	output:
// 		set group, file("${group}.clinvar.vcf") into snpsift

// 	"""
// 	SnpSift -Xmx60g annotate $params.CLINVAR \\
// 		-info CLNSIG,CLNACC,CLNREVSTAT $vcf > ${group}.clinvar.vcf
// 	"""

// }

// Annotating variants with Genmod
// process annotate_genmod {
// 	cpus 2
// 	tag "$group"

// 	input:
// 		set group, file(vcf) from snpsift

// 	output:
// 		set group, file("${group}.genmod.vcf") into genmod

// 	"""
// 	genmod annotate --genome-build 38 --annotate_regions $vcf -o ${group}.genmod.vcf
// 	"""
// }

// # Annotating variant inheritance models:
process inher_models {
	cpus 6
	memory '64 GB'
	tag "$group"

	input:
		set group, file(vcf) from vcfanno_vcf
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
	tag "$group"

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
// process loqdb {
// 	cpus 1
// 	queue 'bigmem'
// 	errorStrategy 'retry'
// 	maxErrors 5
// 	tag "$group"

// 	input:
// 		set group, file(vcf) from mod_vcf

// 	output:
// 		set group, file("${group}.loqdb.vcf") into loqdb_vcf

// 	"""
// 	export PORT_CMDSCOUT2_MONGODB=33002 #TA BORT VÄLDIGT FULT
// 	/opt/bin/loqus_db_filter.pl $vcf PORT_CMDSCOUT2_MONGODB 38 > ${group}.loqdb.vcf
// 	"""
// }

// Marking splice INDELs: 
process mark_splice {
	cpus 1
	tag "$group"

	input:
		set group, file(vcf) from mod_vcf

	output:
		set group, file("${group}.marksplice.vcf") into splice_marked

	"""
	/opt/bin/mark_spliceindels.pl $vcf > ${group}.marksplice.vcf
	"""
}

// Extract all INDELs from VCF for CADD annotation
process extract_indels_for_cadd {
	cpus 1
	tag "$group"

	input:
		set group, file(vcf) from split_cadd
	
	output:
		set group, file("${group}.only_indels.vcf") into indel_cadd_vep

	"""
	bcftools view $vcf -V snps -o ${group}.only_indels.vcf 
	"""    
}

// Annotate Indels with VEP+Gnomad genomes. Filter variants below threshold
process indel_vep {
	cpus 5
	container = '/fs1/resources/containers/ensembl-vep_latest.sif'
	tag "$group"

	input:
		set group, file(vcf) from indel_cadd_vep

	output:
		set group, file("${group}.only_indels.vep.filtered.vcf") into indel_cadd_vcf
	"""
	vep \\
		-i $vcf \\
		-o ${group}.only_indels.vep.vcf \\
		--offline \\
		--cache \\
		--merged \\
		--vcf \\
		-custom $params.GNOMAD_GENOMES,gnomADg,vcf,exact,0,AF \\
		--dir_cache $params.VEP_CACHE \\
		--force_overwrite \\
		--no_stats \\
		--fork ${task.cpus}
	filter_indels.pl ${group}.only_indels.vep.vcf > ${group}.only_indels.vep.filtered.vcf
	"""
}

// Calculate CADD scores for all indels
process calculate_indel_cadd {
	cpus 5
	container = '/home/cadd_worker/container_cadd_v1.5_hg38_20200117.sif'
	containerOptions '--bind /local/ --bind /home/cadd_worker/'
	scratch '/home/cadd_worker/'
	stageInMode 'copy'
	stageOutMode 'copy'
	queue 'bigmem'
	tag "$group"

	input:
		set group, file(vcf) from indel_cadd_vcf

	output:
		set group, file("${group}.indel_cadd.gz") into indel_cadd

	"""
		export TMPDIR='/home/cadd_worker/'
		source activate cadd-env-v1.5
		/opt/cadd/CADD.sh -g GRCh38 -o ${group}.indel_cadd.gz $vcf
	"""
}

// Add the calculated indel CADDs to the vcf
process add_cadd_scores_to_vcf {
	cpus 4
	tag "$group"

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
	tag "$group"

	input:
		set group, file(vcf) from indel_cadd_added

	output:
		set group, file("${group}.scored.vcf") into scored_vcf

	script:
		if (mode == "family") {
			"""
			genmod score -i $group -c $params.rank_model -r $vcf -o ${group}.score1.vcf
			genmod compound ${group}.score1.vcf > ${group}.score2.vcf
			genmod sort -p -f $group ${group}.score2.vcf -o ${group}.scored.vcf
			"""
		}
		else {
			"""
			genmod score -i $group -c $params.rank_model_s -r $vcf -o ${group}.score1.vcf
			genmod sort -p -f $group ${group}.score1.vcf -o ${group}.scored.vcf
			"""
		}

}

// Bgzipping and indexing VCF: 
process vcf_completion {
	cpus 16
	publishDir "${OUTDIR}/vcf", mode: 'copy', overwrite: 'true', pattern: '*.vcf.gz*'
	tag "$group"

	input:
		set group, file(vcf) from scored_vcf

	output:
		set group, file("${group}.scored.vcf.gz"), file("${group}.scored.vcf.gz.tbi") into vcf_peddy, snv_sv_vcf
		file("${group}.INFO") into snv_INFO

	"""
	bgzip -@ ${task.cpus} $vcf -f
	tabix ${vcf}.gz -f
	echo "SNV	${OUTDIR}/vcf/${group}.scored.vcf.gz" > ${group}.INFO
	"""
}


// Running PEDDY: 
process peddy {
	publishDir "${OUTDIR}/ped", mode: 'copy' , overwrite: 'true'
	container = '/fs1/resources/containers/wgs_20200115.sif'
	cpus 6
	tag "$group"

	input:
		file(ped) from ped_peddy
		set group, file(vcf), file(idx) from vcf_peddy

	output:
		set file("${group}.ped_check.csv"),file("${group}.peddy.ped"), file("${group}.sex_check.csv") into peddy_files
		file("${group}.INFO") into peddy_INFO

	"""
	source activate py3-env
	python -m peddy --sites hg38 -p ${task.cpus} $vcf $ped --prefix $group
	echo "PEDDY	${OUTDIR}/ped/${group}.ped_check.csv,${OUTDIR}/ped/${group}.peddy.ped,${OUTDIR}/ped/${group}.sex_check.csv" > ${group}.INFO
	"""
}

// Extract all variants (from whole genome) with a gnomAD af > x%
process fastgnomad {
	cpus 2
	memory '16 GB'
	tag "$group"
	publishDir "${OUTDIR}/vcf", mode: 'copy', overwrite: 'true'

	when:
		!params.onco

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
	tag "$group"

	input:
		set gr, file(vcf) from vcf_upd
		set group, id, sex, mother, father, phenotype, diagnosis, type, assay, clarity_sample_id, ffpe, analysis from meta_upd.filter{ item -> item[7] == 'proband' }

	output:
		file("upd.bed") into upd_plot
		set group, file("upd.sites.bed") into upd_table

	script:
		if( mode == "family" && trio == true ) {
			"""
			upd --vcf $vcf --proband $id --mother $mother --father $father --af-tag GNOMADAF regions > upd.bed
			upd --vcf $vcf --proband $id --mother $mother --father $father --af-tag GNOMADAF sites > upd.sites.bed
			"""
		}
		else {
			"""
			touch upd.bed
			touch upd.sites.bed
			"""
		}
}


process upd_table {
	publishDir "${OUTDIR}/plots", mode: 'copy' , overwrite: 'true'
	tag "$group"

	input:
		set group, file(upd_sites) from upd_table

	output:
		file("${group}.UPDtable.xls")

	when:
		mode == "family" && trio == true

	"""
	upd_table.pl $upd_sites > ${group}.UPDtable.xls
	"""
}


// Call ROH regions from SNP vcf
process roh {
	tag "$group"

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
	tag "$group"
	cpus 2
	memory '32 GB'

	input:
		set id, group, file(bam), file(bai), gr, sex, type from cov_bam.join(meta_gatkcov, by:1)

	output:
		set group, id, type, sex, file("${id}.standardizedCR.tsv"), file("${id}.denoisedCR.tsv") into cov_plot, cov_gens

	when:
		params.gatkcov

	"""
	source activate gatk4-env

	gatk CollectReadCounts \\
		-I $bam -L $params.COV_INTERVAL_LIST \\
		--interval-merging-rule OVERLAPPING_ONLY -O ${bam}.hdf5

	gatk --java-options "-Xmx30g" DenoiseReadCounts \\
		-I ${bam}.hdf5 --count-panel-of-normals ${PON[sex]} \\
		--standardized-copy-ratios ${id}.standardizedCR.tsv \\
		--denoised-copy-ratios ${id}.denoisedCR.tsv

	gatk PlotDenoisedCopyRatios \\
		--standardized-copy-ratios ${id}.standardizedCR.tsv \\
		--denoised-copy-ratios ${id}.denoisedCR.tsv \\
		--sequence-dictionary $params.GENOMEDICT \\
		--minimum-contig-length 46709983 --output . --output-prefix $id
	"""
}


// Plot ROH, UPD and coverage in a genomic overview plot
process overview_plot {
	publishDir "${OUTDIR}/plots", mode: 'copy' , overwrite: 'true'
	tag "$group"

	input:
		file(upd) from upd_plot
		set gr, file(roh) from roh_plot
		set group, id, type, sex, file(cov_stand), file(cov_denoised) from cov_plot.groupTuple()


	output:
		file("${id[proband_idx]}.genomic_overview.png")

	script:
		proband_idx = type.findIndexOf{ it == "proband" }

	"""
	genome_plotter.pl --dict $params.GENOMEDICT \\
		 --sample ${id[proband_idx]} \\
		 --upd $upd \\
		 --roh $roh \\
		 --sex ${sex[proband_idx]} \\
		 --cov ${cov_denoised[proband_idx]} \\
		 --out ${id[proband_idx]}.genomic_overview.png
	"""
}

process generate_gens_data {
	publishDir "${OUTDIR}/plot_data", mode: 'copy' , overwrite: 'true'
	tag "$group"
	cpus 1

	when:
		!params.onco

	input:
		set id, group, file(gvcf), g, type, sex, file(cov_stand), file(cov_denoise) from gvcf_gens.join(cov_gens, by:[1])

	output:
		set file("${id}.cov.bed.gz"), file("${id}.baf.bed.gz"), file("${id}.cov.bed.gz.tbi"), file("${id}.baf.bed.gz.tbi")

	"""
	generate_gens_data.pl $cov_stand $gvcf $id $params.GENS_GNOMAD
	"""
}

process manta {
	cpus = 56
	publishDir "${OUTDIR}/sv_vcf/", mode: 'copy', overwrite: 'true'
	tag "$id"
	time '24h'
	memory '150GB'

	when:
		params.sv && !params.onco

	input:
		set group, id, file(bam), file(bai) from bam_manta

	output:
		set group, id, file("${id}.manta.vcf.gz") into called_manta

	script:
		bams = bam.join('--bam ')

	"""
	configManta.py --bam $bams --reference $genome_file --runDir . 
	python runWorkflow.py -m local -j ${task.cpus}
	mv results/variants/diploidSV.vcf.gz ${id}.manta.vcf.gz
	mv results/variants/diploidSV.vcf.gz.tbi ${id}.manta.vcf.gz.tbi
	"""
}

process manta_panel {
	cpus = 56
	publishDir "${OUTDIR}/sv_vcf/", mode: 'copy', overwrite: 'true'
	tag "$id"
	time '24h'
	memory '50GB'

	when:
		params.sv && params.onco

	input:
		set group, id, file(bam), file(bai) from bam_manta_panel

	output:
		set group, id, file("${id}.manta.vcf.gz") into called_manta_panel


	"""
	configManta.py --bam $bam --reference $genome_file --runDir . --exome --callRegions $params.bedgz --generateEvidenceBam
	python runWorkflow.py -m local -j ${task.cpus}
	mv results/variants/diploidSV.vcf.gz ${id}.manta.vcf.gz
	mv results/variants/diploidSV.vcf.gz.tbi ${id}.manta.vcf.gz.tbi
	"""
}

process delly_panel {
	cpus = 5
	publishDir "${OUTDIR}/sv_vcf/", mode: 'copy', overwrite: 'true'
	tag "$id"
	time '24h'
	memory '50GB'

	when:
		params.sv && params.onco

	input:
		set group, id, file(bam), file(bai) from bam_delly_panel

	output:
		set group, id, file("${id}.vcf.gz") into called_delly_panel


	"""
	delly call -g $genome_file -o ${id}.bcf $bam 
	bcftools view ${id}.bcf > ${id}.vcf
	bgzip -c ${id}.vcf > ${id}.vcf.gz
	"""
}

process cnvkit_panel {
	cpus = 5
	container = '/fs1/resources/containers/twistmyeloid_active.sif'
	publishDir "${OUTDIR}/sv_vcf/", mode: 'copy', overwrite: 'true'
	tag "$id"
	time '24h'
	memory '20GB'

	when:
		params.sv && params.onco

	input:
		set group, id, file(bam), file(bai) from bam_cnvkit_panel
		set id, val(INS_SIZE), val(MEAN_DEPTH), val(COV_DEV) from qc_cnvkit_val
		set group, file(vcf) from vcf_cnvkit
	
	output:
		set group, id, file("${id}.cnvkit_filtered.vcf") into called_cnvkit_panel

	"""
	cnvkit.py batch $bam -r $params.cnvkit_reference -p 5 -d results/
	cnvkit.py call results/*.cns -v $vcf -o ${id}.call.cns
	filter_cnvkit.pl ${id}.call.cns $MEAN_DEPTH > ${id}.filtered
	cnvkit.py export vcf ${id}.filtered > ${id}.cnvkit_filtered.vcf
	"""

}

process svdb_merge_panel {
	cpus 1
	cache 'deep'
	tag "$group"
	publishDir "${OUTDIR}/sv_vcf/merged/", mode: 'copy', overwrite: 'true'

	input:
		set group, id, file(mantaV) from called_manta_panel.groupTuple()
		set group, id, file(dellyV) from called_delly_panel.groupTuple()
		set group, id, file(melt) from melt_vcf.groupTuple()
		set group, id, file(cnvkitV) from called_cnvkit_panel.groupTuple()
				
	output:
		set group, id, file("${group}.merged.filtered.melt.vcf") into vep_sv_panel, annotsv_panel 
		//set group, id, file("${group}.merged.filtered.vcf") into annotsv_panel

	script:
		tmp = mantaV.collect {it + ':manta ' } + dellyV.collect {it + ':delly ' } + cnvkitV.collect {it + ':cnvkit ' }
		vcfs = tmp.join(' ')

		"""
		source activate py3-env
		svdb --merge --vcf $vcfs --no_intra --pass_only --bnd_distance 2500 --overlap 0.7 --priority manta,delly,cnvkit > ${group}.merged.vcf
		filter_panel_cnv.pl ${group}.merged.vcf $params.intersect_bed > ${group}.merged.filtered.vcf
		vcf-concat ${group}.merged.filtered.vcf $melt | vcf-sort -c > ${group}.merged.filtered.melt.vcf
		"""


}

process tiddit {
	cpus = 2
	publishDir "${OUTDIR}/sv_vcf/", mode: 'copy', overwrite: 'true'   
	time '24h'
	tag "$id"
	memory '32GB'

	when:
		params.sv && !params.onco

	input:
		set group, id, file(bam), file(bai) from bam_tiddit

	output:
		set group, id, file("${id}.tiddit.filtered.vcf") into called_tiddit

	"""
	TIDDIT.py --sv -o ${id}.tiddit --bam $bam
	grep -E \"#|PASS\" ${id}.tiddit.vcf > ${id}.tiddit.filtered.vcf
	"""
}


process cnvnator {
	cpus = 24
	container = '/fs1/resources/containers/wgs_cnvnator_2019-09-06.sif'
	scratch '/local/scratch'
	stageInMode 'copy'
	stageOutMode 'copy'
	time '10h'
	tag "$id"
	memory '80GB'

	when:
		params.sv && !params.onco

	input:
		set group, id, file(bam), file(bai) from bam_nator

	output:
		set group, id, file("${id}.cnvnator_calls*") into cnvnator_subchr

	shell:
	'''
	for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 'X' 'Y'; do
	  cnvnator -root !{id}.root.$chr -chrom $chr -tree !{bam} &&
	  cnvnator -root !{id}.root.$chr -chrom $chr -his 100 &&
	  cnvnator -root !{id}.root.$chr -chrom $chr -stat 100 &&
	  cnvnator -root !{id}.root.$chr -chrom $chr -partition 100 &&
	  cnvnator -root !{id}.root.$chr -chrom $chr -call 100  > !{id}.cnvnator_calls_$chr &
	done
	wait
	'''
}


process merge_cnvnator {
	cpus 20
	container = '/fs1/resources/containers/wgs_cnvnator_2019-09-06.sif'
	publishDir "${OUTDIR}/sv_vcf/", mode: 'copy', overwrite: 'true'
	tag "$group"
	memory '32GB'

	input:
		set group, id, file(chr_vcf) from cnvnator_subchr
	output:
		set group, id, file("${id}.cnvnator.merged.vcf") into merged_cnvnator

	script:
	parts = chr_vcf.join(' ')
	"""
	cat $parts >> ${id}.cnvnator.calls
	cnvnator2VCF.pl -prefix ID -reference GRCh38 ${id}.cnvnator.calls $params.split_ref > ${id}.cnvnator.vcf
	mergeCNVnator.pl ${id}.cnvnator.vcf > ${id}.cnvnator.merged.vcf
	"""
}

process post_cnvnator {
	cpus 1
	
	input:
		set group, id, file(vcf) from merged_cnvnator

	output:
		set group, id, file("${id}.cnvnator.merged.renamed.vcf.gz"), file("${id}.cnvnator.merged.renamed.vcf.gz.tbi") into called_cnvnator
	"""
	java -jar /opt/conda/envs/CMD-WGS/share/picard-2.21.2-1/picard.jar \\
	RenameSampleInVcf INPUT=$vcf OUTPUT=${id}.cnvnator.merged.renamed.vcf NEW_SAMPLE_NAME=$id
	bgzip ${id}.cnvnator.merged.renamed.vcf
	tabix ${id}.cnvnator.merged.renamed.vcf.gz
	"""
}

process svdb_merge {
	cpus 1
	cache 'deep'
	tag "$group"
	publishDir "${OUTDIR}/sv_vcf/merged/", mode: 'copy', overwrite: 'true'

	input:
		set group, id, file(mantaV) from called_manta.groupTuple()
		set group, id, file(tidditV) from called_tiddit.groupTuple()
		set group, id, file(natorV), file(natorI) from called_cnvnator.groupTuple()
		
	output:
		set group, id, file("${group}.merged.vcf") into vcf_vep, annotsv_vcf

	script:

		if (mode == "family") {
			vcfs = []
			manta = []
			tiddit = []
			cnvnator = []
			for (i = 1; i <= mantaV.size(); i++) {
				tmp = mantaV[i-1] + ':manta' + "${i}"
				tmp1 = tidditV[i-1] + ':tiddit' + "${i}"
				tmp2 = natorV[i-1] + ':cnvnator' + "${i}"
				vcfs = vcfs + tmp + tmp1 + tmp2
				mt = 'manta' + "${i}"
				tt = 'tiddit' + "${i}"
				ct = 'cnvnator' + "${i}"
				manta = manta + mt
				tiddit = tiddit + tt
				cnvnator = cnvnator + ct
			}
			prio = manta + tiddit + cnvnator
			prio = prio.join(',')
			vcfs = vcfs.join(' ')
			"""
			source activate py3-env
			svdb --merge --vcf $vcfs --no_intra --pass_only --bnd_distance 2500 --overlap 0.7 --priority $prio > ${group}.merged_tmp.vcf
			merge_callsets.pl ${group}.merged_tmp.vcf > ${group}.merged.vcf
			"""
		}

		else {
			tmp = mantaV.collect {it + ':manta ' } + tidditV.collect {it + ':tiddit ' } + natorV.collect {it + ':cnvnator ' }
			vcfs = tmp.join(' ')
			"""
			source activate py3-env
			svdb --merge --vcf $vcfs --no_intra --pass_only --bnd_distance 2500 --overlap 0.7 --priority manta,tiddit,cnvnator > ${group}.merged.vcf
			"""
		}

}

//create AnnotSV tsv file
process annotsv {
	container = '/fs1/resources/containers/annotsv.v2.3.sif'
	cpus 2
	tag "$group"
	publishDir "${OUTDIR}/annotsv/", mode: 'copy', overwrite: 'true'

	input:
		set group, id, file(sv) from annotsv_vcf.mix(annotsv_panel)
		
	output:
		set group, file("${group}_annotsv.tsv") into annotsv

	"""
	export ANNOTSV="/AnnotSV"
	/AnnotSV/bin/AnnotSV -SvinputFile $sv \\
		-typeOfAnnotation full \\
		-outputDir $group \\
		-genomeBuild GRCh38
	mv $group/*.annotated.tsv ${group}_annotsv.tsv
	"""
}

process vep_sv {
	cpus 56
	container = '/fs1/resources/containers/ensembl-vep_latest.sif'
	tag "$group"
	
	input:
		set group, id, file(vcf) from vcf_vep.mix(vep_sv_panel)

	output:
		set group, id, file("${group}.vep.vcf") into vep_vcf

	"""
	vep \\
		-i $vcf \\
		-o ${group}.vep.vcf \\
		--offline \\
		--merged \\
		--everything \\
		--vcf \\
		--no_stats \\
		--fork ${task.cpus} \\
		--force_overwrite \\
		--plugin LoFtool \\
		--fasta $params.VEP_FASTA \\
		--dir_cache $params.VEP_CACHE \\
		--dir_plugins $params.VEP_CACHE/Plugins \\
		--max_sv_size 50000000 \\
		--distance 200 -cache
	"""
}

process postprocess_vep {
	cpus = 1
	tag "$group"

	input:
		set group, id, file(vcf) from vep_vcf

	output:
		set group, file("${group}.vep.clean.merge.omim.vcf") into artefact_vcf
	
	"""
	python /fs1/viktor/nextflow_svwgs/bin/cleanVCF.py --vcf $vcf > ${group}.vep.clean.vcf
	svdb --merge --overlap 0.9 --notag --vcf ${group}.vep.clean.vcf > ${group}.vep.clean.merge.vcf
	sed -i '3 i ##INFO=<ID=set,Number=1,Type=String,Description="Source VCF for the merged record in SVDB">' ${group}.vep.clean.merge.vcf
    sed -i '3 i ##INFO=<ID=VARID,Number=1,Type=String,Description="The variant ID of merged samples">' ${group}.vep.clean.merge.vcf
	add_omim.pl ${group}.vep.clean.merge.vcf > ${group}.vep.clean.merge.omim.vcf
	"""
}

// Query artefact db
process artefact {
	cpus 1
	tag "$group"

	input:
		set group, file(sv) from artefact_vcf

	output:
		set group, file("${group}.artefact.vcf") into manip_vcf

	"""
	source activate py3-env
	svdb \\
	--sqdb $params.svdb --query \\
	--query_vcf $sv --out_occ ACOUNT --out_frq AFRQ > ${group}.artefact.vcf
	"""
}


process prescore {
	cpus 1
	tag "$group"

	input:
		set group, file(sv_artefact) from manip_vcf
		file(ped) from ped_prescore
		set group, file(annotsv) from annotsv

	output:
		set group, file("${group}.annotatedSV.vcf") into annotatedSV

	"""
	prescore_sv.pl \\
	--sv $sv_artefact --ped $ped --annotsv $annotsv --osv ${group}.annotatedSV.vcf
	"""
}

process score_sv {
	cpus 5
	tag "$group $mode"
	publishDir "${OUTDIR}/vcf", mode: 'copy', overwrite: 'true', pattern: '*.vcf.gz*'

	input:
		set group, file(vcf) from annotatedSV

	output:
		set group, file("${group}.sv.scored.sorted.vcf.gz"), file("${group}.sv.scored.sorted.vcf.gz.tbi") into sv_rescore
		file("${group}.INFO") into sv_INFO
		set group, file("${group}.sv.scored.sorted.vcf.gz") into svvcf_bed
				
	script:
	
		if (mode == "family") {
			"""
			genmod score -i $group -c $params.svrank_model -r $vcf -o ${group}.sv.scored_tmp.vcf
			bcftools sort -O v -o ${group}.sv.scored.sorted.vcf ${group}.sv.scored_tmp.vcf 
			bgzip -@ ${task.cpus} ${group}.sv.scored.sorted.vcf -f
			tabix ${group}.sv.scored.sorted.vcf.gz -f
			echo "SV	${OUTDIR}/vcf/${group}.sv.scored.sorted.vcf.gz" > ${group}.INFO
			"""
		}
		else {
			"""
			genmod score -i $group -c $params.svrank_model_s -r $vcf -o ${group}.sv.scored.vcf
			bcftools sort -O v -o ${group}.sv.scored.sorted.vcf ${group}.sv.scored.vcf
			bgzip -@ ${task.cpus} ${group}.sv.scored.sorted.vcf -f
			tabix ${group}.sv.scored.sorted.vcf.gz -f
			echo "SV	${OUTDIR}/vcf/${group}.sv.scored.sorted.vcf.gz" > ${group}.INFO
			"""
		}
}

process compound_finder {
	cpus 5
	tag "$group $mode"
	publishDir "${OUTDIR}/vcf", mode: 'copy', overwrite: 'true', pattern: '*.vcf.gz*'

	when:
		mode == "family"

	input:
		set group, file(vcf), file(tbi) from sv_rescore
		file(ped) from ped_compound
		set group, file(snv), file(tbi) from snv_sv_vcf

	output:
		set group, file("${group}.snv.rescored.sorted.vcf.gz"), file("${group}.snv.rescored.sorted.vcf.gz.tbi"), \
			file("${group}.sv.rescored.sorted.vcf.gz"), file("${group}.sv.rescored.sorted.vcf.gz.tbi") into vcf_yaml
		file("${group}.INFO") into svcompound_INFO
				
	script:
		"""
		compound_finder.pl \\
			--sv $vcf --ped $ped --snv $snv \\
			--osv ${group}.sv.rescored.sorted.vcf \\
			--osnv ${group}.snv.rescored.sorted.vcf 
		bgzip -@ ${task.cpus} ${group}.sv.rescored.sorted.vcf -f
		bgzip -@ ${task.cpus} ${group}.snv.rescored.sorted.vcf -f
		tabix ${group}.sv.rescored.sorted.vcf.gz -f
		tabix ${group}.snv.rescored.sorted.vcf.gz -f
		echo "SVc	${OUTDIR}/vcf/${group}.sv.rescored.sorted.vcf.gz,${OUTDIR}/vcf/${group}.snv.rescored.sorted.vcf.gz" > ${group}.INFO
		"""

}

// Collects $group.INFO files from each process output that should be included in the yaml for scout loading //
// If a new process needs to be added to yaml. It needs to follow this procedure, as well as be handled in create_yml.pl //
bam_INFO
	.mix(snv_INFO,sv_INFO,str_INFO,peddy_INFO,madde_INFO,svcompound_INFO,tissue_INFO)
	.collectFile()
	.set{ yaml_INFO }
process svvcf_to_bed {
	publishDir "${OUTDIR}/bed", mode: 'copy' , overwrite: 'true'
	tag "group"

	when:
		!params.onco

	input:
		set group, file(vcf) from svvcf_bed

	output:
		file("${group}.sv.bed")

	"""
	svvcf_to_bed.pl $vcf > ${group}.sv.bed
	"""
}

process create_yaml {
	queue 'bigmem'
	publishDir "${OUTDIR}/yaml", mode: 'copy' , overwrite: 'true'
	publishDir "${CRONDIR}/scout", mode: 'copy' , overwrite: 'true'
	errorStrategy 'retry'
	maxErrors 5
	tag "$group"

	input:
		file(INFO) from yaml_INFO
		file(ped) from ped_scout
		set group, id, sex, mother, father, phenotype, diagnosis, type, assay, clarity_sample_id, ffpe, analysis from yml_diag

	output:
		set group, file("${group}.yaml") into yaml

	script:

	"""
	export PORT_CMDSCOUT2_MONGODB=33002 #TA BORT VÄLDIGT FULT
	create_yml.pl \\
		--g $group,$clarity_sample_id --d $diagnosis --p PORT_CMDSCOUT2_MONGODB --out ${group}.yaml --ped $ped --files $INFO --assay $assay,$analysis --antype $params.antype
	"""
}
