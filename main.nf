#!/usr/bin/env nextflow

nextflow.enable.dsl=2


workflow {

	// Print startup and conf output dirs and modes.

	// TODO: Params assignment inside workflow block is a temp solution:
	//       outdir and subdir need to be combined in new var, re-setting
	//       params.outdir won't work.
	params.results_output_dir = params.outdir + '/' + params.subdir
	params.cron_output_dir = params.crondir // TODO: switch back to crondir
	params.mode = file(params.csv).countLines() > 2 ? "family" : "single"
	params.trio = file(params.csv).countLines() > 3 ? true : false

	log.info("Hello.")

	def PON = [F: params.GATK_PON_FEMALE, M: params.GATK_PON_MALE]

	// Count lines of input csv, if more than 2(header + 1 ind) then mode is set to family //

	log.info("Input CSV: " + params.csv)
	log.info("mode: " + params.mode)
	log.info("trio analysis: " + params.trio)
	log.info("Results output dir: " + params.results_output_dir)
	log.info("Results subdir: " + params.subdir)
	log.info("CRON output dir: " + params.cron_output_dir)

	// Print commit-version of active deployment
	// TODO: stuff this one into versions too, for good measure.
	file(params.git)
		.readLines()
		.each { println "git commit-hash: "+it }

	// Print active container
	println("container: " + file(params.container).toRealPath())

	ch_versions = Channel.empty()
	Channel
		.fromPath(params.csv)
		.splitCsv(header: true)
		.set { ch_samplesheet }


	NEXTFLOW_WGS(ch_samplesheet)

	NEXTFLOW_WGS.out.versions.view()
	ch_versions.mix(NEXTFLOW_WGS.out.versions)
	// OUTPUT_VERSIONS(ch_versions)

}

workflow NEXTFLOW_WGS {

	take:
	ch_samplesheet

	main:

	// CHANNEL PREP //
	ch_samplesheet.view()
	ch_fastq = ch_samplesheet
		.filter {
			row -> row.read1.endsWith("fastq.gz") && row.read2.endsWith("fastq.gz")
		}
		.map { row ->
		def group = row.group
		def id = row.id
		def fastq_r1 = row.read1
		def fastq_r2 = row.read2
		tuple(group, id, fastq_r1, fastq_r2) // TODO: filter non fq
	}

	ch_fastq.view()

	// TODO: expand and implement across all processes:
	ch_meta = ch_samplesheet.map{ row->
		tuple(row.group, row.id, row.sex, row.type)
	}

	ch_bam_start = ch_samplesheet
		.filter {
			row -> row.read1.endsWith("bam")
		}
		.map {
			row ->
			def group = row.group
			def id = row.id
			def bam = row.read1
			def bai = row.read2
			tuple(group, id, bam, bai)
		}

	copy_bam(ch_bam_start)
	ch_bam_start = copy_bam.out.bam_bai

	ch_bam_bai = Channel.empty()
	ch_bam_bai = ch_bam_bai.mix(ch_bam_start)
	ch_bam_bai.view()

	// PED //
	ch_ped_input = ch_samplesheet
		.filter { row -> row.type == "proband" }
		.map { row ->
			def group = row.group
			def id = row.id
			def type = row.type
			def sex = row.sex
			def father = row.father
			def mother = row.mother
			tuple(group, id, type, sex, mother, father)
			}

	create_ped(ch_ped_input)
	ch_ped_base = create_ped.out.ped_base
	ch_ped_fa = Channel.empty()
	ch_ped_ma = Channel.empty()

	ch_ped_trio = Channel.empty()
	ch_ped_trio = ch_ped_trio.mix(ch_ped_base)
	if(params.mode == "family" && params.assay == "wgs") {

		ch_ped_fa.mix(create_ped.out.ped_fa)
		ch_ped_ma.mix(create_ped.out.ped_ma)

		ch_ped_trio = ch_ped_base.mix(ch_ped_fa, ch_ped_ma)
		madeline(ch_ped_trio) // TODO: fetch info

	}

	// FASTQ //
	if (params.umi) {
		//TODO: versions!
		fastp(ch_fastq)
		ch_fastq = fastp.out.fastq_trimmed_reads
	}

	// ALIGN //
	//TODO: handle false or remove?
	if (params.align) {
		bwa_align(ch_fastq)
		markdup(bwa_align.out.bam_bai)
		ch_bam_bai = ch_bam_bai.mix(markdup.out.dedup_bam_bai)
	}

	bqsr(ch_bam_bai)


	// SNV CALLING //
	dnascope(ch_bam_bai, bqsr.out.dnascope_bqsr)
	gvcf_combine(dnascope.out.gvcf_tbi.groupTuple())

	ch_split_normalize = gvcf_combine.out.combined_vcf
	ch_split_normalize_concat_vcf = Channel.empty()

	// TODO: move antypes and similar to constants?
	if (params.antype == "panel") {
		freebayes(ch_bam_bai)
		ch_split_normalize_concat_vcf = freebayes.out.freebayes_variants
	}


	// MITO
	if (params.antype == "wgs") { // TODO: if params.mito etc ? will probably mess up split_normalize

		fetch_MTseqs(ch_bam_bai)

		// MITO BAM QC
		sentieon_mitochondrial_qc(fetch_MTseqs.out.bam_bai)
		build_mitochondrial_qc_json(sentieon_mitochondrial_qc.out.qc_tsv)

		// SNVs
		ch_mutect2_input = fetch_MTseqs.out.bam_bai.groupTuple()
		run_mutect2(ch_mutect2_input)
		split_normalize_mito(run_mutect2.out.vcf, ch_meta)
		run_hmtnote(split_normalize_mito.out.vcf)

		ch_split_normalize_concat_vcf = run_hmtnote.out.vcf
		run_haplogrep(run_mutect2.out.vcf)

		// SVs
		run_eklipse(fetch_MTseqs.out.bam_bai, ch_meta)
	}

	// SNV ANNOTATION
	if (params.annotate) {
		// SNPs
		split_normalize(ch_split_normalize, ch_split_normalize_concat_vcf)
		annotate_vep(split_normalize.out.intersected_vcf)
		vcfanno(annotate_vep.out.vcf)
		modify_vcf(vcfanno.out.vcf)
		mark_splice(modify_vcf.out.vcf)

		//INDELS
		extract_indels_for_cadd(split_normalize.out.intersected_vcf)
		indel_vep(extract_indels_for_cadd.out.vcf)
		calculate_indel_cadd(indel_vep.out.vcf)
		bgzip_indel_cadd(calculate_indel_cadd.out.cadd_gz)
		add_cadd_scores_to_vcf(mark_splice.out.splice_marked.join(bgzip_indel_cadd.out.cadd_tbi))

		// INHERITANCE MODELS

		ch_inher_models_input = add_cadd_scores_to_vcf.out.vcf
			.cross(ch_ped_trio)
			.map { vcf_tuple, ped_tuple ->
				def group = vcf_tuple[0]
				def vcf = vcf_tuple[1]
				def type = ped_tuple[1]
				def ped = ped_tuple[2]
				tuple(group, vcf, type, ped) // Combine elements as desired
			}
			.view()

		inher_models(ch_inher_models_input)

		// SCORE VARIANTS //
		genmodscore(inher_models.out.vcf)
		vcf_completion(genmodscore.out.scored_vcf)

		ch_peddy_input_vcf = vcf_completion.out.vcf_tbi
			.filter { vcf ->
				def type = vcf[1] // TODO: how to proof against position change?
				type == "proband"
			}

		ch_peddy_input_vcf.view()
		// TODO: Move this guy to QC:
		peddy(ch_peddy_input_vcf, ch_ped_base)

		// fastgnomad
		fastgnomad(split_normalize.out.norm_uniq_dpaf_vcf)
		// upd
		ch_upd_meta = ch_samplesheet
			.filter { row ->
				row.type == "proband"
			}
			.map { row ->
				tuple(row.group, row.id, row.mother, row.father)
			}

		// upd
		upd(fastgnomad.out.vcf, ch_upd_meta)
		upd_table(upd.out.upd_sites)

		// roh
		roh(fastgnomad.out.vcf)

	}

	if (params.sv) {

		// TODO: define elsewhere
		ch_gatk_ref = Channel
			.fromPath(params.gatkreffolders)
			.splitCsv(header:true)
			.map{ row-> tuple(row.i, row.refpart) }


		gatk_coverage(ch_bam_bai)
		ch_gatk_coverage = gatk_coverage.out.coverage_tsv
		gatk_call_ploidy(ch_gatk_coverage)
		ch_gatk_ploidy = gatk_call_ploidy.out.call_ploidy

		// TODO: do the joining and combining outside
		gatk_call_cnv(ch_gatk_coverage.join(ch_gatk_ploidy, by: [0,1]).combine(ch_gatk_ref))
		//postprocessgatk(gatk_call_cnv.out.gatk_calls)

		// TODO: these two processes can be merged.
		//       antype.panel has an additional arg to manta
		//       and different resource allocation
		ch_manta_out = Channel.empty()
		if (params.antype == "wgs") {
			manta(ch_bam_bai)
			ch_manta_out = ch_manta_out.mix(manta.out.vcf)
		}

		if (params.antype == "panel") {
			manta_panel(ch_bam_bai)
			ch_manta_out = ch_manta_out.mix(manta_panel.out.vcf)
		}


	}


	ch_versions = Channel.empty()
	// ch_versions.mix(fastp.out.versions)
	ch_versions.mix(bwa_align.out.versions)
	ch_versions.mix(markdup.out.versions)
	ch_versions.mix(bqsr.out.versions)
	ch_versions.mix(dnascope.out.versions)
	ch_versions.mix(gvcf_combine.out.versions)
	// TODO: yasnippet this above
	ch_versions.view()

	emit:
		versions = ch_versions
}


	// Input channels for alignment, variant calling and annotation //
	// Channel
	// 	.fromPath(params.csv)
	// 	.splitCsv(header:true)
	// 	.map{ row-> tuple(row.group,
	// 					  row.id,
	// 					  (row.containsKey("bam") ? file(row.bam) : (row.containsKey("vcf") ? file(row.vcf) : file(row.read1) ) ),
	// 					  (row.containsKey("bai") ? file(row.bai) : (row.containsKey("idx") ? file(row.idx) : file(row.read2) ) ) ) }
	// 	.set { ch_samplesheet }

	// fastq = Channel.create()


	// fastq_sharded = Channel.create()
	// fastq_umi = Channel.create()
	// annotate_only = Channel.create()

	// If input-files has bam files bypass alignment, otherwise go for fastq-channels => three options for fastq, sharded bwa, normal bwa or umi trimming
	// TODO: move this out into workflow?
	// input_files.view().choice(bam_choice, fastq, fastq_sharded, fastq_umi, annotate_only ) { it[2] =~ /\.bam/ ? 0 : ( it[2] =~ /\.vcf.gz/ ? 4 : (params.shardbwa ? 2 : (params.umi ? 3 : 1) )) }

	// TODO: annotate-only
	// annotate_only.into{
	// 	annotate_only_vep;
	// 	annotate_only_cadd
	// }

	// Channel
	// 	.fromPath(params.csv)
	// 	.splitCsv(header:true)
	// 	.map{ row-> tuple(row.group, row.assay) }
    //     .set{ meta_loqusdb_no_sv_calling }

	// // Input channels for various meta information //
	// Channel
	// 	.fromPath(params.csv)
	// 	.splitCsv(header:true)
	// 	.map{ row-> tuple(row.id, row.diagnosis, row.read1, row.read2) }
	// 	.set{ qc_extra }

	// Channel
	// 	.fromPath(params.csv)
	// 	.splitCsv(header:true)
	// 	.map{ row-> tuple(row.group, row.id, row.sex, row.mother, row.father, row.phenotype, row.diagnosis, row.type, row.assay, row.clarity_sample_id, (row.containsKey("ffpe") ? row.ffpe : false), (row.containsKey("analysis") ? row.analysis : false) ) }
	// 	.set { ped; yml_diag; meta_upd; meta_str }


	// Channel
	// 	.fromPath(params.csv)
	// 	.splitCsv(header:true)
	// 	.map{ row-> tuple(row.group, row.id, row.sex, row.type) }
	// 	.into { meta_gatkcov; meta_exp; meta_svbed; meta_pod; meta_mutect2; meta_eklipse}


	// Channel
	// 	.fromPath(params.gatkreffolders)
	// 	.splitCsv(header:true)
	// 	.map{ row-> tuple(row.i, row.refpart) }
	// 	.into{ gatk_ref; gatk_postprocess }


	// // Check whether genome assembly is indexed //
	// if(params.genome_file) {
	// 	bwaId = Channel
	// 		.fromPath("${params.genome_file}.bwt")
	// 		.ifEmpty { exit 1, "BWA index not found: ${params.genome_file}.bwt" }
	// }


//}

// workflow.onComplete {

// 	def msg = """\
// 		Pipeline execution summary
// 		---------------------------
// 		Completed at: ${workflow.complete}
// 		Duration    : ${workflow.duration}
// 		Success     : ${workflow.success}
// 		scriptFile  : ${workflow.scriptFile}
// 		workDir     : ${workflow.workDir}
// 		csv         : ${params.csv}
// 		exit status : ${workflow.exitStatus}
// 		errorMessage: ${workflow.errorMessage}
// 		errorReport :
// 		"""
// 		.stripIndent()
// 	def error = """\
// 		${workflow.errorReport}
// 		"""
// 		.stripIndent()

// 	def base = csv.getBaseName()
// 	File logFile = new File("${params.crondir}/logs/${base}.complete")
// 	if (!logFile.getParentFile().exists()) {
// 		logFile.getParentFile().mkdirs()
// 	}
// 	logFile.text = msg
// 	logFile.append(error)
// }


process fastp {
	cpus 10
	tag "$id"
	time '1h'
	memory '20 GB'
	scratch true
	stageInMode 'copy'
	stageOutMode 'copy'
	container "${params.container_fastp}"

	input:
		tuple val(group), val(id), path(fq_r1), path(fq_r2)

	output:
		tuple val(group), val(id), path("${id}_R1_a_q_u_trimmed.fq.gz"), path("${id}_R2_a_q_u_trimmed.fq.gz"), emit: fastq_trimmed_reads
		path("*versions.yml"), emit: versions

	when:
		params.umi

	script:
		"""
		fastp -i $fq_r1 -I $fq_r2 --stdout \\
			-U --umi_loc=per_read --umi_len=3 \\
			-w ${task.cpus} \\
		| fastp --stdin --interleaved_in -f 2 -F 2 \\
			-o ${id}_R1_a_q_u_trimmed.fq.gz \\
			-O ${id}_R2_a_q_u_trimmed.fq.gz \\
			-l 30 \\
			-w ${task.cpus}
		
		${fastp_version(task)}
		"""

	stub:
		"""
		touch "${id}_R1_a_q_u_trimmed.fq.gz"
		touch "${id}_R2_a_q_u_trimmed.fq.gz"

		${fastp_version(task)}
		"""
}
def fastp_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    fastp: \$(echo \$(fastp -v 2>&1) | cut -f 2 -d " ")
	END_VERSIONS
	"""
}


process bwa_align {
	cpus 50
	memory '100 GB' 	// 64 GB peak giab //
	scratch true
	stageInMode 'copy'
	stageOutMode 'copy'
	tag "$id"
	container  "${params.container_sentieon}"

	input:
		tuple val(group), val(id), path(fastq_r1), path(fastq_r2)

	output:
		tuple val(group), val(id), path("${id}_merged.bam"), path("${id}_merged.bam.bai"), emit: bam_bai
		path "*versions.yml", emit: versions

	when:
		params.align

	script:
		"""
		sentieon bwa mem \\
			-M \\
			-K ${params.bwa_K_size} \\
			-R '@RG\\tID:${id}\\tSM:${id}\\tPL:illumina' \\
			-t ${task.cpus} \\
			${params.genome_file} $fastq_r1 $fastq_r2 \\
			| sentieon util sort \\
			-r ${params.genome_file} \\
			-o ${id}_merged.bam \\
			-t ${task.cpus} --sam2bam -i -

		${bwa_align_versions(task)}
		"""

	stub:
		"""
		touch "${id}_merged.bam"
		touch "${id}_merged.bam.bai"

		${bwa_align_versions(task)}
		"""

}
def bwa_align_versions(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    sentieon: \$(echo \$(sentieon util --version 2>&1) | sed -e "s/sentieon-genomics-//g")
	    bwa: \$(echo \$(sentieon bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
	END_VERSIONS
	"""
}

process markdup {
	cpus 40
	errorStrategy 'retry'
	maxErrors 5
	tag "$id"
	memory '50 GB' // 12GB peak GIAB
	time '3h'
	container  "${params.container_sentieon}"
	publishDir "${params.results_output_dir}/bam", mode: 'copy' , overwrite: 'true', pattern: '*_dedup.bam*'

	input:
		tuple val(group), val(id), path(bam), path(bai)

	output:
		tuple val(group), val(id), path("${id}_dedup.bam"), path("${id}_dedup.bam.bai"), emit: dedup_bam_bai
		tuple val(group), val(id), path("dedup_metrics.txt"), emit: dedup_metrics
		tuple val(group), path("${group}_bam.INFO"), emit: dedup_bam_INFO
		path "*versions.yml", emit: versions

	script:
		"""
		sentieon driver \\
			--temp_dir /local/scratch/ \\
			-t ${task.cpus} \\
			-i $bam --shard 1:1-248956422 --shard 2:1-242193529 --shard 3:1-198295559 --shard 4:1-190214555 --shard 5:1-120339935 --shard 5:120339936-181538259 --shard 6:1-170805979 --shard 7:1-159345973 --shard 8:1-145138636 --shard 9:1-138394717 --shard 10:1-133797422 --shard 11:1-135086622 --shard 12:1-56232327 --shard 12:56232328-133275309 --shard 13:1-114364328 --shard 14:1-107043718 --shard 15:1-101991189 --shard 16:1-90338345 --shard 17:1-83257441 --shard 18:1-80373285 --shard 19:1-58617616 --shard 20:1-64444167 --shard 21:1-46709983 --shard 22:1-50818468 --shard X:1-124998478 --shard X:124998479-156040895 --shard Y:1-57227415 --shard M:1-16569 \\
			--algo LocusCollector \\
			--fun score_info ${id}.score

		sentieon driver \\
			--temp_dir /local/scratch/ \\
			-t ${task.cpus} \\
			-i $bam \\
			--algo Dedup --score_info ${id}.score \\
			--metrics dedup_metrics.txt \\
			--rmdup ${id}_dedup.bam

		# TODO: To some separate process?
		echo "BAM	$id	/access/${params.subdir}/bam/${id}_dedup.bam" > ${group}_bam.INFO

		${markdup_versions(task)}
		"""

	stub:
		"""
		touch "${id}_dedup.bam"
		touch "${id}_dedup.bam.bai"
		touch "dedup_metrics.txt"
		touch "${group}_bam.INFO"

 		${markdup_versions(task)}
 		"""
}
def markdup_versions(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
	END_VERSIONS
	"""
}

process copy_bam {

	tag "$id"
	cpus 1
	memory '2GB'
	time '1h'

	input:
		tuple val(group), val(id), path(bam), path(bai)

	output:
		tuple val(group), val(id), path("${id}_dedup.bam"), path("${id}_dedup.bam.bai"), emit: bam_bai
	script:
		"""
		ionice -c 2 -n 7 cp ${bam} "${id}_dedup.copy.bam"
		ionice -c 2 -n 7 cp ${bai} "${id}_dedup.copy.bam.bai"
		"""

	stub:
		"""
		touch "${id}_dedup.copy.bam"
		touch "${id}_dedup.copy.bam.bai"
		"""
}


// // TODO: no
// // For melt to work if started
// process dedupdummy {
// 	input:
// 		tuple val(group), val(id), path(bam), path(bai)
// 	output:
// 		tuple id, path("dummy"), emit: dedup_dummy
// 	when:
// 		params.run_melt

// 	script:
// 	"""
// 	echo test > dummy
// 	"""
// }


process bqsr {
	cpus 40
	errorStrategy 'retry'
	maxErrors 5
	tag "$id"
	memory '30 GB'
	// 12gb peak giab //
	time '5h'
	container  "${params.container_sentieon}"
	publishDir "${params.results_output_dir}/bqsr", mode: 'copy' , overwrite: 'true', pattern: '*.table'

	input:
		tuple val(group), val(id), path(bam), path(bai)

	output:
		tuple val(group), val(id), path("${id}.bqsr.table"), emit: dnascope_bqsr
		path "*versions.yml", emit: versions

	script:
		"""
		sentieon driver -t ${task.cpus} \\
			-r ${params.genome_file} -i $bam \\
			--algo QualCal ${id}.bqsr.table \\
			-k $params.KNOWN_SITES

		${bqsr_version(task)}
		"""

	stub:
		"""
		touch "${id}.bqsr.table"
		${bqsr_version(task)}
		"""
}
def bqsr_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
	END_VERSIONS
	"""
}

// When rerunning sample
process dnascope {
	cpus 54
	memory '100 GB'
	// 12 GB peak giab //
	time '4h'
	tag "$id"
	container  "${params.container_sentieon}"

	input:
		tuple val(group), val(id), val(bam), val(bai)
		tuple val(group2), val(id2), val(bqsr)

	output:
		tuple val(group), val(id), path("${id}.dnascope.gvcf.gz"), path("${id}.dnascope.gvcf.gz.tbi"), emit: gvcf_tbi
		path "*versions.yml", emit: versions

	when:
		params.varcall

	script:
		"""
		sentieon driver \\
			-t ${task.cpus} \\
			-r ${params.genome_file} \\
			-q $bqsr \\
			-i $bam \\
			--shard 1:1-248956422  \\
			--shard 2:1-242193529  \\
			--shard 3:1-198295559  \\
			--shard 4:1-190214555  \\
			--shard 5:1-120339935  \\
			--shard 5:120339936-181538259  \\
			--shard 6:1-170805979  \\
			--shard 7:1-159345973  \\
			--shard 8:1-145138636  \\
			--shard 9:1-138394717  \\
			--shard 10:1-133797422  \\
			--shard 11:1-135086622  \\
			--shard 12:1-56232327  \\
			--shard 12:56232328-133275309  \\
			--shard 13:1-114364328  \\
			--shard 14:1-107043718  \\
			--shard 15:1-101991189  \\
			--shard 16:1-90338345  \\
			--shard 17:1-83257441  \\
			--shard 18:1-80373285  \\
			--shard 19:1-58617616  \\
			--shard 20:1-64444167  \\
			--shard 21:1-46709983  \\
			--shard 22:1-50818468  \\
			--shard X:1-124998478  \\
			--shard X:124998479-156040895  \\
			--shard Y:1-57227415  \\
			--shard M:1-16569 \\
			--algo DNAscope --emit_mode GVCF ${id}.dnascope.gvcf.gz

		${dnascope_version(task)}
		"""

	stub:
		"""
		touch "${id}.dnascope.gvcf.gz"
		touch "${id}.dnascope.gvcf.gz.tbi"

		${dnascope_version(task)}
		"""
}
def dnascope_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
	END_VERSIONS
	"""
}


// //Collect various QC data:
// process sentieon_qc {
// 	cpus 52
// 	memory '30 GB'
// 	tag "$id"
// 	time '2h'
// 	container  "${params.container_sentieon}"

// 	input:
// 		tuple val(id), val(group), path(bam), path(bai)

// 	output:
// 		tuple val(id), val(group), path("mq_metrics.txt"), path("qd_metrics.txt"), path("gc_summary.txt"),
// 		path("gc_metrics.txt"), path("aln_metrics.txt"), path("is_metrics.txt"), path("assay_metrics.txt"),
// 		path("cov_metrics.txt"), path("cov_metrics.txt.sample_summary"), emit: ch_sentieon_qc_metrics
// 		path "*versions.yml", emit: versions

// 	script:
// 		target = ""
// 		// A bit of cheating here - these are really optional arguments
// 		panel_command = "touch cov_metrics.txt cov_metrics.txt.sample_summary"
// 		cov = "WgsMetricsAlgo assay_metrics.txt"

// 		if (params.onco || params.exome) {
// 			target = "--interval $params.intervals"
// 			cov = "CoverageMetrics --cov_thresh 1 --cov_thresh 10 --cov_thresh 30 --cov_thresh 100 --cov_thresh 250 --cov_thresh 500 cov_metrics.txt"
// 			panel_command = "sentieon driver -r ${params.genome_file} -t ${task.cpus} -i ${bam} --algo HsMetricAlgo --targets_list ${params.intervals} --baits_list ${params.intervals} assay_metrics.txt"
// 		}

// 		"""
// 		sentieon driver \\
// 			-r ${params.genome_file} $target \\
// 			-t ${task.cpus} \\
// 			-i $bam \\
// 			--algo MeanQualityByCycle mq_metrics.txt \\
// 			--algo QualDistribution qd_metrics.txt \\
// 			--algo GCBias --summary gc_summary.txt gc_metrics.txt \\
// 			--algo AlignmentStat aln_metrics.txt \\
// 			--algo InsertSizeMetricAlgo is_metrics.txt \\
// 			--algo $cov
// 		$panel_command

// 		${sentieon_qc_version(task)}
// 		"""

// 	stub:
// 		"""
// 		touch "assay_metrics.txt"
// 		touch "mq_metrics.txt"
// 		touch "qd_metrics.txt"
// 		touch "gc_summary.txt"
// 		touch "gc_metrics.txt"
// 		touch "aln_metrics.txt"
// 		touch "is_metrics.txt"
// 		touch "cov_metrics.txt"
// 		touch "cov_metrics.txt.sample_summary"
// 		${sentieon_qc_version(task)}
// 		"""
// }
// def sentieon_qc_version(task) {
// 	"""
// 	cat <<-END_VERSIONS > ${task.process}_versions.yml
// 	${task.process}:
// 	    sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
// 	END_VERSIONS
// 	"""
// }

// process sentieon_qc_postprocess {
// 	cpus 2
// 	memory '1 GB'
// 	tag "$id"
// 	time '2h'

// 	input:
// 		tuple val(id), path(dedup)
// 		tuple val(id), val(group), path(mq_metrics), path(qd_metrics), path(gc_summary), path(gc_metrics), path(aln_metrics)
// 		path (is_metrics), path(assay_metrics), path(cov_metrics), path(cov_metrics_sample_summary)

// 	output:
// 		tuple val(group), val(id), path("${id}_qc.json"), emit: qc_cdm
// 		tuple val(group), val(id), path("${id}_qc.json"), emit: qc_melt

// 	script:

// 		assay = (params.onco || params.exome) ? "panel" : "wgs"
// 		"""
// 		qc_sentieon.pl \\
// 			--SID ${id} \\
// 			--type ${assay} \\
// 			--align_metrics_file ${aln_metrics} \\
// 			--insert_file ${is_metrics} \\
// 			--dedup_metrics_file ${dedup} \\
// 			--metrics_file ${assay_metrics} \\
// 			--gcsummary_file ${gc_summary} \\
// 			--coverage_file ${cov_metrics} \\
// 			--coverage_file_summary ${cov_metrics_sample_summary} \\
// 			> ${id}_qc.json

// 		"""
// }

// process d4_coverage {
// 	cpus 16
// 	memory '10 GB'
// 	publishDir "${params.results_output_dir}/cov", mode: 'copy', overwrite: 'true', pattern: '*.d4'
// 	tag "$id"
// 	container  "${params.container_d4tools}"

// 	input:
// 		tuple val(group), val(id), path(bam), path(bai)

// 	output:
// 		path("${id}_coverage.d4")
// 		tuple val(group), val(id), path("${id}_coverage.d4"), emit: ch_final_d4
// 		path "*versions.yml", emit: versions
// 		tuple val(group), path("${group}_d4.INFO"), emit: d4_INFO

// 	when:
// 		params.run_chanjo2

// 	script:
// 	"""
// 	d4tools create \\
// 		--threads ${task.cpus} \\
// 		"${bam}" \\
// 		"${id}_coverage.d4"

// 	echo "D4	$id	/access/${params.subdir}/cov/${id}_coverage.d4" > ${group}_d4.INFO

// 	${d4_coverage_version(task)}
// 	"""

// 	stub:
// 	"""
// 	touch "${id}_coverage.d4"
// 	touch "${group}_d4.INFO"

// 	${d4_coverage_version(task)}
// 	"""
// }
// def d4_coverage_version(task) {
// 	"""
// 	cat <<-END_VERSIONS > ${task.process}_versions.yml
// 	${task.process}:
// 	    d4tools: \$(echo \$( d4tools 2>&1 | head -1 ) | sed "s/.*version: //" | sed "s/)//" )
// 	END_VERSIONS
// 	"""
// }

// process verifybamid2 {
// 	cpus 16
// 	memory '10 GB'
// 	// publishDir "${params.results_output_dir}/contamination", mode: 'copy', overwrite: 'true', pattern: '*.selfSM'
// 	tag "$id"
// 	container  "${params.container_verifybamid2}"

// 	input:
// 		tuple val(group), val(id), path(bam), path(bai)

// 	output:
// 		path("${id}.result.selfSM")
// 		path("${id}.result.Ancestry")
// 		path "*versions.yml", emit: versions

// 	script:
// 		if ( params.antype == "wgs") {
// 			"""
// 			verifybamid2 \
// 				--SVDPrefix ${params.verifybamid2_svdprefix} \
// 				--Reference ${params.genome_file} \
// 				--BamFile ${bam}

// 				mv result.selfSM ${id}.result.selfSM
// 				mv result.Ancestry ${id}.result.Ancestry
// 			${verifybamid2_version(task)}
// 			"""
// 		}
// 		else {
// 			"""
// 			verifybamid2 \
// 				--DisableSanityCheck \
// 				--SVDPrefix ${params.verifybamid2_svdprefix} \
// 				--Reference ${params.genome_file} \
// 				--BamFile ${bam}

// 				mv result.selfSM ${id}.result.selfSM
// 				mv result.Ancestry ${id}.result.Ancestry
// 			${verifybamid2_version(task)}
// 			"""
// 		}


// 	stub:
// 		"""
// 		touch "${id}.result.selfSM"
// 		touch "${id}.result.Ancestry"

// 		${verifybamid2_version(task)}
// 		"""
// }
// def verifybamid2_version(task) {
// 	"""
// 	cat <<-END_VERSIONS > ${task.process}_versions.yml
// 	${task.process}:
// 	    VerifyBamID2: \$( echo \$( verifybamid2 --help 2>&1 | grep Version ) | sed "s/^.*Version://" )
// 	END_VERSIONS
// 	"""
// }

// // Calculate coverage for paneldepth
// process depth_onco {
// 	cpus 2
// 	time '1h'
// 	memory '10 GB'
// 	publishDir "${params.results_output_dir}/cov", mode: 'copy', overwrite: 'true'
// 	tag "$id"
// 	input:
// 		tuple val(group), val(id), path(bam), path(bai)

// 	output:
// 		path("${id}.lowcov.overlapping.bed"), emit: cov_onco


// 	when:
// 		params.assay == "swea"

// 	script:
// 		"""
// 		panel_depth.pl $bam $params.scoutbed > ${id}.lowcov.bed
// 		overlapping_genes.pl ${id}.lowcov.bed $params.gene_regions > ${id}.lowcov.overlapping.bed
// 		"""

// 	stub:
// 		"""
// 		touch "${id}.lowcov.overlapping.bed"
// 		"""
// }

// process SMNCopyNumberCaller {
// 	cpus 10
// 	memory '25GB'
// 	time '2h'
// 	publishDir "${params.results_output_dir}/plots/SMNcnc", mode: 'copy' , overwrite: 'true', pattern: '*.pdf*'
// 	tag "$id"

// 	input:
// 		tuple val(group), val(id), path(bam), path(bai)

// 	output:
// 		path("*.tsv"), emit: smn_tsv
// 		tuple path("*.pdf"), path("*.json")
// 		tuple val(group), path("${group}_smn.INFO"), emit: smn_INFO
// 		path "*versions.yml", emit: versions

// 	when:
// 		params.antype == "wgs"

// 	script:
// 		"""
// 		samtools view -H $bam | \\
// 			sed -e 's/SN:1/SN:chr1/' | sed -e 's/SN:2/SN:chr2/' |  \\
// 			sed -e 's/SN:3/SN:chr3/' | sed -e 's/SN:4/SN:chr4/' |  \\
// 			sed -e 's/SN:5/SN:chr5/' | sed -e 's/SN:6/SN:chr6/' |  \\
// 			sed -e 's/SN:7/SN:chr7/' | sed -e 's/SN:8/SN:chr8/' |  \\
// 			sed -e 's/SN:9/SN:chr9/' | sed -e 's/SN:10/SN:chr10/' | \\
// 			sed -e 's/SN:11/SN:chr11/' | sed -e 's/SN:12/SN:chr12/' |  \\
// 			sed -e 's/SN:13/SN:chr13/' | sed -e 's/SN:14/SN:chr14/' |  \\
// 			sed -e 's/SN:15/SN:chr15/' | sed -e 's/SN:16/SN:chr16/' |  \\
// 			sed -e 's/SN:17/SN:chr17/' | sed -e 's/SN:18/SN:chr18/' |  \\
// 			sed -e 's/SN:19/SN:chr19/' | sed -e 's/SN:20/SN:chr20/' |  \\
// 			sed -e 's/SN:21/SN:chr21/' | sed -e 's/SN:22/SN:chr22/' |  \\
// 			sed -e 's/SN:X/SN:chrX/' | sed -e 's/SN:Y/SN:chrY/' |   \\
// 			sed -e 's/SN:MT/SN:chrM/' | \\
// 			samtools reheader - $bam > ${id}.bam
// 		samtools index -b ${id}.bam -@ ${task.cpus}
// 		echo ${id}.bam > manifest.txt
// 		smn_caller.py --manifest manifest.txt --genome 38 --prefix ${id} --outDir . --threads ${task.cpus}
// 		rm ${id}.bam
// 		source activate py3-env
// 		python /SMNCopyNumberCaller/smn_charts.py -s ${id}.json -o .
// 		mv ${id}.tsv ${group}_SMN.tsv
// 		echo "SMN ${params.accessdir}/smn/${group}_SMN.tsv" > ${group}_smn.INFO

// 		${smn_copy_number_caller_version(task)}
// 		"""

// 	stub:
// 		"""
// 		touch "${id}.bam"
// 		touch "${id}.tsv"
// 		touch "${id}.pdf"
// 		touch "${id}.json"
// 		touch "${group}_SMN.tsv"
// 		touch "${group}_smn.INFO"

// 		${smn_copy_number_caller_version(task)}
// 		"""
// }
// // collects each individual's SMNCNC-tsv and creates one tsv-file
// smn_tsv
// 	.collectPath(keepHeader: true, storeDir: "${params.results_output_dir}/smn/")
// def smn_copy_number_caller_version(task) {
// 	"""
// 	cat <<-END_VERSIONS > ${task.process}_versions.yml
// 	${task.process}:
// 	    samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
// 	    smn-copy-number-caller: 1.1.2
// 	END_VERSIONS
// 	"""
// }


// ////////////////////////////////////////////////////////////////////////
// ////////////////////////// EXPANSION HUNTER ////////////////////////////
// ////////////////////////////////////////////////////////////////////////

// // call STRs using ExpansionHunter, and plot alignments with GraphAlignmentViewer
// process expansionhunter {
// 	tag "$group"
// 	cpus 2
// 	time '10h'
// 	memory '40 GB'

// 	input:
// 		tuple val(group), val(id), path(bam), path(bai), sex, type \


// 	output:
// 		tuple val(group), val(id), path("${group}.eh.vcf"), emit: expansionhunter_vcf
// 		tuple val(group), val(id), path("${group}.eh_realigned.sort.bam"), path("${group}.eh_realigned.sort.bam.bai"), path("${group}.eh.vcf"), emit: reviewer
// 		path "*versions.yml", emit: versions

// 	when:
// 		params.str

// 	script:
// 		"""
// 		source activate htslib10
// 		ExpansionHunter \
// 			--reads $bam \
// 			--reference ${params.genome_file} \
// 			--variant-catalog $params.expansionhunter_catalog \
// 			--output-prefix ${group}.eh
// 		samtools sort ${group}.eh_realigned.bam -o ${group}.eh_realigned.sort.bam
// 		samtools index ${group}.eh_realigned.sort.bam

// 		${expansionhunter_version(task)}
// 		"""

// 	stub:
// 		"""
// 		source activate htslib10
// 		touch "${group}.eh.vcf"
// 		touch "${group}.eh_realigned.sort.bam"
// 		touch "${group}.eh_realigned.sort.bam.bai"

// 		${expansionhunter_version(task)}
// 		"""
// }
// def expansionhunter_version(task) {
// 	"""
// 	cat <<-END_VERSIONS > ${task.process}_versions.yml
// 	${task.process}:
// 	    expansionhunter: \$(echo \$(ExpansionHunter --version 2>&1) | sed 's/.*ExpansionHunter v// ; s/]//')
// 	    samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
// 	END_VERSIONS
// 	"""
// }

// // annotate expansionhunter vcf
// process stranger {
// 	tag "$group"
// 	memory '1 GB'
// 	time '10m'
// 	cpus 2
// 	container  "${params.container_stranger}"

// 	input:
// 		tuple val(group), val(id), path(eh_vcf)

// 	output:
// 		tuple val(group), val(id), path("${group}.fixinfo.eh.stranger.vcf"), emit: expansionhunter_vcf_anno
// 		path "*versions.yml", emit: versions

// 	script:
// 		"""
// 		stranger ${eh_vcf} -f $params.expansionhunter_catalog > ${group}.eh.stranger.vcf
// 		grep ^# ${group}.eh.stranger.vcf > ${group}.fixinfo.eh.stranger.vcf
// 		grep -v ^# ${group}.eh.stranger.vcf | sed 's/ /_/g' >> ${group}.fixinfo.eh.stranger.vcf

// 		${stranger_version(task)}
// 		"""

// 	stub:
// 		"""
// 		touch "${group}.fixinfo.eh.stranger.vcf"
// 		${stranger_version(task)}
// 		"""
// }
// def stranger_version(task) {
// 	"""
// 	cat <<-END_VERSIONS > ${task.process}_versions.yml
// 	${task.process}:
// 	    stranger: \$( stranger --version )
// 	END_VERSIONS
// 	"""
// }

// //for i in $( ls *.svg | cut -f 2 -d "." ); do echo "STR_IMG $i /access/!{params.subdir}/plots/reviewer/!{group}/!{group}.${i}.svg" >> !{group}_rev.INFO; done
// process reviewer {
// 	tag "$group"
// 	cpus 2
// 	time '1h'
// 	memory '1 GB'
// 	errorStrategy 'ignore'
// 	container  "${params.container_reviewer}"
// 	publishDir "${params.results_output_dir}/plots/reviewer/${group}", mode: 'copy' , overwrite: 'true', pattern: '*.svg'

// 	input:
// 		tuple val(group), val(id), path(bam), path(bai), path(vcf)

// 	output:
// 		path("*svg")
// 		//tuple val(group), path("${group}_rev.INFO"), emit: reviewer_INFO
// 		path "*versions.yml", emit: versions

// 	shell:
// 		version_str = reviewer_version(task)
// 		'''
// 		grep LocusId !{params.expansionhunter_catalog} | sed 's/[",^ ]//g' | cut -d':' -f2 | perl -na -e 'chomp; \
// 		system("REViewer --reads !{bam} \
// 		--vcf !{vcf} \
// 		--reference !{genome_file} \
// 		--catalog !{params.expansionhunter_catalog} \
// 		--locus $_ \
// 		--output-prefix !{id}");'

// 		echo "!{version_str}" > "!{task.process}_versions.yml"
// 		'''

// 	stub:
// 		version_str = reviewer_version(task)
// 		"""
// 		touch "${id}.svg"
// 		echo "${version_str}" > "${task.process}_versions.yml"
// 		"""
// }
// def reviewer_version(task) {
// 	// This docstring looks different
// 	// If spaces similarly to the others, this leads to additional whitespace above and below the version text
// 	"""${task.process}:
// 	    reviewer: \$(echo \$(REViewer --version 2>&1) | sed 's/^.*REViewer v//')"""
// }

// // split multiallelic sites in expansionhunter vcf
// // FIXME: Use env variable for picard path...
// process vcfbreakmulti_expansionhunter {
// 	publishDir "${params.results_output_dir}/vcf", mode: 'copy' , overwrite: 'true', pattern: '*.vcf.gz'
// 	tag "$group"
// 	time '1h'
// 	memory '50 GB'

// 	input:
// 		tuple val(group), val(id), path(eh_vcf_anno), val(sex), val(mother), val(father), val(phenotype), val(diagnosis), val(type), val(assay), val(clarity_sample_id), val(ffpe), val(analysis)

// 	output:
// 		path("${group}.expansionhunter.vcf.gz"), emit: expansionhunter_scout
// 		tuple val(group), path("${group}_str.INFO"), emit: str_INFO
// 		path "*versions.yml", emit: versions

// 	script:
// 		if (father == "") { father = "null" }
// 		if (mother == "") { mother = "null" }
// 		if (params.mode == "family") {
// 			"""
// 			java -jar /opt/conda/envs/CMD-WGS/share/picard-2.21.2-1/picard.jar RenameSampleInVcf INPUT=${eh_vcf_anno} OUTPUT=${eh_vcf_anno}.rename.vcf NEW_SAMPLE_NAME=${id}
// 			vcfbreakmulti ${eh_vcf_anno}.rename.vcf > ${group}.expansionhunter.vcf.tmp
// 			familyfy_str.pl --vcf ${group}.expansionhunter.vcf.tmp --mother $mother --father $father --out ${group}.expansionhunter.vcf
// 			bgzip ${group}.expansionhunter.vcf
// 			tabix ${group}.expansionhunter.vcf.gz
// 			echo "STR	${params.accessdir}/vcf/${group}.expansionhunter.vcf.gz" > ${group}_str.INFO

// 			${vcfbreakmulti_expansionhunter_version(task)}
// 			"""
// 		}
// 		else {
// 			"""
// 			java -jar /opt/conda/envs/CMD-WGS/share/picard-2.21.2-1/picard.jar RenameSampleInVcf INPUT=${eh_vcf_anno} OUTPUT=${eh_vcf_anno}.rename.vcf NEW_SAMPLE_NAME=${id}
// 			vcfbreakmulti ${eh_vcf_anno}.rename.vcf > ${group}.expansionhunter.vcf
// 			bgzip ${group}.expansionhunter.vcf
// 			tabix ${group}.expansionhunter.vcf.gz
// 			echo "STR	${params.accessdir}/vcf/${group}.expansionhunter.vcf.gz" > ${group}_str.INFO

// 			${vcfbreakmulti_expansionhunter_version(task)}
// 			"""

// 		}

// 	stub:
// 		"""
// 		touch "${group}.expansionhunter.vcf.gz"
// 		touch "${group}_str.INFO"

// 		${vcfbreakmulti_expansionhunter_version(task)}
// 		"""
// }
// def vcfbreakmulti_expansionhunter_version(task) {
// 	"""
// 	cat <<-END_VERSIONS > ${task.process}_versions.yml
// 	${task.process}:
// 	    vcflib: 1.0.9
// 	    rename-sample-in-vcf: \$(echo \$(java -jar /opt/conda/envs/CMD-WGS/share/picard-2.21.2-1/picard.jar RenameSampleInVcf --version 2>&1) | sed 's/-SNAPSHOT//')
// 	    tabix: \$(echo \$(tabix --version 2>&1) | sed 's/^.*(htslib) // ; s/ Copyright.*//')
// 	END_VERSIONS
// 	"""
// }

// //////////////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////////////
// /////////////////////////////////////////////////////////////////////////
// process melt_qc_val {
// 	tag "$id"
// 	time '20m'
// 	memory '50 MB'
// 	input:
// 		tuple val(group), val(id), qc

// 	output:
// 		tuple id, val(INS_SIZE), val(MEAN_DEPTH), val(COV_DEV), emit: qc_melt_val
// 		tuple val(group), val(id), val(INS_SIZE), val(MEAN_DEPTH), val(COV_DEV), emit: qc_cnvkit_val


// 	when:
// 		params.run_melt

// 	script:
// 		// Collect qc-data if possible
// 		def ins_dev
// 		def coverage
// 		def ins_size
// 		qc.readLines().each{
// 			if (it =~ /\"(ins_size_dev)\" : \"(\S+)\"/) {
// 				ins_dev = it =~ /\"(ins_size_dev)\" : \"(\S+)\"/
// 			}
// 			if (it =~ /\"(mean_coverage)\" : \"(\S+)\"/) {
// 				coverage = it =~ /\"(mean_coverage)\" : \"(\S+)\"/
// 			}
// 			if (it =~ /\"(ins_size)\" : \"(\S+)\"/) {
// 				ins_size = it =~ /\"(ins_size)\" : \"(\S+)\"/
// 			}
// 		}
// 		// might need to be defined for -resume to work "def INS_SIZE" and so on....
// 		INS_SIZE = ins_size[0][2]
// 		MEAN_DEPTH = coverage[0][2]
// 		COV_DEV = ins_dev[0][2]
// 		"""
// 		echo hej > hej
// 		"""

// 	stub:
// 		INS_SIZE = 0
// 		MEAN_DEPTH = 0
// 		COV_DEV = 0
// 		"""
// 		echo $INS_SIZE $MEAN_DEPTH $COV_DEV > qc.val
// 		"""
// }

// // MELT always give VCFs for each type of element defined in mei_list
// // If none found -> 0 byte vcf. merge_melt.pl merges the three, if all empty
// // it creates a vcf with only header
// // merge_melt.pl gives output ${id}.melt.merged.vcf
// process melt {
// 	cpus 3
// 	errorStrategy 'retry'
// 	container  "${params.container_melt}"
// 	tag "$id"
// 	// memory seems to scale with less number of reads?
// 	memory '70 GB'
// 	time '3h'
// 	publishDir "${params.results_output_dir}/vcf", mode: 'copy' , overwrite: 'true', pattern: '*.vcf'

// 	input:
// 		tuple val(id), val(group), path(bam), path(bai), val(INS_SIZE), val(MEAN_DEPTH), val(COV_DEV)

// 	output:
// 		tuple val(group), val(id), path("${id}.melt.merged.vcf"), emit: melt_vcf_nonfiltered
// 		path "*versions.yml", emit: versions

// 	when:
// 		params.run_melt

// 	script:
// 		"""
// 		java -jar /opt/MELTv2.2.2/MELT.jar Single \\
// 			-bamfile $bam \\
// 			-r 150 \\
// 			-h ${params.genome_file} \\
// 			-n /opt/MELTv2.2.2/add_bed_files/Hg38/Hg38.genes.bed \\
// 			-z 500000 \\
// 			-d 50 \\
// 			-t $params.mei_list \\
// 			-w . \\
// 			-c $MEAN_DEPTH \\
// 			-e $INS_SIZE \\
// 			-exome
// 		merge_melt.pl $params.meltheader $id

// 		${melt_version(task)}
// 		"""

// 	stub:
// 		"""
// 		touch "${id}.melt.merged.vcf"
// 		${melt_version(task)}
// 		"""
// }
// def melt_version(task) {
// 	"""
// 	cat <<-END_VERSIONS > ${task.process}_versions.yml
// 	${task.process}:
// 	    melt: \$(echo \$(java -jar /opt/MELTv2.2.2/MELT.jar -h | grep "^MELTv" | cut -f1 -d" " | sed "s/MELTv//" ) )
// 	END_VERSIONS
// 	"""
// }

// process intersect_melt {
// 	cpus 2
// 	tag "$id"
// 	memory '2 GB'
// 	time '1h'
// 	publishDir "${params.results_output_dir}/vcf", mode: 'copy' , overwrite: 'true', pattern: '*.vcf'

// 	input:
// 		tuple val(group), val(id), path(vcf)

// 	output:
// 		tuple val(group), val(id), path("${id}.melt.merged.intersected.vcf"), emit: ch_melt_vcf
// 		path "*versions.yml", emit: versions

// 	when:
// 		params.run_melt

// 	script:
// 		"""
// 		bedtools intersect -a $vcf -b $params.intersect_bed -header > ${id}.melt.merged.intersected.vcf
// 		${intersect_melt_version(task)}
// 		"""

// 	stub:
// 		"""
// 		touch "${id}.melt.merged.intersected.vcf"
// 		${intersect_melt_version(task)}
// 		"""
// }
// def intersect_melt_version(task) {
// 	"""
// 	cat <<-END_VERSIONS > ${task.process}_versions.yml
// 	${task.process}:
// 	    bedtools: \$(echo \$(bedtools --version 2>&1) | sed -e "s/^.*bedtools v//" )
// 	END_VERSIONS
// 	"""
// }


// process bamtoyaml {
// 	cpus 1
// 	time "5m"
// 	memory "2MB"

// 	input:
// 		tuple val(group), val(id), bam, bai

// 	output:
// 		tuple val(group), path("${group}_bamstart.INFO"), emit: bamchoice_INFO

// 	script:
// 		"""
// 		echo "BAM	$id	/access/${params.subdir}/bam/${bam.getName()}" > ${group}_bamstart.INFO
// 		"""

// 	stub:
// 		"""
// 		touch "${group}_bamstart.INFO"
// 		"""
// }


process gvcf_combine {
	cpus 16
	tag "$group"
	memory '5 GB'
	time '5h'
	container  "${params.container_sentieon}"

	input:
		tuple val(group), val(id), path(gvcfs), path(gvcf_idxs)

	output: // Off to split_normalize, together with other stuff
		tuple val(group), val(id), path("${group}.combined.vcf"), path("${group}.combined.vcf.idx"), emit: combined_vcf
		path "*versions.yml", emit: versions

	script:
		all_gvcfs = gvcfs.collect { it.toString() }.sort().join(' -v ')
		println(all_gvcfs)
		"""
		sentieon driver \\
			-t ${task.cpus} \\
			-r ${params.genome_file} \\
			--algo GVCFtyper \\
			-v $all_gvcfs ${group}.combined.vcf

		${gvcf_combine_version(task)}
		"""

	stub:
		all_gvcfs = gvcfs.collect { it.toString() }.sort().join(' -v ')
		println(all_gvcfs)
		"""
		touch "${group}.combined.vcf"
		touch "${group}.combined.vcf.idx"

		${gvcf_combine_version(task)}
		"""
}
def gvcf_combine_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
	END_VERSIONS
	"""
}

// Create ped
process create_ped {
	tag "$group"
	time '20m'
	publishDir "${params.results_output_dir}/ped", mode: 'copy' , overwrite: 'true'
	memory '1 GB'

	input:
	tuple val(group), val(id), val(type), val(sex), val(mother), val(father)

	output:
		tuple val(group), val(type), path("${group}_base.ped"), emit: ped_base
		tuple val(group), val(type_ma), path("${group}_ma.ped"), emit: ped_ma, optional: true
		tuple val(group), val(type_fa), path("${group}_fa.ped"), emit: ped_fa, optional: true

	script:
		if ( father == "" ) {
			father = "0"
		}
		if ( mother == "" ) {
			mother = "0"
		}
		type_fa = "fa"
		type_ma = "ma"
		"""
		create_ped.pl --mother $mother --father $father --group $group --id $id --sex $sex
		"""

	stub:
		type_fa = "fa"
		type_ma = "ma"
		"""
		touch "${group}_base.ped"
		touch "${group}_ma.ped"
		touch "${group}_fa.ped"

        echo $type_fa $type_ma > type.val
		"""
}

//madeline ped, run if family mode
process madeline {
	publishDir "${params.results_output_dir}/ped", mode: 'copy' , overwrite: 'true', pattern: '*.xml'
	memory '1 GB'
	time '1h'
	cpus 2
	container  "${params.container_madeline}"

	input:
		tuple val(group), val(type), path(ped)

	output:
		path("${ped}.madeline.xml")
		tuple val(group), path("${group}_madde.INFO"), emit: madde_INFO
		path "*versions.yml", emit: versions

	script:
		"""
		source activate tools
		ped_parser \\
			-t ped $ped \\
			--to_madeline \\
			-o ${ped}.madeline
		madeline2 \\
			-L "IndividualId" ${ped}.madeline \\
			-o ${ped}.madeline \\
			-x xml
		echo "MADDE	$type ${params.accessdir}/ped/${ped}.madeline.xml" > ${group}_madde.INFO

		${madeline_version(task)}
		"""

	stub:
		"""
		source activate tools
		touch "${group}_madde.INFO"
		touch "${ped}.madeline.xml"

		${madeline_version(task)}
		"""
}
def madeline_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    ped-parser: \$(echo \$(ped_parser --version 2>&1) | sed -e "s/^.*ped_parser version: //")
	    madeline: \$(echo \$(madeline2 --version 2>&1) | grep : | sed -e"s/^.*Madeline //; s/PDE : 1.*//")
	END_VERSIONS
	"""
}

process freebayes {
	cpus 1
	time '2h'
	memory '10 GB'
	container  "${params.container_twist_myeloid}"

	input:
		tuple val(group), val(id), path(bam), path(bai)

	output:
		tuple val(group), path("${id}.pathfreebayes.vcf_no_header.tsv"), emit: freebayes_variants
		path "*versions.yml", emit: versions


	when:
		params.antype == "panel"


	script:
		if (params.onco) {
			"""
			freebayes -f ${params.genome_file} --pooled-continuous --pooled-discrete -t $params.intersect_bed --min-repeat-entropy 1 -F 0.03 $bam > ${id}.freebayes.vcf
			vcfbreakmulti ${id}.freebayes.vcf > ${id}.freebayes.multibreak.vcf
			bcftools norm -m-both -c w -O v -f ${params.genome_file} -o ${id}.freebayes.multibreak.norm.vcf ${id}.freebayes.multibreak.vcf
			vcfanno_linux64 -lua $params.VCFANNO_LUA $params.vcfanno ${id}.freebayes.multibreak.norm.vcf > ${id}.freebayes.multibreak.norm.anno.vcf
			grep ^# ${id}.freebayes.multibreak.norm.anno.vcf > ${id}.freebayes.multibreak.norm.anno.path.vcf
			grep -v ^# ${id}.freebayes.multibreak.norm.anno.vcf | grep -i pathogenic > ${id}.freebayes.multibreak.norm.anno.path.vcf2
			cat ${id}.freebayes.multibreak.norm.anno.path.vcf ${id}.freebayes.multibreak.norm.anno.path.vcf2 > ${id}.freebayes.multibreak.norm.anno.path.vcf3
			filter_freebayes.pl ${id}.freebayes.multibreak.norm.anno.path.vcf3 > "${id}.pathfreebayes.vcf_no_header.tsv"

			${freebayes_version(task)}
			"""
		}
		else {
			"""
			touch "${id}.pathfreebayes.vcf_no_header.tsv"

			${freebayes_version(task)}
			"""
		}

	stub:
		"""
		touch "${id}.pathfreebayes.vcf_no_header.tsv"

		${freebayes_version(task)}
		"""
}
def freebayes_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    freebayes: \$(echo \$(freebayes --version 2>&1) | sed 's/version:\s*v//g')
	    vcflib: 1.0.9
	    bcftools: \$(echo \$(bcftools --version 2>&1) | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
	    vcfanno: \$(echo \$(vcfanno_linux64 2>&1 | grep version | cut -f3 -d' ')  )
	END_VERSIONS
	"""
}

/////////////// MITOCHONDRIA SNV CALLING ///////////////
///////////////                          ///////////////

// create an MT BAM file
process fetch_MTseqs {
	cpus 2
	memory '10GB'
	time '1h'
	tag "$id"
	publishDir "${params.results_output_dir}/bam", mode: 'copy', overwrite: 'true', pattern: '*.bam*'

	input:
		tuple val(group), val(id), path(bam), path(bai)

    output:
        tuple val(group), val(id), file ("${id}_mito.bam"), path("${id}_mito.bam.bai"), emit: bam_bai
		tuple val(group), path("${group}_mtbam.INFO"), emit: mtBAM_INFO
		path "*versions.yml", emit: versions

	script:
		"""
		sambamba view -f bam $bam M > ${id}_mito.bam
		samtools index -b ${id}_mito.bam
		echo "mtBAM	$id	/access/${params.subdir}/bam/${id}_mito.bam" > ${group}_mtbam.INFO

		${fetch_MTseqs_version(task)}
		"""

	stub:
		"""
		touch "${id}_mito.bam"
		touch "${id}_mito.bam.bai"
		touch "${group}_mtbam.INFO"

		${fetch_MTseqs_version(task)}
		"""
}
def fetch_MTseqs_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    sambamba: \$(echo \$(sambamba --version 2>&1) | awk '{print \$2}' )
	    samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
	END_VERSIONS
	"""
}


process sentieon_mitochondrial_qc {

    // Fetch mitochondrial coverage statistics
    // Calculate mean_coverage and pct_above_500x

    cpus 30
    memory '20 GB'
	tag "$id"
	time '2h'
	container  "${params.container_sentieon}"

	input:
        tuple val(group), val(id), path(bam), path(bai)

	output:
    	tuple val(group), val(id), path("${id}_mito_coverage.tsv"), emit: qc_tsv
		path "*versions.yml", emit: versions

	when:
	    params.antype == "wgs"

	script:
		"""
		sentieon driver \\
			-r ${params.genome_file} \\
			-t ${task.cpus} \\
			-i $bam \\
			--algo CoverageMetrics \\
			--cov_thresh 500 \\
			mt_cov_metrics.txt

		head -1 mt_cov_metrics.txt.sample_interval_summary > "${id}_mito_coverage.tsv"
		grep "^M" mt_cov_metrics.txt.sample_interval_summary >> "${id}_mito_coverage.tsv"
		${sentieon_mitochondrial_qc_version(task)}
		"""

	stub:
		"""
		touch "${id}_mito_coverage.tsv"
		${sentieon_mitochondrial_qc_version(task)}
		"""
}
def sentieon_mitochondrial_qc_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
	END_VERSIONS
	"""
}

process build_mitochondrial_qc_json {
    memory '1 GB'
    cpus 2
    tag "$id"
    time "1h"

    input:
        tuple val(group), val(id), path(mito_qc_file)
    output:
        tuple val(group), val(id), path("${id}_mito_qc.json"), emit: qc_json

	script:
		"""
		mito_tsv_to_json.py ${mito_qc_file} > "${id}_mito_qc.json"
		"""
	stub:
		"""
		touch "${id}_mito_qc.json"
		"""
}


// gatk FilterMutectCalls in future if FPs overwhelms tord/sofie/carro
process run_mutect2 {
	cpus 4
	memory '50 GB'
	time '1h'
	tag "$group"
	publishDir "${params.results_output_dir}/vcf", mode: 'copy', overwrite: 'true', pattern: '*.vcf'

	input:
		tuple val(group), val(id), path(bam), path(bai)

	output:
		tuple val(group), val(id), path("${group}.mutect2.vcf"), emit: vcf
		path "*versions.yml", emit: versions

	when:
		!params.onco

	script:
		bams = bam.join(' -I ')

		"""
		source activate gatk4-env
		gatk Mutect2 \
		--mitochondria-mode \
		-R $params.genome_file \
		-L M \
		-I $bams \
		-O ${group}.mutect2.vcf

		${run_mutect2_version(task)}
		"""

	stub:
		bams = bam.join(' -I ')
		println(bams)
		"""
		source activate gatk4-env
		touch "${group}.mutect2.vcf"

		${run_mutect2_version(task)}
		"""
}
def run_mutect2_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    gatk: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$// ; s/-SNAPSHOT//')
	END_VERSIONS
	"""
}

// split and left-align variants
process split_normalize_mito {
	cpus 2
	memory '1GB'
	time '1h'

	input:
		tuple val(group), val(id), path(mito_snv_vcf)
		tuple val(g2), val(id2), val(sex), val(type)

	output:
		tuple val(group), path("${group}.mutect2.breakmulti.filtered5p.0genotyped.proband.vcf"), emit: vcf
		path "*versions.yml", emit: versions

	script:
		proband_idx = type.findIndexOf{ it == "proband" }

		"""
		grep -vP "^M\\s+955" $mito_snv_vcf > ${mito_snv_vcf}.fix
		bcftools norm -m-both -o ${mito_snv_vcf}.breakmulti ${mito_snv_vcf}.fix
		bcftools sort ${mito_snv_vcf}.breakmulti | bgzip > ${mito_snv_vcf}.breakmulti.fix
		tabix -p vcf ${mito_snv_vcf}.breakmulti.fix
		bcftools norm -f $params.rCRS_fasta -o ${mito_snv_vcf.baseName}.adjusted.vcf ${mito_snv_vcf}.breakmulti.fix
		bcftools view -i 'FMT/AF[*]>0.05' ${mito_snv_vcf.baseName}.adjusted.vcf -o ${group}.mutect2.breakmulti.filtered5p.vcf
		bcftools filter -S 0 --exclude 'FMT/AF[*]<0.05' ${group}.mutect2.breakmulti.filtered5p.vcf -o ${group}.mutect2.breakmulti.filtered5p.0genotyped.vcf
		filter_mutect2_mito.pl ${group}.mutect2.breakmulti.filtered5p.0genotyped.vcf ${id2[proband_idx]} > ${group}.mutect2.breakmulti.filtered5p.0genotyped.proband.vcf

		${split_normalize_mito_version(task)}
		"""

	stub:
		"""
		touch "${group}.mutect2.breakmulti.filtered5p.0genotyped.proband.vcf"
		${split_normalize_mito_version(task)}
		"""
}
def split_normalize_mito_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    vcflib: 1.0.9
	    bcftools: \$(echo \$(bcftools --version 2>&1) | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
	    tabix: \$(echo \$(tabix --version 2>&1) | sed 's/^.*(htslib) // ; s/ Copyright.*//')
	END_VERSIONS
	"""
}

// use python tool HmtNote for annotating vcf
// future merging with diploid genome does not approve spaces in info-string
// TODO: what is this future merging issue and does it still apply?
process run_hmtnote {
	cpus 2
	memory '5GB'
	time '1h'

	input:
		tuple val(group), path(vcf)

	output:
		tuple val(group), path("${group}.fixinfo.vcf"), emit: vcf
		path "*versions.yml", emit: versions

	script:
		"""
		source activate tools
		hmtnote annotate ${vcf} ${group}.hmtnote --offline
		grep ^# ${group}.hmtnote > ${group}.fixinfo.vcf
		grep -v ^# ${group}.hmtnote | sed 's/ /_/g' >> ${group}.fixinfo.vcf

		${run_hmtnote_version(task)}
		"""

	stub:
		"""
		source activate tools
		touch "${group}.fixinfo.vcf"

		${run_hmtnote_version(task)}
		"""
}
def run_hmtnote_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    hmtnote: \$(echo \$(hmtnote --version 2>&1) | sed 's/^.*hmtnote, version //; s/Using.*\$//' )
	END_VERSIONS
	"""
}

// run haplogrep 2 on resulting vcf
process run_haplogrep {
	time '1h'
	memory '50 GB'
	cpus 2
	publishDir "${params.results_output_dir}/plots/mito", mode: 'copy', overwrite: 'true', pattern: '*.png'

	input:
		tuple val(group), val(id), path(mito_snv_vcf)

	output:
		path "${group}.haplogrep.png"
		tuple val(group), path("${group}_haplo.INFO"), emit: haplogrep_INFO
		path "*versions.yml", emit: versions

	script:
		version_str = run_haplogrep_version(task)
		"""
		for sample in \$(bcftools query -l "${mito_snv_vcf}"); do

			bcftools view -c1 -Oz -s "\$sample" -o "\${sample}.vcf.gz" "${mito_snv_vcf}"
			java  -Xmx16G -Xms16G -jar /opt/bin/haplogrep.jar classify \
			--in "\${sample}.vcf.gz" \\
			--out "\${sample}.hg2.vcf" \\
			--format vcf \\
			--lineage 1

			dot "\${sample}.hg2.vcf.dot" -Tps2 > "\${sample}.hg2.vcf.ps"

			gs -dSAFER -dBATCH -dNOPAUSE -sDEVICE=png16m -dGraphicsAlphaBits=4 -r1200 -dDownScaleFactor=3 -sOutputFile=\${sample}.hg2.vcf.png \${sample}.hg2.vcf.ps

		done
		montage -mode concatenate -tile 3x1 *.png ${group}.haplogrep.png
		echo "IMG haplogrep ${params.accessdir}/plots/mito/${group}.haplogrep.png" > "${group}_haplo.INFO"

		echo "${version_str}" > "${task.process}_versions.yml"
		"""

	stub:
		version_str = run_haplogrep_version(task)
		"""
		touch "${group}.haplogrep.png"
		touch "${group}_haplo.INFO"

		echo "${version_str}" > "${task.process}_versions.yml"
		"""
}
def run_haplogrep_version(task) {
	// This docstring looks different
	// If spaces similarly to the others, this leads to additional whitespace above and below the version text
	"""${task.process}:
	    haplogrep: \$(echo \$(java -jar /opt/bin/haplogrep.jar classify 2>&1) | sed "s/htt.*Classify v// ; s/ .*//")
	    montage: \$(echo \$(gm -version 2>&1) | head -1 | sed -e "s/GraphicsMagick //" | cut -d" " -f1 )"""
}

// use eKLIPse for detecting mitochondrial deletions
process run_eklipse {

	tag "$id"
	cpus 2
	// in rare cases with samples above 50 000x this can peak at 500+ GB of VMEM. Add downsampling!
	memory '100GB'
	time '60m'
	publishDir "${params.results_output_dir}/plots/mito", mode: 'copy', overwrite: 'true', pattern: '*.txt'
	publishDir "${params.results_output_dir}/plots/mito", mode: 'copy', overwrite: 'true', pattern: '*.png'

	input:
		tuple val(group), val(id), path(bam), path(bai)
		tuple val(group2), val(id2), val(sex), val(type)

	output:
		tuple path("*.png"), path("${id}.hetplasmid_frequency.txt")
		tuple val(group), path("${id}_eklipse.INFO"), emit: eklipse_INFO, optional: true
		path "*versions.yml", emit: versions

	script:
		yml_info_command = ""
		if (type == "proband") {
			yml_info_command = "echo 'IMG eklipse ${params.accessdir}/plots/mito/${id}_eklipse.png' > ${id}_eklipse.INFO"
		}
		"""
		source activate htslib10
		echo "${bam}\tsample" > infile.txt
		python /eKLIPse/eKLIPse.py \
		-in infile.txt \
		-ref /eKLIPse/data/NC_012920.1.gb
		mv eKLIPse_*/eKLIPse_deletions.csv ./${id}_deletions.csv
		mv eKLIPse_*/eKLIPse_genes.csv ./${id}_genes.csv
		mv eKLIPse_*/eKLIPse_sample.png ./${id}_eklipse.png
		hetplasmid_frequency_eKLIPse.pl --bam ${bam} --in ${id}_deletions.csv
		mv hetplasmid_frequency.txt ${id}.hetplasmid_frequency.txt
		$yml_info_command

		${run_eklipse_version(task)}
		"""

	stub:
		yml_info_command = ""
		if (type == "proband") {
			yml_info_command = "echo 'IMG eklipse ${params.accessdir}/plots/mito/${id}_eklipse.png' > ${id}_eklipse.INFO"
		}
		"""
		source activate htslib10
		touch "${id}.hetplasmid_frequency.txt"
		touch "${id}.png"
		$yml_info_command

		${run_eklipse_version(task)}
		"""
}
def run_eklipse_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    eklipse: 1.8
	END_VERSIONS
	"""
}

//eklipseM_INFO.collectPath(name: "eklipse.INFO").set{ eklipse_INFO }

// Splitting & normalizing variants, merging with Freebayes/Mutect2, intersecting against exome/clinvar introns
process split_normalize {
	cpus 2
	publishDir "${params.results_output_dir}/vcf", mode: 'copy', overwrite: 'true', pattern: '*.vcf'
	tag "$group"
	memory '50 GB'
	time '1h'
	input:
		tuple val(group), val(ids), path(vcf), path(idx) // is ids supposed to be tuple?
		tuple val(group2), path(vcfconcat)

	output:
		tuple val(group), path("${group}.norm.uniq.DPAF.vcf"), emit: norm_uniq_dpaf_vcf
		tuple val(group), val(id), path("${group}.intersected.vcf"), emit: intersected_vcf
		path "*versions.yml", emit: versions

	script:
	id = ids[0]
	// rename M to MT because genmod does not recognize M
	if (params.onco || params.assay == "modycf") {
		"""
		cat $vcf $vcfconcat > ${id}.concat.freebayes.vcf
		vcfbreakmulti ${id}.concat.freebayes.vcf > ${group}.multibreak.vcf
		bcftools norm -m-both -c w -O v -f ${params.genome_file} -o ${group}.norm.vcf ${group}.multibreak.vcf
		bcftools sort ${group}.norm.vcf | vcfuniq > ${group}.norm.uniq.vcf
		wgs_DPAF_filter.pl ${group}.norm.uniq.vcf > ${group}.norm.uniq.DPAF.vcf
		bedtools intersect \
			-a ${group}.norm.uniq.DPAF.vcf \\
			-b $params.intersect_bed \\
			-u -header > ${group}.intersected.vcf

		${split_normalize_version(task)}
		"""
	}

	else {
		"""
		vcfbreakmulti ${vcf} > ${group}.multibreak.vcf
		bcftools norm -m-both -c w -O v -f ${params.genome_file} -o ${group}.norm.vcf ${group}.multibreak.vcf
		bcftools sort ${group}.norm.vcf | vcfuniq > ${group}.norm.uniq.vcf
		wgs_DPAF_filter.pl ${group}.norm.uniq.vcf > ${group}.norm.uniq.DPAF.vcf
		bedtools intersect \\
			-a ${group}.norm.uniq.DPAF.vcf \\
			-b $params.intersect_bed \\
			-u -header > ${group}.intersected_diploid.vcf
		java -jar /opt/conda/envs/CMD-WGS/share/picard-2.21.2-1/picard.jar MergeVcfs \
		I=${group}.intersected_diploid.vcf I=$vcfconcat O=${group}.intersected.vcf
		sed 's/^M/MT/' -i ${group}.intersected.vcf
		sed 's/ID=M,length/ID=MT,length/' -i ${group}.intersected.vcf

		${split_normalize_version(task)}
		"""
	}

	stub:
	id = ids[0]
		"""
		touch "${group}.norm.uniq.DPAF.vcf"
		touch "${group}.intersected.vcf"
		touch "${group}.multibreak.vcf"

		${split_normalize_version(task)}
		"""
}
def split_normalize_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    vcflib: 1.0.9
	    bcftools: \$(echo \$(bcftools --version 2>&1) | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
	    bedtools: \$(echo \$(bedtools --version 2>&1) | sed -e "s/^.*bedtools v//" )
	    merge-vcfs: \$(echo \$(java -jar /opt/conda/envs/CMD-WGS/share/picard-2.21.2-1/picard.jar MergeVcfs --version 2>&1 | sed 's/-SNAPSHOT//'))
	END_VERSIONS
	"""
}

// /////////////// Collect QC, emit: single file ///////////////

// process merge_qc_json {
//     cpus 2
//     errorStrategy 'retry'
//     maxErrors 5
//     publishDir "${params.results_output_dir}/qc", mode: 'copy' , overwrite: 'true', pattern: '*.QC'
//     tag "$id"
//     time '1h'
// 	memory '1 GB'

//     input:
//         tuple val(group), val(id), path(qc)

//     output:
//         tuple id, path("${id}.QC"), emit: qc_cdm_merged

//     script:
//         qc_json_files = qc.join(' ')
// 		"""
// 		merge_json_files.py ${qc_json_files} > ${id}.QC
// 		"""

// 	stub:
// 		"""
// 		touch "${id}.QC"
// 		"""
// }

// // Load QC data, emit: CDM (via middleman)
// process qc_to_cdm {
// 	cpus 2
// 	errorStrategy 'retry'
// 	maxErrors 5
// 	publishDir "${params.crondir}/qc", mode: 'copy' , overwrite: 'true'
// 	tag "$id"
// 	time '1h'

// 	input:
// 		tuple id, path(qc), diagnosis, r1, r2

// 	output:
// 		path("${id}.cdm"), emit: cdm_done


// 	when:
// 		!params.noupload


// 	script:
// 		parts = r1.split('/')
// 		idx =  parts.findIndexOf {it ==~ /......_......_...._........../}
// 		rundir = parts[0..idx].join("/")
// 		"""
// 		echo "--run-folder $rundir --sample-id $id --subassay $diagnosis --assay $params.assay --qc ${params.results_output_dir}/qc/${id}.QC" > ${id}.cdm
// 		"""
// }


process annotate_vep {
	container  "${params.container_vep}"
	cpus 30
	tag "$group"
	memory '50 GB'
	time '5h'

	input:
		tuple val(group), val(id), path(vcf)

	output:
		tuple val(group), path("${group}.vep.vcf"), emit: vcf
		path "*versions.yml", emit: versions

	script:
		"""
		vep \\
			-i ${vcf} \\
			-o ${group}.vep.vcf \\
			--offline \\
			--sift b --polyphen b --ccds --hgvs --symbol --numbers --domains --regulatory --canonical --protein --biotype --af --af_1kg --max_af --pubmed --uniprot --mane --tsl --appris --variant_class --gene_phenotype --mirna \\
			--merged \\
			--vcf \\
			--no_stats \\
			--synonyms $params.VEP_SYNONYMS \\
			--fork ${task.cpus} \\
			--force_overwrite \\
			--fasta $params.VEP_FASTA \\
			--dir_cache $params.VEP_CACHE \\
			--dir_plugins $params.VEP_PLUGINS \\
			--distance $params.VEP_TRANSCRIPT_DISTANCE \\
			-cache \\
			--plugin CADD,$params.CADD \\
			--plugin LoFtool \\
			--plugin MaxEntScan,$params.MAXENTSCAN,SWA,NCSS \\
			--plugin dbNSFP,$params.DBNSFP,transcript_match=1,REVEL_score,REVEL_rankscore \\
			-custom $params.GNOMAD_EXOMES,gnomADe,vcf,exact,0,AF_popmax,AF,popmax \\
			-custom $params.GNOMAD_GENOMES,gnomADg,vcf,exact,0,AF_popmax,AF,popmax \\
			-custom $params.GNOMAD_MT,gnomAD_mt,vcf,exact,0,AF_hom,AF_het \\
			-custom $params.PHYLOP,phyloP100way,bigwig \\
			-custom $params.PHASTCONS,phastCons,bigwig

		${annotate_vep_version(task)}
		"""

	stub:
		"""
		touch "${group}.vep.vcf"
		${annotate_vep_version(task)}
		"""
}
def annotate_vep_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    vep: \$( echo \$(vep --help 2>&1) | sed 's/^.*Versions:.*ensembl-vep : // ; s/ .*\$//')
	END_VERSIONS
	"""
}


// // gene, clinvar, loqusdb, enigma(onco)
process vcfanno {
	memory '1GB'
	time '1h'
	errorStrategy 'retry'
	maxErrors 5
	cpus 2

	input:
		tuple val(group), path(vcf)

	output:
		tuple val(group), path("${group}.clinvar.loqusdb.gene.vcf"), emit: vcf
		path "*versions.yml", emit: versions

	script:
		"""
		vcfanno_linux64 -lua $params.VCFANNO_LUA $params.vcfanno $vcf > ${group}.clinvar.loqusdb.gene.vcf
		${vcfanno_version(task)}
		"""

	stub:
		"""
		touch "${group}.clinvar.loqusdb.gene.vcf"
		${vcfanno_version(task)}
		"""
}
def vcfanno_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    vcfanno: \$(echo \$(vcfanno_linux64 2>&1 | grep version | cut -f3 -d' ')  )
	END_VERSIONS
	"""
}


// // Extracting most severe consequence:
// // Modifying annotations by VEP-plugins, and adding to info-field:
// // Modifying CLNSIG field to allow it to be used by genmod score properly:
// TODO: give process better name
process modify_vcf {
	cpus 2
	tag "$group"
	memory '1 GB'
	time '1h'

	input:
		tuple val(group), path(vcf)

	output:
		tuple val(group), path("${group}.mod.vcf"), emit: vcf

	script:
		"""
		modify_vcf_scout.pl $vcf > ${group}.mod.vcf
		"""

	stub:
		"""
		touch "${group}.mod.vcf"
		"""
}


// Marking splice INDELs:
process mark_splice {
	cpus 2
	tag "$group"
	memory '1 GB'
	time '1h'

	input:
		tuple val(group), path(vcf)

	output:
		tuple val(group), path("${group}.marksplice.vcf"), emit: splice_marked

	script:
		"""
		/opt/bin/mark_spliceindels.pl $vcf > ${group}.marksplice.vcf
		"""

	stub:
		"""
		touch "${group}.marksplice.vcf"
		"""
}

// Extract all INDELs
process extract_indels_for_cadd {
	cpus 2
	tag "$group"
	memory '1 GB'
	time '1h'

	input:
		tuple val(group), val(id), path(vcf)

	output:
		tuple val(group), path("${group}.only_indels.vcf"), emit: vcf
		path "*versions.yml", emit: versions

	script:
		"""
		bcftools view $vcf -V snps -o ${group}.only_indels.vcf
		${extract_indels_for_cadd_version(task)}
		"""

	stub:
		"""
		touch "${group}.only_indels.vcf"
		${extract_indels_for_cadd_version(task)}
		"""
}
def extract_indels_for_cadd_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    bcftools: \$(echo \$(bcftools --version 2>&1) | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
	END_VERSIONS
	"""
}

// Annotate Indels with VEP+Gnomad genomes. Filter variants below threshold
process indel_vep {
	cpus 5
	container  "${params.container_vep}"
	tag "$group"
	memory '10 GB'
	time '3h'

	input:
		tuple val(group), path(vcf)

	output:
		tuple val(group), path("${group}.only_indels.vep.filtered.vcf"), emit: vcf
		path "*versions.yml", emit: versions

	script:
		"""
		vep \\
			-i $vcf \\
			-o ${group}.only_indels.vep.vcf \\
			--offline \\
			--cache \\
			--merged \\
			--vcf \\
			--synonyms $params.VEP_SYNONYMS \\
			--fasta $params.VEP_FASTA \\
			-custom $params.GNOMAD_GENOMES,gnomADg,vcf,exact,0,AF \\
			-custom $params.GNOMAD_MT,gnomAD_mt,vcf,exact,0,AF_hom,AF_het \\
			--dir_cache $params.VEP_CACHE \\
			--force_overwrite \\
			--no_stats \\
			--fork ${task.cpus}
		filter_indels.pl ${group}.only_indels.vep.vcf > ${group}.only_indels.vep.filtered.vcf
		${indel_vep_version(task)}
		"""

	stub:
		"""
		touch "${group}.only_indels.vep.filtered.vcf"
		${indel_vep_version(task)}
		"""
}
def indel_vep_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    vep: \$( echo \$(vep --help 2>&1) | sed 's/^.*Versions:.*ensembl-vep : // ; s/ .*\$//')
	END_VERSIONS
	"""
}

// Calculate CADD scores for all indels
process calculate_indel_cadd {
	cpus 2
	container  "${params.container_cadd}"
	tag "$group"
	memory '20 GB'
	time '3h'

	input:
		tuple val(group), path(vcf)

	output:
		tuple val(group), path("${group}.indel_cadd.gz"), emit: cadd_gz
		path "*versions.yml", emit: versions

	//TODO: does this output a VCF? if so fix file ending.
	script:
		"""
		/CADD-scripts/CADD.sh -c ${task.cpus} -g GRCh38 -o ${group}.indel_cadd.gz $vcf
		${calculate_indel_cadd_version(task)}
		"""

	stub:
		"""
		touch "${group}.indel_cadd.gz"
		${calculate_indel_cadd_version(task)}
		"""
}
def calculate_indel_cadd_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    cadd: \$(echo \$(/CADD-scripts/CADD.sh -v 2>&1) | sed -e "s/^.*CADD-v// ; s/ (c).*//")
	END_VERSIONS
	"""
}

process bgzip_indel_cadd {

	tag "$group"
	cpus 4
	memory '1 GB'
	time '5m'
	container  "${params.container_bcftools}"

	input:
		tuple val(group), path(cadd_scores)

	output:
		tuple val(group), path("${group}.cadd.gz"), path("${group}.cadd.gz.tbi"), emit: cadd_tbi
		path "*versions.yml", emit: versions

	script:
		"""
		gunzip -c ${cadd_scores} > "${group}.cadd"
		bgzip -@ ${task.cpus} "${group}.cadd"
		tabix -p vcf "${group}.cadd.gz"
		${bgzip_indel_cadd_version(task)}
		"""

	stub:
		"""
		touch "${group}.cadd.gz"
		touch "${group}.cadd.gz.tbi"
		${bgzip_indel_cadd_version(task)}
		"""
}
def bgzip_indel_cadd_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    tabix: \$(echo \$(tabix --version 2>&1) | sed 's/^.*(htslib) // ; s/ Copyright.*//')
	    bcftools: \$(echo \$(bcftools --version 2>&1) | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
	END_VERSIONS
	"""
}

// Add the calculated indel CADDs to the vcf
process add_cadd_scores_to_vcf {
	cpus 4
	tag "$group"
	memory '1 GB'
	time '5m'
	container  "${params.container_genmod}"

	input:
		tuple val(group), path(vcf), path(cadd_scores), path(cadd_scores_tbi)

	output:
		tuple val(group), path("${group}.cadd.vcf"), emit: vcf
		path "*versions.yml", emit: versions

	script:
		"""
		genmod annotate --cadd-file ${cadd_scores} ${vcf} > ${group}.cadd.vcf

		${add_cadd_scores_to_vcf_version(task)}
		"""

	stub:
		"""
		touch "${group}.cadd.vcf"
		${add_cadd_scores_to_vcf_version(task)}
		"""
}
def add_cadd_scores_to_vcf_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    genmod: \$(echo \$(genmod --version 2>&1) | sed -e "s/^.*genmod version: //")
	END_VERSIONS
	"""
}


// # Annotating variant inheritance models:
process inher_models {
	tag "$group"
	cpus 3
	memory '80 GB'
	time '1h'
	container  "${params.container_genmod}"

	input:
		tuple val(group), path(vcf), val(type), path(ped)

	output:
		tuple val(group), val(type), path("${group}.models.vcf"), emit: vcf
		path "*versions.yml", emit: versions

	script:
		"""
		genmod models $vcf -p ${task.cpus} -f $ped > ${group}.models.vcf
		${inher_models_version(task)}
		"""

	stub:
		"""
		touch "${group}.models.vcf"
		${inher_models_version(task)}
		"""
}
def inher_models_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    genmod: \$(echo \$(genmod --version 2>&1) | sed -e "s/^.*genmod version: //")
	END_VERSIONS
	"""
}


// Scoring variants:
// Adjusting compound scores:
// Sorting VCF according to score:
process genmodscore {
	tag "$group"
	cpus 2
	memory '20 GB'
	time '1h'
	container  "${params.container_genmod}"

	input:
		tuple val(group), val(type), path(vcf)

	output:
		tuple val(group), val(type), path("${group_score}.scored.vcf"), emit: scored_vcf
		path "*versions.yml", emit: versions

	script:
		group_score = ( type == "ma" || type == "fa" ) ? "${group}_${type}" : group

		if ( params.mode == "family" && params.antype == "wgs" ) {
			"""
			genmod score -i $group_score -c $params.rank_model -r $vcf -o ${group_score}.only_rankscore.vcf
			genmod compound \
				--threshold ${params.genmod_compound_trio_threshold} \
				--penalty ${params.genmod_compound_trio_penalty} \
				-o ${group_score}.with_compounds.vcf \
				${group_score}.only_rankscore.vcf
			sed 's/RankScore=${group}:/RankScore=${group_score}:/g' -i ${group_score}.with_compounds.vcf
			genmod sort -p -f $group_score ${group_score}.with_compounds.vcf -o ${group_score}.scored.vcf

			${genmodscore_version(task)}
			"""
		}
		else {
			"""
			genmod score -i $group_score -c $params.rank_model_s -r $vcf -o ${group_score}.only_rankscore.vcf

			# To get compounds without applying rank score penalty
			genmod compound \
				--penalty 0 \
				-o ${group_score}.with_compounds.vcf \
				${group_score}.only_rankscore.vcf

			genmod sort \
				-p \
				-f $group_score \
				-o ${group_score}.scored.vcf \
				${group_score}.with_compounds.vcf

			${genmodscore_version(task)}
			"""
		}

	stub:
		group_score = group
		"""
		touch "${group_score}.scored.vcf"
		${genmodscore_version(task)}
		"""
}
def genmodscore_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    genmod: \$(echo \$(genmod --version 2>&1) | sed -e "s/^.*genmod version: //")
	END_VERSIONS
	"""
}

// Bgzipping and indexing VCF:
process vcf_completion {
	cpus 16
	publishDir "${params.results_output_dir}/vcf", mode: 'copy', overwrite: 'true', pattern: '*.vcf.gz*'
	tag "$group"
	time '1h'
	memory '5 GB'

	input:
		tuple val(group), val(type), path(vcf)

	output:
		tuple val(group), val(type), path("${group_score}.scored.vcf.gz"), path("${group_score}.scored.vcf.gz.tbi"), emit: vcf_tbi
		tuple val(group), path("${group}_snv.INFO"), emit: snv_INFO
		path "*versions.yml", emit: versions

	script:
		group_score = ( type == "ma" || type == "fa" ) ? "${group}_${type}" : group

		"""
		sed 's/^MT/M/' -i $vcf
		sed 's/ID=MT,length/ID=M,length/' -i $vcf
		bgzip -@ ${task.cpus} $vcf -f
		tabix ${vcf}.gz -f
		echo "SNV	$type	${params.accessdir}/vcf/${group_score}.scored.vcf.gz" > ${group}_snv.INFO

		${vcf_completion_version(task)}
		"""

	stub:
		group_score = ( type == "ma" || type == "fa" ) ? "${group}_${type}" : group

		"""
		touch "${group_score}.scored.vcf.gz"
		touch "${group_score}.scored.vcf.gz.tbi"
		touch "${group}_snv.INFO"

		${vcf_completion_version(task)}
		"""
}
def vcf_completion_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    tabix: \$(echo \$(tabix --version 2>&1) | sed 's/^.*(htslib) // ; s/ Copyright.*//')
	END_VERSIONS
	"""
}

process peddy {

	publishDir "${params.results_output_dir}/ped", mode: 'copy' , overwrite: 'true', pattern: '*.ped'
	publishDir "${params.results_output_dir}/ped", mode: 'copy' , overwrite: 'true', pattern: '*.csv'

	cpus 4
	tag "$group"
	time '1h'
	memory '20GB'

	input:
		tuple val(group), val(type), path(vcf), val(idx)
		tuple val(group1), val(type1), path(ped)

	output:
		tuple path("${group}.ped_check.csv"),path("${group}.peddy.ped"), path("${group}.sex_check.csv"), emit: peddy_files
		tuple val(group), path("${group}_peddy.INFO"), emit: peddy_INFO
		path "*versions.yml", emit: versions

	when:
		!params.annotate_only && params.run_peddy

	script:
		"""
		source activate py3-env
		python -m peddy --sites hg38 -p ${task.cpus} $vcf $ped --prefix $group
		echo "PEDDY	${params.accessdir}/ped/${group}.ped_check.csv,${params.accessdir}/ped/${group}.peddy.ped,${params.accessdir}/ped/${group}.sex_check.csv" > ${group}_peddy.INFO

		${peddy_version(task)}
		"""

	stub:
		"""
		source activate py3-env
		touch "${group}.ped_check.csv"
		touch "${group}.peddy.ped"
		touch "${group}.sex_check.csv"
		touch "${group}_peddy.INFO"

		${peddy_version(task)}
		"""
}
def peddy_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    peddy: \$(echo \$(python -m peddy --version 2>&1) | sed 's/^.*peddy, version //')
	END_VERSIONS
	"""
}

// Extract all variants (
process fastgnomad {
	cpus 2
	memory '40 GB'
	tag "$group"
	publishDir "${params.results_output_dir}/vcf", mode: 'copy', overwrite: 'true', pattern: '*.vcf'
	time '2h'

	input:
		tuple val(group), path(vcf)

	output:
		tuple val(group), path("${group}.SNPs.vcf"), emit: vcf

	when:
		params.antype == "wgs"

	script:
		"""
		gzip -c $vcf > ${vcf}.gz
		annotate -g $params.FASTGNOMAD_REF -i ${vcf}.gz > ${group}.SNPs.vcf
		"""

	stub:
		"""
		touch "${group}.SNPs.vcf"
		"""
}


// Call UPD regions
process upd {
	tag "$group"
	time '1h'
	memory '1 GB'
	cpus 2

	input:
		tuple val(group), path(vcf)
		tuple val(group2), val(id), val(mother), val(father)

	output:
		path("upd.bed"), emit: upd_bed
		tuple val(group), path("upd.sites.bed"), emit: upd_sites
		path "*versions.yml", emit: versions

	script:
		if( params.mode == "family" && params.trio ) {
			"""
			upd --vcf $vcf --proband $id --mother $mother --father $father --af-tag GNOMADAF regions > upd.bed
			upd --vcf $vcf --proband $id --mother $mother --father $father --af-tag GNOMADAF sites > upd.sites.bed

			${upd_version(task)}
			"""
		}
		else {
			"""
			touch "upd.bed"
			touch "upd.sites.bed"

			${upd_version(task)}
			"""
		}

	stub:
		"""
		touch "upd.bed"
		touch "upd.sites.bed"

		${upd_version(task)}
		"""
}
def upd_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    upd: \$(echo \$(upd --version 2>&1))
	END_VERSIONS
	"""
}


process upd_table {
	publishDir "${params.results_output_dir}/plots", mode: 'copy' , overwrite: 'true'
	tag "$group"
	time '1h'
	memory '1 GB'
	cpus 2

	input:
		tuple val(group), path(upd_sites)

	output:
		path("${group}.UPDtable.xls")

	when:
		params.mode == "family" && params.trio

	script:
		"""
		upd_table.pl $upd_sites > ${group}.UPDtable.xls
		"""

	stub:
		"""
		touch "${group}.UPDtable.xls"
		"""
}


// Call ROH regions
process roh {
	tag "$group"
	time '1h'
	memory '1 GB'
	cpus 2

	input:
		tuple val(group), path(vcf)

	output:
		tuple val(group), path("roh.txt"), emit: roh_plot
		path "*versions.yml", emit: versions

	script:
		"""
		bcftools roh --rec-rate 1e-9 --AF-tag GNOMADAF ${vcf} -o roh.txt
		${roh_version(task)}
		"""

	stub:
		"""
		touch "roh.txt"
		${roh_version(task)}
		"""
}
def roh_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    bcftools: \$(echo \$(bcftools --version 2>&1) | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
	END_VERSIONS
	"""
}

// // Create coverage profile using GATK
// process gatkcov {
// 	publishDir "${params.results_output_dir}/cov", mode: 'copy' , overwrite: 'true', pattern: '*.tsv'
// 	tag "$group"
// 	cpus 2
// 	memory '80 GB'
// 	time '5h'

// 	input:
// 		tuple val(id), val(group), path(bam), path(bai), gr, sex, type

// 	output:
// 		tuple val(group), val(id), val(type), sex, path("${id}.standardizedCR.tsv"), path("${id}.denoisedCR.tsv"), emit: cov_plot, cov_gens
// 		path "*versions.yml", emit: versions

// 	when:
// 		params.gatkcov

// 	script:
// 		"""
// 		source activate gatk4-env

// 		gatk CollectReadCounts \\
// 			-I $bam -L $params.COV_INTERVAL_LIST \\
// 			--interval-merging-rule OVERLAPPING_ONLY -O ${bam}.hdf5

// 		gatk --java-options "-Xmx30g" DenoiseReadCounts \\
// 			-I ${bam}.hdf5 --count-panel-of-normals ${PON[sex]} \\
// 			--standardized-copy-ratios ${id}.standardizedCR.tsv \\
// 			--denoised-copy-ratios ${id}.denoisedCR.tsv

// 		gatk PlotDenoisedCopyRatios \\
// 			--standardized-copy-ratios ${id}.standardizedCR.tsv \\
// 			--denoised-copy-ratios ${id}.denoisedCR.tsv \\
// 			--sequence-dictionary $params.GENOMEDICT \\
// 			--minimum-contig-length 46709983 --output . --output-prefix $id

// 		${gatkcov_version(task)}
// 		"""

// 	stub:
// 		"""
// 		source activate gatk4-env
// 		touch "${id}.standardizedCR.tsv"
// 		touch "${id}.denoisedCR.tsv"

// 		${gatkcov_version(task)}
// 		"""
// }
// def gatkcov_version(task) {
// 	"""
// 	cat <<-END_VERSIONS > ${task.process}_versions.yml
// 	${task.process}:
// 	    gatk: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$// ; s/-SNAPSHOT//')
// 	END_VERSIONS
// 	"""
// }


// // Plot ROH, UPD and coverage in a genomic overview plot
// process overview_plot {

// 	cpus 2
// 	tag "$group"
// 	time '1h'
// 	memory '5 GB'
// 	publishDir "${params.results_output_dir}/plots", mode: 'copy' , overwrite: 'true', pattern: "*.png"

// 	input:
// 		path(upd)
// 		tuple val(group), path(roh)
// 		tuple val(group), val(id), val(type), sex, path(cov_stand), path(cov_denoised)

// 	output:
// 		path("${group}.genomic_overview.png")
// 		tuple val(group), path("${group}_oplot.INFO"), emit: oplot_INFO

// 	script:
// 		proband_idx = type.findIndexOf{ it == "proband" }
// 		"""
// 		genome_plotter.pl --dict $params.GENOMEDICT \\
// 			--sample ${id[proband_idx]} \\
// 			--upd $upd \\
// 			--roh $roh \\
// 			--sex ${sex[proband_idx]} \\
// 			--cov ${cov_denoised[proband_idx]} \\
// 			--out ${group}.genomic_overview.png
// 		echo "IMG overviewplot	${params.accessdir}/plots/${group}.genomic_overview.png" > ${group}_oplot.INFO
// 		"""

// 	stub:
// 		"""
// 		touch "${group}.genomic_overview.png"
// 		touch "${group}_oplot.INFO"
// 		"""
// }

// process generate_gens_data {
// 	publishDir "${params.results_output_dir}/plot_data", mode: 'copy' , overwrite: 'true', pattern: "*.gz*"
// 	publishDir "${params.crondir}/gens", mode: 'copy', overwrite: 'true', pattern: "*.gens"
// 	tag "$group"
// 	cpus 1
// 	time '3h'
// 	memory '5 GB'

// 	input:
// 		tuple val(id), val(group), path(gvcf), g, val(type), sex, path(cov_stand), path(cov_denoise)

// 	output:
// 		tuple path("${id}.cov.bed.gz"), path("${id}.baf.bed.gz"), path("${id}.cov.bed.gz.tbi"), path("${id}.baf.bed.gz.tbi"), path("${id}.overview.json.gz")
// 		path("${id}.gens"), emit: gens_middleman

// 	when:
// 		params.prepare_gens_data

// 	script:
// 		"""
// 		generate_gens_data.pl $cov_stand $gvcf $id $params.GENS_GNOMAD
// 		echo "gens load sample --sample-id $id --genome-build 38 --baf ${params.gens_accessdir}/${id}.baf.bed.gz --coverage ${params.gens_accessdir}/${id}.cov.bed.gz --overview-json ${params.gens_accessdir}/${id}.overview.json.gz" > ${id}.gens
// 		"""

// 	stub:
// 		"""
// 		touch "${id}.cov.bed.gz"
// 		touch "${id}.baf.bed.gz"
// 		touch "${id}.cov.bed.gz.tbi"
// 		touch "${id}.baf.bed.gz.tbi"
// 		touch "${id}.overview.json.gz"
// 		touch "${id}.gens"
// 		"""
// }

// SV-calling //

// GATK panel+wgs //

process gatk_coverage {
	cpus 2
	memory '50GB'
	time '2h'
	container  "${params.container_gatk}"
	tag "$id"
	input:
		tuple val(group), val(id), path(bam), path(bai)

	output:
		tuple val(group), val(id), path("${id}.tsv"), emit: coverage_tsv
		path "*versions.yml", emit: versions


	when:
		params.gatkcnv

	script:
		"""
		export THEANO_FLAGS="base_compiledir=."
		export MKL_NUM_THREADS=${task.cpus}
		export OMP_NUM_THREADS=${task.cpus}
		set +u
		source activate gatk
		gatk --java-options "-Xmx20g" CollectReadCounts \\
			-L $params.gatk_intervals \\
			-R $params.genome_file \\
			-imr OVERLAPPING_ONLY \\
			-I $bam \\
			--format TSV -O ${id}.tsv

		${gatk_coverage_version(task)}
		"""

	stub:
		"""
		export THEANO_FLAGS="base_compiledir=."
		export MKL_NUM_THREADS=${task.cpus}
		export OMP_NUM_THREADS=${task.cpus}
		set +u
		source activate gatk
		touch "${id}.tsv"

		${gatk_coverage_version(task)}
		"""
}
def gatk_coverage_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    gatk: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$// ; s/-SNAPSHOT//')
	END_VERSIONS
	"""
}

process gatk_call_ploidy {
	cpus 10
	memory '50GB'
	time '2h'
	container  "${params.container_gatk}"
	tag "$id"

	input:
		tuple val(group), val(id), path(coverage_tsv)

	output:
		tuple val(group), val(id), path("ploidy.tar"), emit: call_ploidy
		path "*versions.yml", emit: versions

	script:
		"""
		export THEANO_FLAGS="base_compiledir=."
		export MKL_NUM_THREADS=${task.cpus}
		export OMP_NUM_THREADS=${task.cpus}
		set +u
		source activate gatk
		gatk --java-options "-Xmx20g" DetermineGermlineContigPloidy \\
			--model $params.ploidymodel \\
			-I $coverage_tsv \\
			-O ploidy/ \\
			--output-prefix $group
		tar -cvf ploidy.tar ploidy/

		${gatk_call_ploidy_version(task)}
		"""

	stub:
		"""
		export THEANO_FLAGS="base_compiledir=."
		export MKL_NUM_THREADS=${task.cpus}
		export OMP_NUM_THREADS=${task.cpus}
		set +u
		source activate gatk
		touch "ploidy.tar"

		${gatk_call_ploidy_version(task)}
		"""
}
def gatk_call_ploidy_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    gatk: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$// ; s/-SNAPSHOT//')
	END_VERSIONS
	"""
}

process gatk_call_cnv {
	cpus 8
	memory '50GB'
	time '3h'
	container  "${params.container_gatk}"
	tag "$id"

	input:
	tuple val(group), val(id), path(tsv), path(ploidy), val(i), val(refpart) //TODO: is reffart a path or val


	output:
	//TODO: wtf is i
	tuple val(group), val(id), val(i), path("${group}_${i}.tar"), emit: gatk_calls
	path "*versions.yml", emit: versions

	script:
		"""
		export THEANO_FLAGS="base_compiledir=."
		set +u
		source activate gatk
		export HOME=/local/scratch
		export MKL_NUM_THREADS=${task.cpus}
		export OMP_NUM_THREADS=${task.cpus}
		tar -xvf ploidy.tar
		mkdir ${group}_${i}
		gatk --java-options "-Xmx25g" GermlineCNVCaller \\
			--run-mode CASE \\
			-I $tsv \\
			--contig-ploidy-calls ploidy/${group}-calls/ \\
			--model ${refpart} \\
			--output ${group}_${i}/ \\
			--output-prefix ${group}_${i}
		tar -cvf ${group}_${i}.tar ${group}_${i}/

		${gatk_call_cnv_version(task)}
		"""

	stub:
		"""
		export THEANO_FLAGS="base_compiledir=."
		set +u
		source activate gatk
		export HOME=/local/scratch
		export MKL_NUM_THREADS=${task.cpus}
		export OMP_NUM_THREADS=${task.cpus}
		source activate gatk
		touch "${group}_${i}.tar"

		${gatk_call_cnv_version(task)}
		"""
}
def gatk_call_cnv_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    gatk: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$// ; s/-SNAPSHOT//')
	END_VERSIONS
	"""
}

process postprocessgatk {
	cpus 5
	memory '50GB'
	time '3h'
	container  "${params.container_gatk}"
	publishDir "${params.results_output_dir}/sv_vcf/", mode: 'copy', overwrite: 'true', pattern: '*.vcf.gz'
	tag "$id"

	input:
	// TODO: wtf is i
	tuple val(group), val(id), val(i), path(tar), path(ploidy), val(shard_no), val(shard)

	output:
		tuple val(group), val(id),path("genotyped-intervals-${group}-vs-cohort30.vcf.gz"), path("genotyped-segments-${group}-vs-cohort30.vcf.gz"), path("denoised-${group}-vs-cohort30.vcf.gz"), emit: called_gatk
		path "*versions.yml", emit: versions


	script:
		def modelshards = shard.join(' --model-shard-path ') // join each reference shard
		def caseshards = []
		// TODO: put into func
		// // join each shard(n) that's been called
	    i.each { shard_name ->
        	def shard_path = group + '_' + shard_name + '/' + group + '_' + shard_name + '-calls'
        	caseshards << shard_path
    	}
		caseshards = caseshards.join( ' --calls-shard-path ')
		version_str = postprocessgatk_version(task)
		"""
		THEANO_FLAGS="base_compiledir=/fs1/resources/theano"

		for model in ${tar}; do
			tar -xvf \$model
		done

		tar -xvf ${ploidy}

		set +u
		source activate gatk
		export MKL_NUM_THREADS=${task.cpus}
		export OMP_NUM_THREADS=${task.cpus}

		gatk --java-options "-Xmx25g" PostprocessGermlineCNVCalls \
			--allosomal-contig X --allosomal-contig Y \
			--contig-ploidy-calls ploidy/${group}-calls/ \
			--sample-index 0 \\
			--output-genotyped-intervals genotyped-intervals-${group}-vs-cohort30.vcf.gz \
			--output-genotyped-segments genotyped-segments-${group}-vs-cohort30.vcf.gz \
			--output-denoised-copy-ratios denoised-${group}-vs-cohort30.vcf.gz \
			--sequence-dictionary ${params.GENOMEDICT} \
			--calls-shard-path ${caseshards} \
			--model-shard-path ${modelshards}

		echo "${version_str}" > "${task.process}_versions.yml"
		"""

	stub:
		def modelshards = shard.join(' --model-shard-path ') // join each reference shard
		def caseshards = []
		// TODO: lsp complains about indexing var
		// // join each shard(n) that's been called
	    i.each { shard_name ->
        	def shard_path = group + '_' + shard_name + '/' + group + '_' + shard_name + '-calls'
        	caseshards << shard_path
    	}
		caseshards = caseshards.join( ' --calls-shard-path ')
		version_str = postprocessgatk_version(task)
		"""
		THEANO_FLAGS="base_compiledir=/fs1/resources/theano"
		set +u
		source activate gatk
		export MKL_NUM_THREADS=${task.cpus}
		export OMP_NUM_THREADS=${task.cpus}
		source activate gatk
		touch "genotyped-intervals-${group}-vs-cohort30.vcf.gz"
		touch "genotyped-segments-${group}-vs-cohort30.vcf.gz"
		touch "denoised-${group}-vs-cohort30.vcf.gz"

		echo "${modelshards}"
		echo "${caseshards}"

		echo "${version_str}" > "${task.process}_versions.yml"
		"""
}
def postprocessgatk_version(task) {
	// This docstring looks different
	// If spaces similarly to the others, this leads to additional whitespace above and below the version text
	"""${task.process}:
	    gatk: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$// ; s/-SNAPSHOT//')"""
}


process filter_merge_gatk {
	cpus 2
	tag "$group"
	time '2h'
	memory '1 GB'
	publishDir "${params.results_output_dir}/sv_vcf", mode: 'copy', overwrite: 'true'

	input:
		tuple val(group), val(id), path(inter), path(gatk), path(denoised)

	output:
		tuple val(group), val(id), path("${id}.gatk.filtered.merged.vcf"), emit: merged_gatk

	script:
		"""
		filter_gatk.pl $gatk > ${id}.gatk.filtered.vcf
		mergeGATK.pl ${id}.gatk.filtered.vcf > ${id}.gatk.filtered.merged.vcf
		"""

	stub:
		"""
		touch "${id}.gatk.filtered.merged.vcf"
		"""
}

process manta {
	cpus  56
	publishDir "${params.results_output_dir}/sv_vcf/", mode: 'copy', overwrite: 'true', pattern: '*.vcf.gz'
	tag "$id"
	time '10h'
	memory '150 GB'
	input:
		tuple val(group), val(id), path(bam), path(bai)

	output:
		tuple val(group), val(id), path("${id}.manta.vcf.gz"), emit: vcf
		path "*versions.yml", emit: versions

	script:
		bams = bam.join('--bam ')

		"""
		configManta.py --bam $bam --reference ${params.genome_file} --runDir .
		python runWorkflow.py -m local -j ${task.cpus}
		mv results/variants/diploidSV.vcf.gz ${id}.manta.vcf.gz
		mv results/variants/diploidSV.vcf.gz.tbi ${id}.manta.vcf.gz.tbi

		${manta_version(task)}
		"""

	stub:
		"""
		touch "${id}.manta.vcf.gz"
		${manta_version(task)}
		"""
}
def manta_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    manta: \$( configManta.py --version )
	END_VERSIONS
	"""
}

process manta_panel {
	cpus  20
	publishDir "${params.results_output_dir}/sv_vcf/", mode: 'copy', overwrite: 'true', pattern: '*.vcf.gz'
	tag "$id"
	time '1h'
	memory '50 GB'


	input:
		tuple val(group), val(id), path(bam), path(bai)

	output:
		tuple val(group), val(id), path("${id}.manta.vcf.gz"), emit: vcf
		path "*versions.yml", emit: versions

	when:
		params.sv && params.antype == "panel"

	script:
		"""
		configManta.py --bam $bam --reference ${params.genome_file} --runDir . --exome --callRegions $params.bedgz --generateEvidenceBam
		python runWorkflow.py -m local -j ${task.cpus}
		mv results/variants/diploidSV.vcf.gz ${id}.manta.vcf.gz
		mv results/variants/diploidSV.vcf.gz.tbi ${id}.manta.vcf.gz.tbi

		${manta_panel_version(task)}
		"""

	stub:
		"""
		touch "${id}.manta.vcf.gz"
		${manta_panel_version(task)}
		"""
}
def manta_panel_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	    manta: \$( configManta.py --version )
	END_VERSIONS
	"""
}


// process cnvkit_panel {
// 	cpus  5
// 	container  "${params.container_twist_myeloid}"
// 	publishDir "${params.results_output_dir}/sv_vcf/", mode: 'copy', overwrite: 'true', pattern: '*.vcf'
// 	publishDir "${params.results_output_dir}/plots/", mode: 'copy', overwrite: 'true', pattern: '*.png'
// 	tag "$id"
// 	time '1h'
// 	memory '20 GB'
// 	input:
// 		tuple val(group), val(id), path(bam), path(bai), path(vcf), path(multi), val(INS_SIZE), val(MEAN_DEPTH), val(COV_DEV)
// 		//tuple id, val(INS_SIZE), val(MEAN_DEPTH), val(COV_DEV)
// 		//tuple val(group), val(id), path(vcf)

// 	output:
// 		tuple val(group), val(id), path("${id}.cnvkit_filtered.vcf"), emit: called_cnvkit_panel
// 		path("${id}.call.cns"), emit: unfiltered_cns
// 		path("${group}.genomic_overview.png")
// 		tuple val(group), path("${group}_oplot.INFO"), emit: cnvkit_INFO
// 		path "*versions.yml", emit: versions


// 	when:
// 		params.sv && params.antype == "panel"

// 	script:
// 		"""
// 		cnvkit.py batch $bam -r $params.cnvkit_reference -p 5 -d results/
// 		cnvkit.py call results/*.cns -v $vcf -o ${id}.call.cns
// 		filter_cnvkit.pl ${id}.call.cns $MEAN_DEPTH > ${id}.filtered
// 		cnvkit.py export vcf ${id}.filtered -i "$id" > ${id}.cnvkit_filtered.vcf
// 		cnvkit.py scatter -s results/*dedup.cn{s,r} -o ${group}.genomic_overview.png -v $vcf -i $id
// 		echo "IMG overviewplot	${params.accessdir}/plots/${group}.genomic_overview.png" > ${group}_oplot.INFO

// 		${cnvkit_panel_version(task)}
// 		"""

// 	stub:
// 		"""
// 		touch "${id}.cnvkit_filtered.vcf"
// 		touch "${id}.call.cns"
// 		touch "${group}.genomic_overview.png"
// 		touch "${group}_oplot.INFO"

// 		${cnvkit_panel_version(task)}
// 		"""
// }
// def cnvkit_panel_version(task) {
// 	"""
// 	cat <<-END_VERSIONS > ${task.process}_versions.yml
// 	${task.process}:
// 	    cnvkit: \$(cnvkit.py version | sed -e 's/cnvkit v//g')
// 	END_VERSIONS
// 	"""
// }

// process svdb_merge_panel {
// 	container  "${params.container_svdb}"
// 	cpus 2
// 	cache 'deep'
// 	tag "$group"
// 	publishDir "${params.results_output_dir}/sv_vcf/merged/", mode: 'copy', overwrite: 'true', pattern: '*.vcf'
// 	time '1h'
// 	memory '1 GB'
// 	input:
// 		tuple val(group), val(id), path(vcfs)

// 	output:
// 		tuple val(group), val(id), path("${group}.merged.vcf"), emit: ch_postprocess_merged_panel_sv
// 		path "*versions.yml", emit: versions


// 	when:
// 		params.antype == "panel"

// 	script:
// 		if (vcfs.size() > 1) {
// 			// for each sv-caller add idx, find vcf and find priority, add in priority order! //
// 			// index of vcfs added
// 			manta_idx = vcfs.findIndexOf{ it =~ 'manta' }
// 			cnvkit_idx = vcfs.findIndexOf{ it =~ 'cnvkit' }
// 			gatk_idx = vcfs.findIndexOf{ it =~ 'gatk' }

// 			// find vcfs //
// 			manta = manta_idx >= 0 ? vcfs[manta_idx].collect {it + ':manta ' } : null
// 			cnvkit = cnvkit_idx >= 0 ? vcfs[cnvkit_idx].collect {it + ':cnvkit ' } : null
// 			gatk = gatk_idx >= 0 ? vcfs[gatk_idx].collect {it + ':gatk ' } : null
// 			tmp = manta + gatk + cnvkit
// 			tmp = tmp - null
// 			vcfs_svdb = tmp.join(' ')

// 			// find priorities //
// 			mantap = manta_idx >= 0 ? 'manta' : null
// 			gatkp = gatk_idx >= 0 ? 'gatk' : null
// 			cnvkitp = cnvkit_idx >= 0 ? 'cnvkit' : null
// 			tmpp = [mantap, gatkp, cnvkitp]
// 			tmpp = tmpp - null
// 			priority = tmpp.join(',')

// 			"""
// 			svdb \\
// 			  --merge \\
// 			  --vcf $vcfs_svdb \\
// 			  --no_intra \\
// 			  --pass_only \\
// 			  --bnd_distance 2500 \\
// 			  --overlap 0.7 \\
// 			  --priority $priority \\
// 			  --ins_distance 0 > ${group}.merged.tmp


// 			# copy callers out of INFO.tuple to INFO.SCOUT_CUSTOM
// 			add_callers_to_scout_custom.py \\
// 				--callers $priority \\
// 				--merged_vcf ${group}.merged.tmp > ${group}.merged.callers.tmp

// 			add_vcf_header_info_records.py \\
// 				--vcf ${group}.merged.callers.tmp \\
// 				--info SCOUT_CUSTOM . String "Custom annotations for scout" '' '' \\
// 				--output ${group}.merged.vcf

// 			${svdb_merge_panel_version(task)}
// 			"""
// 		}
// 		else {
// 			"""
// 			mv $vcf ${group}.merged.vcf
// 			${svdb_merge_panel_version(task)}
// 			"""
// 		}

// 	stub:
// 		"""
// 		touch "${group}.merged.vcf"
// 		${svdb_merge_panel_version(task)}
// 		"""
// }
// def svdb_merge_panel_version(task) {
// 	"""
// 	cat <<-END_VERSIONS > ${task.process}_versions.yml
// 	${task.process}:
// 	    svdb: \$( echo \$(svdb) | head -1 | sed 's/usage: SVDB-\\([0-9]\\.[0-9]\\.[0-9]\\).*/\\1/' )
// 	END_VERSIONS
// 	"""
// }

// process postprocess_merged_panel_sv_vcf {
// 	cpus 2
// 	tag "$group"
// 	publishDir "${params.results_output_dir}/sv_vcf/merged/", mode: 'copy', overwrite: 'true', pattern: '*.vcf'
// 	time '1h'
// 	memory '1 GB'


// 	input:
// 		tuple val(group), val(id), path(merged_vcf)
// 		tuple val(group), val(id), path(melt_vcf)

// 	output:
// 		tuple val(group), val(id), path("${group}.merged.bndless.genotypefix.melt.vcf"), emit: vep_sv_panel, annotsv_panel
// 		tuple val(group), path("${group}.merged.bndless.genotypefix.melt.vcf"), emit: loqusdb_sv_panel
// 		path "*versions.yml", emit: versions


// 	script:
// 		"""
// 		# Remove BNDs
// 		grep -v "BND" $merged_vcf > ${group}.merged.bndless.vcf

// 		# Any 0/0 GT -> 0/1, otherwise loqus will reject them.
// 		modify_cnv_genotypes_for_loqusdb.pl --merged_panel_sv_vcf ${group}.merged.bndless.vcf > ${group}.merged.bndless.genotypefix.vcf

// 		# Add MELT data to info vars:
// 		add_vcf_header_info_records.py \\
// 			--vcf ${group}.merged.bndless.genotypefix.vcf \\
// 			--info MELT_RANK . String "Evidence level 1-5, 5 - highest" '' '' \\
// 			--info MELT_QC . String "Quality of call" '' '' \\
// 			--output ${group}.merged.bndless.genotypefix.headers.vcf

// 		# Combine with MELT:
// 		vcf-concat  ${group}.merged.bndless.genotypefix.headers.vcf $melt_vcf | vcf-sort -c > ${group}.merged.bndless.genotypefix.melt.vcf
// 		${postprocess_merged_panel_sv_version(task)}
// 		"""

// 	stub:
// 		"""
// 		touch ${group}.merged.bndless.genotypefix.melt.vcf
// 		${postprocess_merged_panel_sv_version(task)}
// 		"""

// }
// def postprocess_merged_panel_sv_version(task) {
// 	"""
// 	cat <<-END_VERSIONS > ${task.process}_versions.yml
// 	${task.process}:
// 	    vcftools: \$(vcftools --version | cut -f 2 -d " " | tr -d "()")
// 	END_VERSIONS
// 	"""
// }

// process tiddit {
// 	cpus  2
// 	publishDir "${params.results_output_dir}/sv_vcf/", mode: 'copy', overwrite: 'true', pattern: '*.vcf'
// 	time '10h'
// 	tag "$id"
// 	memory '15 GB'

// 	input:
// 		tuple val(group), val(id), path(bam), path(bai)

// 	output:
// 		tuple val(group), val(id), path("${id}.tiddit.filtered.vcf"), emit: called_tiddit
// 		path "*versions.yml", emit: versions


// 	when:
// 		params.sv && params.antype == "wgs"

// 	script:
// 		"""
// 		TIDDIT.py --sv -o ${id}.tiddit --bam $bam
// 		grep -E \"#|PASS\" ${id}.tiddit.vcf > ${id}.tiddit.filtered.vcf
// 		${tiddit_version(task)}
// 		"""

// 	stub:
// 		"""
// 		touch "${id}.tiddit.filtered.vcf"
// 		${tiddit_version(task)}
// 		"""
// }
// def tiddit_version(task) {
// 	"""
// 	cat <<-END_VERSIONS > ${task.process}_versions.yml
// 	${task.process}:
// 	    tiddit: \$(echo \$(TIDDIT.py 2>&1) | sed 's/^.*TIDDIT-//; s/ .*\$//')
// 	END_VERSIONS
// 	"""
// }

// process svdb_merge {
// 	cpus 2
// 	container  "${params.container_svdb}"
// 	tag "$group"
// 	publishDir "${params.results_output_dir}/sv_vcf/merged/", mode: 'copy', overwrite: 'true', pattern: '*.vcf'
// 	time '2h'
// 	memory '1 GB'

// 	input:
// 		tuple val(group), val(id), path(mantaV)
// 		tuple val(group), val(id), path(tidditV)
// 		tuple val(group), val(id), path(gatkV)

// 	output:
// 		tuple val(group), val(id), path("${group}.merged.bndless.vcf"), emit: vep_sv, annotsv_vcf
// 		tuple val(group), path("${group}.merged.vcf"), emit: loqusdb_sv
// 		path "*versions.yml", emit: versions

// 	script:
// 		if (params.mode == "family") {
// 			vcfs = []
// 			manta = []
// 			tiddit = []
// 			gatk = []

// 			/*
// 			 Order in which VCFs are merged matters when the merged SV
// 			 is annotated with final position/length, which affects
// 			 artefact matching in loqusdb.

// 			 A possibly better way to sort here would be to sort the
// 			 file by familial-relation (e.g. always sort proband-mother-father)
// 			 this would ensure the same merge-order regardless of sample-id
// 			 */

// 			mantaV = mantaV.collect { it.toString() }.sort()
// 			gatkV = gatkV.collect { it.toString() }.sort()
// 			tidditV = tidditV.collect { it.toString() }.sort()

// 		//TODO: lsp complains about for loop?
// 			// for (i = 1; i <= mantaV.size(); i++) {
// 			// 	tmp = mantaV[i-1] + ':manta' + "${i}"
// 			// 	tmp1 = tidditV[i-1] + ':tiddit' + "${i}"
// 			// 	tmp2 = gatkV[i-1] + ':gatk' + "${i}"
// 			// 	vcfs = vcfs + tmp + tmp1 + tmp2
// 			// 	mt = 'manta' + "${i}"
// 			// 	tt = 'tiddit' + "${i}"
// 			// 	ct = 'gatk' + "${i}"
// 			// 	manta = manta + mt
// 			// 	tiddit = tiddit + tt
// 			// 	gatk = gatk + ct
// 			// }

// 			prio = manta + tiddit + gatk
// 			prio = prio.join(',')
// 			vcfs = vcfs.join(' ')
// 			"""
// 			svdb \\
// 				--merge \\
// 				--vcf $vcfs \\
// 				--no_intra \\
// 				--pass_only \\
// 				--bnd_distance 2500 \\
// 				--overlap 0.7 \\
// 				--priority $prio \\
// 				--ins_distance 0 > ${group}.merged.vcf

// 			grep -v BND ${group}.merged.vcf > ${group}.merged.bndless.vcf

// 			${svdb_merge_version(task)}
// 			"""
// 		}

// 		else {
// 			tmp = mantaV.collect {it + ':manta ' } + tidditV.collect {it + ':tiddit ' } + gatkV.collect {it + ':gatk ' }
// 			vcfs = tmp.join(' ')
// 			"""
// 			svdb \\
// 				--merge \\
// 				--vcf $vcfs \\
// 				--no_intra \\
// 				--pass_only \\
// 				--bnd_distance 2500 \\
// 				--overlap 0.7 \\
// 				--priority manta,tiddit,gatk \\
// 				--ins_distance 0 > ${group}.merged.vcf

// 			grep -v BND ${group}.merged.vcf > ${group}.merged.bndless.vcf

// 			${svdb_merge_version(task)}
// 			"""
// 		}

// 	stub:
// 		"""
// 		touch "${group}.merged.vcf"
// 		touch "${group}.merged.bndless.vcf"

// 		${svdb_merge_version(task)}
// 		"""
// }
// def svdb_merge_version(task) {
// 	"""
// 	cat <<-END_VERSIONS > ${task.process}_versions.yml
// 	${task.process}:
// 	    svdb: \$( echo \$(svdb) | head -1 | sed 's/usage: SVDB-\\([0-9]\\.[0-9]\\.[0-9]\\).*/\\1/' )
// 	END_VERSIONS
// 	"""
// }

// process dummy_svvcf_for_loqusdb {

// 	// add_to_loqusb won't run if no svvcf is generated
// 	// this process creates dummy svvcf for no-SV runs
// 	// assay input only exists to disable nextflow warning
// 	// for channels emitting with less than two input elements

// 	cpus 1
// 	tag "$group"
// 	memory '10 MB'
// 	time '10m'

// 	input:
// 		tuple val(group), val(assay)

// 	output:
// 		tuple val(group), path("${group}.dummy.sv.vcf"), emit: dummy_svvcf_ch

// 	when:
// 		!params.sv

// 	script:
// 		"""
// 		touch ${group}.dummy.sv.vcf
// 		"""

// 	stub:
// 		"""
// 		touch ${group}.dummy.sv.vcf
// 		"""
// }

// process add_to_loqusdb {
// 	cpus 1
// 	publishDir "${params.crondir}/loqus", mode: 'copy' , overwrite: 'true'
// 	tag "$group"
// 	memory '100 MB'
// 	time '25m'
// 	input:
// 		tuple val(group), val(type), path(vcf), path(tbi), val(type), path(ped)
// 		tuple val(group), path(svvcf)

// 	output:
// 		path("${group}*.loqus"), emit: loqusdb_done


// 	when:
// 		!params.noupload && !params.reanalyze

// 	script:
// 		"""
// 		sv_variants=""
// 		nbr_svvcf_records=\$(grep -v '^#' ${svvcf} | wc -l)

// 		if (( \$nbr_svvcf_records > 0 )); then
// 			sv_variants="--sv-variants ${params.accessdir}/sv_vcf/merged/${svvcf}"
// 		fi

// 		echo "-db $params.loqusdb load -f ${params.accessdir}/ped/${ped} --variant-file ${params.accessdir}/vcf/${vcf} \$sv_variants" > ${group}.loqus
// 		"""

// 	stub:
// 		"""
// 		touch "${group}.loqus"
// 		"""
// }

// process annotsv {
// 	container  "${params.container_annotsv}"
// 	cpus 2
// 	tag "$group"
// 	publishDir "${params.results_output_dir}/annotsv/", mode: 'copy', overwrite: 'true', pattern: '*.tsv'
// 	time '5h'
// 	memory '20 GB'

// 	input:
// 		tuple val(group), val(id), path(sv)

// 	output:
// 		tuple val(group), path("${group}_annotsv.tsv"), emit: annotsv, annotsv_ma, annotsv_fa
// 		path "*versions.yml", emit: versions

// 	shell:
// 		version_str = annotsv_version(task)
// 		'''
// 		export ANNOTSV="/AnnotSV"
// 		/AnnotSV/bin/AnnotSV -SvinputFile !{sv} \\
// 			-typeOfAnnotation full \\
// 			-outputDir !{group} \\
// 			-genomeBuild GRCh38
// 		if [ -f !{group}/*.annotated.tsv ]; then
// 			mv !{group}/*.annotated.tsv !{group}_annotsv.tsv
// 		else
// 		    echo "1\n" > !{group}_annotsv.tsv
// 		fi
// 		echo "!{version_str}" > "!{task.process}_versions.yml"
// 		'''

// 	stub:
// 		version_str = annotsv_version(task)
// 		"""
// 		export ANNOTSV="/AnnotSV"
// 		touch "${group}_annotsv.tsv"

// 		echo "${version_str}" > "${task.process}_versions.yml"
// 		"""
// }
// def annotsv_version(task) {
// 	"""${task.process}:
// 	    annotsv: \$( echo \$(/AnnotSV/bin/AnnotSV --version) | sed -e "s/AnnotSV //g ; s/Copyright.*//" )"""
// }

// process vep_sv {
// 	cpus 10
// 	container  "${params.container_vep}"
// 	tag "$group"
// 	memory '50 GB'
// 	time '1h'

// 	input:
// 		tuple val(group), val(id), path(vcf)

// 	output:
// 		tuple val(group), val(id), path("${group}.vep.vcf"), emit: vep_sv_vcf
// 		path "*versions.yml", emit: versions

// 	script:
// 		"""

// 		# Temporary fix for VEP 111.0 annotation bug, where certain MANTA indels are being skipped by VEP
// 		# See: https://github.com/Ensembl/ensembl-vep/issues/1631#issuecomment-1985973568
// 		# Edit 2024-11-01: Fixed in VEP 112.0

// 		sed 's/SVTYPE=/BAZBAZ=/' $vcf > ${group}.vep111-workaround.vcf

// 		vep \\
// 			-i ${group}.vep111-workaround.vcf \\
// 			-o ${group}.vep.vcf \\
// 			--offline \\
// 			--merged \\
// 			--sift b --polyphen b --ccds --hgvs --symbol --numbers --domains --regulatory --canonical --protein --biotype --af --af_1kg --max_af --pubmed --uniprot --mane --tsl --appris --variant_class --gene_phenotype --mirna \\
// 			--synonyms $params.VEP_SYNONYMS \\
// 			--vcf \\
// 			--no_stats \\
// 			--fork ${task.cpus} \\
// 			--force_overwrite \\
// 			--plugin LoFtool \\
// 			--fasta $params.VEP_FASTA \\
// 			--dir_cache $params.VEP_CACHE \\
// 			--dir_plugins $params.VEP_PLUGINS \\
// 			--max_sv_size $params.VEP_MAX_SV_SIZE \\
// 			--distance $params.VEP_TRANSCRIPT_DISTANCE \\
// 			-cache \\
// 			--format vcf

// 		# Re-enable SVTYPE:
// 		sed -i 's/BAZBAZ=/SVTYPE=/' ${group}.vep.vcf
// 		${vep_sv_version(task)}
// 		"""

// 	stub:
// 		"""
// 		touch "${group}.vep.vcf"
// 		${vep_sv_version(task)}
// 		"""
// }
// def vep_sv_version(task) {
// 	"""
// 	cat <<-END_VERSIONS > ${task.process}_versions.yml
// 	${task.process}:
// 	    vep: \$( echo \$(vep --help 2>&1) | sed 's/^.*Versions:.*ensembl-vep : // ; s/ .*\$//')
// 	END_VERSIONS
// 	"""
// }

// process postprocess_vep_sv {
// 	cpus  2
// 	memory '10GB'
// 	time '1h'
// 	tag "$group"
// 	container  "${params.container_svdb}"

// 	input:
// 		tuple val(group), val(id), path(vcf)

// 	output:
// 		tuple val(group), path("${group}.vep.clean.merge.vcf"), emit: add_omim_vcf
// 		path "*versions.yml", emit: versions

// 	script:
// 		"""
// 		# Filter variants with FILTER != . or PASS and variants missing CSQ field.
// 		postprocess_vep_vcf.py $vcf > ${group}.vep.clean.vcf
// 		svdb --merge --overlap 0.9 --notag --vcf ${group}.vep.clean.vcf --ins_distance 0 > ${group}.vep.clean.merge.tmp.vcf

// 		# --notag above will remove set
// 		add_vcf_header_info_records.py \\
// 			--vcf ${group}.vep.clean.merge.tmp.vcf \\
// 			--info tuple 1 String "Source VCF for the merged record in SVDB" '' '' \\
// 			--info VARID 1 String "The variant ID of merged samples" '' '' \\
// 			--output ${group}.vep.clean.merge.headers.tmp.vcf

// 		# Prepare annotations for scout:
// 		normalize_caller_names_in_svdb_fields.py ${group}.vep.clean.merge.headers.tmp.vcf --callers manta gatk tiddit > ${group}.vep.clean.merge.vcf
// 		${postprocess_vep_sv_version(task)}
// 		"""
// 	stub:
// 		"""
// 		touch "${group}.vep.clean.merge.vcf"
// 		${postprocess_vep_sv_version(task)}
// 		"""
// }
// def postprocess_vep_sv_version(task) {
// 	"""
// 	cat <<-END_VERSIONS > ${task.process}_versions.yml
// 	${task.process}:
// 	    svdb: \$( echo \$(svdb) | head -1 | sed 's/usage: SVDB-\\([0-9]\\.[0-9]\\.[0-9]\\).*/\\1/' )
// 	END_VERSIONS
// 	"""
// }

// process add_omim {
// 	cpus 2
// 	memory '10GB'
// 	time '1h'
// 	tag "$group"

// 	input:
// 		tuple val(group), path(vcf)

// 	output:
// 		tuple val(group), path("${group}.vep.clean.merge.omim.vcf"), emit: artefact_vcf

// 	script:
// 		"""
// 		add_omim.pl $params.OMIM_GENES < $vcf > ${group}.vep.clean.merge.omim.vcf
// 		"""
// 	stub:
// 		"""
// 		touch "${group}.vep.clean.merge.omim.vcf"
// 		"""
// }



// // Query artefact db
// process artefact {
// 	cpus 2
// 	tag "$group"
// 	time '10h'
// 	memory '10 GB'
// 	container  "${params.container_svdb}"


// 	input:
// 		tuple val(group), path(sv)
// 	output:
// 		tuple val(group), path("${group}.artefact.vcf"), emit: manip_vcf,manip_vcf_ma,manip_vcf_fa
// 		path "*versions.yml", emit: versions

// 	script:
// 		// use loqusdb dump not svdb database //
// 		if (params.gatkcnv) {
// 			"""
// 			svdb \\
// 			--query --bnd_distance 25000 --overlap 0.7 --in_occ Obs --out_occ ACOUNT --in_frq Frq --out_frq AFRQ  \\
// 			--db $params.svdb \\
// 			--ins_distance 0 \\
// 			--query_vcf $sv > ${group}.artefact.vcf

// 			${artefact_version(task)}
// 			"""
// 		}
// 		// for oncov1-0 still use svdb database remove in future//
// 		else {
// 			"""
// 			svdb \\
// 			--sqdb $params.svdb \\
// 			--query \\
// 			--query_vcf $sv \\
// 			--out_occ ACOUNT \\
// 			--ins_distance 0 \\
// 			--out_frq AFRQ > ${group}.artefact.vcf

// 			${artefact_version(task)}
// 			"""
// 		}

// 	stub:
// 		"""
// 		touch "${group}.artefact.vcf"
// 		${artefact_version(task)}
// 		"""
// }
// def artefact_version(task) {
// 	"""
// 	cat <<-END_VERSIONS > ${task.process}_versions.yml
// 	${task.process}:
// 	    svdb: \$( echo \$(svdb) | head -1 | sed 's/usage: SVDB-\\([0-9]\\.[0-9]\\.[0-9]\\).*/\\1/' )
// 	END_VERSIONS
// 	"""
// }


// process prescore {
// 	cpus 2
// 	tag "$group"
// 	memory '10 GB'
// 	time '1h'

// 	input:
// 		tuple val(group), path(sv_artefact), val(type), path(ped), path(annotsv)

// 	output:
// 		tuple val(group), val(type), path("${group}.annotatedSV.vcf"), emit: annotatedSV

// 	script:
// 		"""
// 		prescore_sv.pl \\
// 		--sv $sv_artefact --ped $ped --annotsv $annotsv --osv ${group}.annotatedSV.vcf
// 		"""

// 	stub:
// 		"""
// 		touch "${group}.annotatedSV.vcf"
// 		"""
// }

// process score_sv {
// 	tag "$group $mode"
// 	cpus 2
// 	memory '10 GB'
// 	time '2h'
// 	container  "${params.container_genmod}"

// 	input:
// 		tuple val(group), val(type), path(in_vcf)

// 	output:
// 		tuple val(group), val(type), path("${group_score}.sv.scored.vcf"), emit: ch_scored_sv
// 		path "*versions.yml", emit: versions

// 	script:
// 		def model = (params.mode == "family" && params.antype == "wgs") ? params.svrank_model : params.svrank_model_s
// 		def group_score = ( type == "ma" || type == "fa" ) ? "${group}_${type}" : group
// 		"""
// 		genmod score --family_id ${group_score} --score_config ${model} --rank_results --outfile "${group_score}.sv.scored.vcf" ${in_vcf}

// 		${score_sv_version(task)}
// 		"""

// 	stub:
// 		group_score = group
// 		"""
// 		touch "${group_score}.sv.scored.vcf"

// 		${score_sv_version(task)}
// 		"""
// }
// def score_sv_version(task) {
// 	"""
// 	cat <<-END_VERSIONS > ${task.process}_versions.yml
// 	${task.process}:
// 	    genmod: \$(echo \$(genmod --version 2>&1) | sed -e "s/^.*genmod version: //")
// 	END_VERSIONS
// 	"""
// }

// process bgzip_scored_genmod {
// 	tag "$group"
// 	cpus 4
// 	memory '1 GB'
// 	time '5m'
// 	publishDir "${params.results_output_dir}/vcf", mode: 'copy', overwrite: 'true', pattern: '*.vcf.gz'
// 	publishDir "${params.results_output_dir}/vcf", mode: 'copy', overwrite: 'true', pattern: '*.vcf.gz.tbi'
// 	container  "${params.container_bcftools}"

// 	input:
// 		tuple val(group), val(type), path(scored_sv_vcf)

// 	output:
// 		tuple val(group), val(type), path("${group_score}.sv.scored.sorted.vcf.gz"), path("${group_score}.sv.scored.sorted.vcf.gz.tbi"), emit: sv_rescore, sv_rescore_ma, sv_rescore_fa
// 		tuple val(group), path("${group_score}.sv.scored.sorted.vcf.gz"), emit: svvcf_bed, svvcf_pod
// 		tuple val(group), path("${group}_sv.INFO"), emit: sv_INFO
// 		path "*versions.yml", emit: versions

// 	script:
// 		group_score = ( type == "ma" || type == "fa" ) ? "${group}_${type}" : group
// 		"""
// 			bcftools sort -O v -o ${group_score}.sv.scored.sorted.vcf ${scored_sv_vcf}
// 			bgzip -@ ${task.cpus} ${group_score}.sv.scored.sorted.vcf -f
// 			tabix ${group_score}.sv.scored.sorted.vcf.gz -f
// 			echo "SV\t$type\t${params.accessdir}/vcf/${group_score}.sv.scored.sorted.vcf.gz" > ${group}_sv.INFO

// 			${bgzip_score_sv_version(task)}
// 		"""
// 	stub:
// 		group_score = ( type == "ma" || type == "fa" ) ? "${group}_${type}" : group
// 		"""
// 			touch "${group_score}.sv.scored.sorted.vcf.gz"
// 			touch "${group_score}.sv.scored.sorted.vcf.gz.tbi"
// 			touch "${group}_sv.INFO"

// 			${bgzip_score_sv_version(task)}
// 		"""
// }
// def bgzip_score_sv_version(task) {
// 	"""
// 	cat <<-END_VERSIONS > ${task.process}_versions.yml
// 	${task.process}:
// 	    tabix: \$(echo \$(tabix --version 2>&1) | sed 's/^.*(htslib) // ; s/ Copyright.*//')
// 	    bcftools: \$(echo \$(bcftools --version 2>&1) | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
// 	END_VERSIONS
// 	"""
// }

// process compound_finder {
// 	cpus 2
// 	tag "$group ${params.mode}"
// 	publishDir "${params.results_output_dir}/vcf", mode: 'copy', overwrite: 'true', pattern: '*.vcf.gz*'
// 	memory '10 GB'
// 	time '2h'

// 	input:
// 		tuple val(group), val(type), path(sv_vcf), path(sv_tbi), path(ped), path(snv_vcf), path(snv_tbi)

// 	output:
// 		tuple val(group), path("${group_score}.snv.rescored.sorted.vcf.gz"), path("${group_score}.snv.rescored.sorted.vcf.gz.tbi"), emit: vcf_yaml
// 		tuple val(group), path("${group}_svp.INFO"), emit: svcompound_INFO
// 		path "*versions.yml", emit: versions


// 	when:
// 		params.mode == "family" && params.assay == "wgs"


// 	script:
// 		group_score = ( type == "ma" || type == "fa" ) ? "${group}_${type}" : group

// 		"""
// 		compound_finder.pl \\
// 			--sv $sv_vcf --ped $ped --snv $snv_vcf \\
// 			--osv ${group_score}.sv.rescored.sorted.vcf \\
// 			--osnv ${group_score}.snv.rescored.sorted.vcf \\
// 			--skipsv
// 		bgzip -@ ${task.cpus} ${group_score}.snv.rescored.sorted.vcf -f
// 		tabix ${group_score}.snv.rescored.sorted.vcf.gz -f
// 		echo "SVc	$type	${params.accessdir}/vcf/${group_score}.sv.scored.sorted.vcf.gz,${params.accessdir}/vcf/${group_score}.snv.rescored.sorted.vcf.gz" > ${group}_svp.INFO

// 		${compound_finder_version(task)}
// 		"""

// 	stub:
// 		group_score = group
// 		"""
// 		touch "${group_score}.snv.rescored.sorted.vcf.gz"
// 		touch "${group_score}.snv.rescored.sorted.vcf.gz.tbi"
// 		touch "${group}_svp.INFO"

// 		${compound_finder_version(task)}
// 		"""
// }
// def compound_finder_version(task) {
// 	"""
// 	cat <<-END_VERSIONS > ${task.process}_versions.yml
// 	${task.process}:
// 	    tabix: \$(echo \$(tabix --version 2>&1) | sed 's/^.*(htslib) // ; s/ Copyright.*//')
// 	END_VERSIONS
// 	"""
// }


// process output_files {
// 	cpus 2
// 	memory '1GB'
// 	time '1h'

// 	input:
// 		tuple val(group), files

// 	output:
// 		tuple val(group), path("${group}.INFO"), emit: yaml_INFO

// 	script:
// 		files = files.join( ' ' )

// 		"""
// 		cat $files > ${group}.INFO
// 		"""

// 	stub:
// 		"""
// 		touch "${group}.INFO"
// 		"""
// }


// process svvcf_to_bed {
// 	publishDir "${params.results_output_dir}/bed", mode: 'copy' , overwrite: 'true'
// 	tag "group"
// 	memory '1 GB'
// 	time '1h'
// 	cpus 2

// 	input:
// 		tuple val(group), path(vcf)
// 		tuple val(group), val(id), sex, type

// 	output:
// 		path("${group}.sv.bed")

// 	when:
// 		params.antype != "panel"


// 	script:
// 		"""
// 		cnv2bed.pl --cnv $vcf --pb $id > ${group}.sv.bed
// 		"""

// 	stub:
// 		"""
// 		touch "${group}.sv.bed"
// 		"""
// }

// process plot_pod {
// 	container  "${params.container_pod}"
// 	publishDir "${params.results_output_dir}/pod", mode: 'copy' , overwrite: 'true'
// 	tag "$group"
// 	time '1h'
// 	memory '1 GB'
// 	cpus 2

// 	input:
// 		tuple val(group), path(snv)
// 		tuple val(group), path(cnv), val(type), path(ped)
// 		tuple val(group), val(id), val(sex), val(type)

// 	output:
// 		tuple path("${id}_POD_karyotype.pdf"), path("${id}_POD_results.html")

// 	when:
// 		params.mode == "family" && params.trio

// 	script:
// 		"""
// 		parental_origin_of_duplication.pl --snv $snv --cnv $cnv --proband $id --ped $ped
// 		"""

// 	stub:
// 		"""
// 		touch "${id}_POD_karyotype.pdf"
// 		touch "${id}_POD_results.html"
// 		"""
// }

// process create_yaml {
// 	publishDir "${params.results_output_dir}/yaml", mode: 'copy' , overwrite: 'true', pattern: '*.yaml'
// 	publishDir "${params.results_output_dir}/yaml/alt_affect", mode: 'copy' , overwrite: 'true', pattern: '*.yaml.*a'
// 	publishDir "${params.crondir}/scout", mode: 'copy' , overwrite: 'true', pattern: '*.yaml'
// 	errorStrategy 'retry'
// 	maxErrors 5
// 	tag "$group"
// 	time '5m'
// 	memory '1 GB'

// 	input:
// 		tuple val(group), val(id), val(sex), val(mother), val(father), val(phenotype), val(diagnosis), val(type), val(assay), val(clarity_sample_id), val(ffpe), val(analysis), val(type), path(ped), path(INFO)

// 	output:
// 		tuple val(group), path("${group}.yaml*"), emit: scout_yaml

// 	script:
// 		"""
// 		create_yml.pl \\
// 			--g $group,$clarity_sample_id \\
// 			--d $diagnosis \\
// 			--panelsdef $params.panelsdef \\
// 			--out ${group}.yaml \\
// 			--ped $ped \\
// 			--files $INFO \\
// 			--assay $assay,$analysis \\
// 			--antype $params.antype \\
// 			--extra_panels $params.extra_panels
// 		"""

// 	stub:
// 		"""
// 		touch "${group}.yaml"
// 		"""
// }

// process combine_versions {
// 	publishDir "${params.results_output_dir}/versions", mode: 'copy', overwrite: 'true', pattern: '*.versions.yml'

// 	// The point of "first" here is that when a process is present in multiple instances
// 	// there is no need to include more than one instance of the versions
// 	input:
// 		tuple val(group), versions

// 	output:
// 		path("${group}.versions.yml")

// 	script:
// 		// versions_joined = versions.sort( my_it -> my_it.name ).join(" ")
// 		"""
// 		cat $versions_joined > ${group}.versions.yml
// 		"""

// 	stub:
// 		// versions_joined = versions.sort( my_it -> my_it.name ).join(" ")
// 		"""
// 		cat $versions_joined > ${group}.versions.yml
// 		"""
// }
