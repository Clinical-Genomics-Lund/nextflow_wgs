#!/usr/bin/env nextflow

// GENERAL PATHS //
OUTDIR = "${params.resultsdir}/${params.subdir}"
// OUTDIR = params.outdir+'/'+params.subdir
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
println("mode: "+mode)
println("trio: "+trio)
// Print commit-version of active deployment
file(params.git)
	.readLines()
	.each { println "git commit-hash: "+it }
// Print active container
container = file(params.container).toRealPath()
println("container: "+container)

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
	logFile = file("${CRONDIR}/outlog.complete")
	logFile.text = msg
	logFile.append(error)
}

// Input channels for alignment, variant calling and annotation //
Channel
	.fromPath(params.csv)
	.splitCsv(header:true)
	.map{ row-> tuple(row.group, 
				row.id, 
				(row.containsKey("bam") ? file(row.bam) : (row.containsKey("vcf") ? file(row.vcf) : file(row.read1) ) ), 
				(row.containsKey("bai") ? file(row.bai) : (row.containsKey("idx") ? file(row.idx) : file(row.read2) ) ) ) }
	.set { input_files }

fastq = Channel.create()
bam_choice = Channel.create()
//vcf_choice = Channel.create()
fastq_sharded = Channel.create()
fastq_umi = Channel.create()
annotate_only = Channel.create()


// If input-files has bam files bypass alignment, otherwise go for fastq-channels => three options for fastq, sharded bwa, normal bwa or umi trimming
input_files.view().choice(bam_choice, fastq, fastq_sharded, fastq_umi, annotate_only ) { it[2] =~ /\.bam/ ? 0 : ( it[2] =~ /\.vcf.gz/ ? 4 : (params.shardbwa ? 2 : (params.umi ? 3 : 1) )) }

annotate_only.into{
	annotate_only_vep;
	annotate_only_cadd
}

bam_choice.into{ 
	expansionhunter_bam_choice; 
	dnascope_bam_choice;
	bampath_start;
	chanjo_bam_choice; 
	yaml_bam_choice; 
	cov_bam_choice; 
	bam_manta_choice; 
	bam_nator_choice; 
	bam_tiddit_choice; 
	bam_mito_choice; 
	bam_SMN_choice; 
	bam_freebayes_choice;
	bam_mantapanel_choice;
	bam_cnvkitpanel_choice;
	bam_dellypanel_choice;
	bam_melt_choice;
	bam_qc_choice;
	dedup_dummy_choice;
	bam_bqsr_choice;
	bam_gatk_choice }

// vcf_choice.into{
// 	split_cadd_choice;
// 	split_vep_choice;
// }

// For melt to work if started from bam-file.
process dedupdummy {
	when:
		params.onco

	input:
		set id, group, file(bam), file(bai) from dedup_dummy_choice
	output:
		set id, file("dummy") into dedup_dummy
	"""
	echo test > dummy
	"""
}
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
	.into { meta_gatkcov; meta_exp; meta_svbed; meta_pod; meta_mutect2}


Channel
	.fromPath(params.gatkreffolders)
	.splitCsv(header:true)
	.map{ row-> tuple(row.i, row.refpart) }
	.into{ gatk_ref; gatk_postprocess }


// Check whether genome assembly is indexed //
if(genome_file ){
	bwaId = Channel
			.fromPath("${genome_file}.bwt")
			.ifEmpty { exit 1, "BWA index not found: ${genome_file}.bwt" }
}

process fastp {
	publishDir "/${OUTDIR}/versions", mode: 'copy' , overwrite: 'true', pattern: '*versions.yml'
	cpus 10
	tag "$id"
	time '1h'
	memory '20 GB'
	scratch true
	stageInMode 'copy'
	stageOutMode 'copy'

	when:
		params.umi

	input:
		set group, val(id), r1, r2 from fastq_umi

	output:
		set group, val(id), file("${id}_R1_a_q_u_trimmed.fq.gz"),file("${id}_R2_a_q_u_trimmed.fq.gz") into fastq_trimmed
		path "*versions.yml"

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
		
		${fastp_version(task)}
		"""

	stub:
		"""
		touch ${id}_R1_a_q_u_trimmed.fq.gz
		touch ${id}_R2_a_q_u_trimmed.fq.gz

		${fastp_version(task)}
		"""
}
def fastp_version(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	fastp: 
	 version: \$(echo \$(fastp -v 2>&1) | cut -f 2 -d " ")
	 container: ${task.container}
	END_VERSIONS
	"""
}


// Align fractions of fastq files with BWA
process bwa_align_sharded {
	publishDir "/${OUTDIR}/versions", mode: 'copy' , overwrite: 'true', pattern: '*versions.yml'
	cpus 50
	memory '120 GB'
	tag "$id $shard"
	time '5h'

	input:
		set val(shard), val(group), val(id), r1, r2 from bwa_shards.combine(fastq_sharded)

	output:
		set val(id), group, file("${id}_${shard}.bwa.sort.bam"), file("${id}_${shard}.bwa.sort.bam.bai") into bwa_shards_ch
		path "*versions.yml"

	when:
		params.align && params.shardbwa

	script:
		"""
		sentieon bwa mem -M \\
			-R '@RG\\tID:${id}\\tSM:${id}\\tPL:illumina' \\
			-K $K_size \\
			-t ${task.cpus} \\
			-p $genome_file '<sentieon fqidx extract -F $shard/$bwa_num_shards -K $K_size $r1 $r2' | sentieon util sort \\
			-r $genome_file \\
			-o ${id}_${shard}.bwa.sort.bam \\
			-t ${task.cpus} --sam2bam -i -
		
		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 Sentieon UTIL: 
		  version: \$(echo \$(sentieon util --version 2>&1) | sed -e "s/sentieon-genomics-//g")
		  container: ${task.container}
		 Sentieon BWA-MEM: 
		  version: \$(echo \$(sentieon bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
		  container: ${task.container}
		END_VERSIONS
		"""

	stub:
		"""
		touch ${id}_${shard}.bwa.sort.bam
		touch ${id}_${shard}.bwa.sort.bam.bai

		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 Sentieon UTIL: 
		  version: \$(echo \$(sentieon util --version 2>&1) | sed -e "s/sentieon-genomics-//g")
		  container: ${task.container}
		 Sentieon BWA-MEM: 
		  version: \$(echo \$(sentieon bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
		  container: ${task.container}
		END_VERSIONS
		"""
}

// Merge the fractioned bam files
process bwa_merge_shards {
	publishDir "/${OUTDIR}/versions", mode: 'copy' , overwrite: 'true', pattern: '*versions.yml'
	cpus 50
	tag "$id"
	time '1h'
	memory '120 GB'

	input:
		set val(id), group, file(shard), file(shard_bai) from bwa_shards_ch.groupTuple(by: [0,1])

	output:
		set id, group, file("${id}_merged.bam"), file("${id}_merged.bam.bai") into merged_bam_locusc
		set id, file("${id}_merged.bam"), file("${id}_merged.bam.bai") into merged_bam_dedup
		path "*versions.yml"
	when:
		params.shardbwa
	
	script:
		bams = shard.sort(false) { a, b -> a.getBaseName() <=> b.getBaseName() } .join(' ')

		"""
		sentieon util merge -o ${id}_merged.bam ${bams}

		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 Sentieon UTIL: 
		  version: \$(echo \$(sentieon util --version 2>&1) | sed -e "s/sentieon-genomics-//g")
		  container: ${task.container}
		END_VERSIONS
		"""

	stub:
		"""
		touch ${id}_merged.bam
		touch ${id}_merged.bam.bai

		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 Sentieon UTIL: 
		  version: \$(echo \$(sentieon util --version 2>&1) | sed -e "s/sentieon-genomics-//g")
		  container: ${task.container}
		END_VERSIONS
		"""
}

// ALTERNATIVE PATH: Unsharded BWA, utilize local scratch space.
process bwa_align {
	publishDir "/${OUTDIR}/versions", mode: 'copy' , overwrite: 'true', pattern: '*versions.yml'
	cpus 50
	memory '80 GB'
	// 64 GB peak giab //
	scratch true
	stageInMode 'copy'
	stageOutMode 'copy'
	tag "$id"
	container = "/fs1/resources/containers/sentieon_202112.sif"

	input:
		set val(group), val(id), file(r1), file(r2) from fastq.mix(fastq_trimmed)

	output:
		set id, group, file("${id}_merged.bam"), file("${id}_merged.bam.bai") into bam_locusc, bam_markdup
		// set id, file("${id}_merged.bam"), file("${id}_merged.bam.bai") into bam_dedup remnant of distri
		path "*versions.yml"

	when:
		params.align && !params.shardbwa

	script:
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
		
		${bwa_align_versions(task)}
		"""

	stub:
		"""
		touch ${id}_merged.bam
		touch ${id}_merged.bam.bai

		${bwa_align_versions(task)}
		"""

}
def bwa_align_versions(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	 Sentieon UTIL: 
	  version: \$(echo \$(sentieon util --version 2>&1) | sed -e "s/sentieon-genomics-//g")
	  container: ${task.container}
	 Sentieon BWA-MEM: 
	  version: \$(echo \$(sentieon bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
	  container: ${task.container}
	END_VERSIONS
	"""
}


process markdup {
	cpus 40
	errorStrategy 'retry'
	maxErrors 5
	tag "$id"
	memory '30 GB'
	// 12gb peak giab //
	time '3h'
	scratch true
	stageInMode 'copy'
	stageOutMode 'copy'
	container = "/fs1/resources/containers/sentieon_202112.sif"
	publishDir "/${OUTDIR}/bam", mode: 'copy' , overwrite: 'true', pattern: '*_dedup.bam*'
	publishDir "/${OUTDIR}/versions", mode: 'copy' , overwrite: 'true', pattern: '*versions.yml'

	input:
		set id, group, file(bam), file(bai) from bam_markdup.mix(merged_bam_dedup)

	output:
		set group, id, file("${id}_dedup.bam"), file("${id}_dedup.bam.bai") into complete_bam, chanjo_bam, expansionhunter_bam, yaml_bam, cov_bam, bam_manta, bam_nator, bam_tiddit, bam_manta_panel, bam_delly_panel, bam_cnvkit_panel, bam_freebayes, bam_mito, smncnc_bam, bam_gatk, depth_onco
		set id, group, file("${id}_dedup.bam"), file("${id}_dedup.bam.bai") into qc_bam, bam_melt, bam_bqsr
		set val(id), file("dedup_metrics.txt") into dedupmet_sentieonqc
		set group, file("${group}_bam.INFO") into bam_INFO
		path "*versions.yml"

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
		echo "BAM	$id	/access/${params.subdir}/bam/${id}_dedup.bam" > ${group}_bam.INFO

		${markdup_versions(task)}
		"""

	stub:
		"""
		touch ${id}_dedup.bam
		touch ${id}_dedup.bam.bai
		touch "dedup_metrics.txt"
		touch "${group}_bam.INFO"

		${markdup_versions(task)}
		"""
}
def markdup_versions(task) {
	"""
	cat <<-END_VERSIONS > ${task.process}_versions.yml
	${task.process}:
	 Sentieon DRIVER: 
	  version: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
	  container: ${task.container}
	END_VERSIONS
	"""
}


process bqsr {
	cpus 40
	errorStrategy 'retry'
	maxErrors 5
	tag "$id"
	memory '30 GB'
	// 12gb peak giab //
	time '5h'
	scratch true
	stageInMode 'copy'
	stageOutMode 'copy'
	container = "/fs1/resources/containers/sentieon_202112.sif"
	publishDir "/${OUTDIR}/bqsr", mode: 'copy' , overwrite: 'true', pattern: '*table'
	publishDir "/${OUTDIR}/versions", mode: 'copy' , overwrite: 'true', pattern: '*versions.yml'

	input:
		set id, group, file(bam), file(bai) from bam_bqsr.mix(bam_bqsr_choice)

	output:
		set group, id, file("${id}.bqsr.table") into dnascope_bqsr
		path "*versions.yml"

	script:
		"""
		sentieon driver -t ${task.cpus} \\
			-r $genome_file -i $bam \\
			--algo QualCal ${id}.bqsr.table \\
			-k $params.KNOWN

		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 Sentieon DRIVER: 
		  version: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
		  container: ${task.container}
		END_VERSIONS
		"""

	stub:
		"""
		touch ${id}.bqsr.table

		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 Sentieon DRIVER: 
		  version: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
		  container: ${task.container}
		END_VERSIONS
		"""
}

//Collect various QC data: 
process sentieon_qc {
	publishDir "/${OUTDIR}/versions", mode: 'copy' , overwrite: 'true', pattern: '*versions.yml'
	cpus 52
	memory '30 GB'
	tag "$id"
	time '2h'
	scratch true
	stageInMode 'copy'
	stageOutMode 'copy'
	container = "/fs1/resources/containers/sentieon_202112.sif"

	input:
		set id, group, file(bam), file(bai), file(dedup) from qc_bam.mix(bam_qc_choice).join(dedupmet_sentieonqc.mix(dedup_dummy))

	output:
		set group, id, file("${id}_qc.json") into qc_cdm
		set group, id, file("${id}_qc.json") into qc_melt
		file("*.txt")
		path "*versions.yml"

	script:
		target = ""
		panel = ""
		cov = "WgsMetricsAlgo wgs_metrics.txt"
		assay = "wgs"
		if( params.onco || params.exome) {
			target = "--interval $params.intervals"
			panel = params.panelhs + "$bam" + params.panelhs2 
			cov = "CoverageMetrics --cov_thresh 1 --cov_thresh 10 --cov_thresh 30 --cov_thresh 100 --cov_thresh 250 --cov_thresh 500 cov_metrics.txt"
			assay = "panel"
		}
		
		"""
		sentieon driver \\
			-r $genome_file $target \\
			-t ${task.cpus} \\
			-i $bam \\
			--algo MeanQualityByCycle mq_metrics.txt \\
			--algo QualDistribution qd_metrics.txt \\
			--algo GCBias --summary gc_summary.txt gc_metrics.txt \\
			--algo AlignmentStat aln_metrics.txt \\
			--algo InsertSizeMetricAlgo is_metrics.txt \\
			--algo $cov
		$panel
		qc_sentieon.pl $id $assay > ${id}_qc.json

		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 Sentieon DRIVER: 
		  version: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
		  container: ${task.container}
		END_VERSIONS
		"""

	stub:
		"""
		touch ${id}_qc.json
		touch ${id}.txt
		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 Sentieon DRIVER: 
		  version: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
		  container: ${task.container}
		END_VERSIONS
		"""

}

   

// Calculate coverage for chanjo
process chanjo_sambamba {
	cpus 16
	memory '10 GB'
	publishDir "/${OUTDIR}/cov", mode: 'copy', overwrite: 'true', pattern: '*.cov'
	publishDir "/${OUTDIR}/versions", mode: 'copy' , overwrite: 'true', pattern: '*versions.yml'
	tag "$id"
	scratch true
	stageInMode 'copy'
	stageOutMode 'copy'

	when:
		params.varcall

	input:	
		set group, id, file(bam), file(bai) from chanjo_bam.mix(chanjo_bam_choice)

	output:
		file("${id}.bwa.chanjo.cov") into chanjocov

	script:
		"""
		sambamba depth region -t ${task.cpus} -L $params.scoutbed -T 10 -T 15 -T 20 -T 50 -T 100 $bam > ${id}.bwa.chanjo.cov

		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 Sambamba: 
		  version: \$(echo \$(sambamba --version 2>&1) | awk '{print \$2}' )
		  container: ${task.container}
		END_VERSIONS
		"""

	stub:
		"""
		touch ${id}.bwa.chanjo.cov

		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 Sambamba: 
		  version: \$(echo \$(sambamba --version 2>&1) | awk '{print \$2}' )
		  container: ${task.container}
		END_VERSIONS
		"""
}																																																								

// Calculate coverage for paneldepth
process depth_onco {
	cpus 2
	memory '10 GB'
	publishDir "/${OUTDIR}/cov", mode: 'copy', overwrite: 'true'
	tag "$id"
	scratch true
	stageInMode 'copy'
	stageOutMode 'copy'

	when:
		params.assay == "swea"

	input:	
		set group, id, file(bam), file(bai) from depth_onco

	output:
		file("${id}.lowcov.overlapping.bed") into cov_onco

	script:
		"""
		panel_depth.pl $bam $params.scoutbed > ${id}.lowcov.bed
		overlapping_genes.pl ${id}.lowcov.bed $params.gene_regions > ${id}.lowcov.overlapping.bed
		"""

	stub:
		"""
		touch ${id}.lowcov.overlapping.bed
		"""
}

process SMNCopyNumberCaller {
	cpus 10
	memory '25GB'
	time '2h'
	publishDir "/${OUTDIR}/plots/SMNcnc", mode: 'copy' , overwrite: 'true', pattern: '*.pdf*'
	publishDir "/${OUTDIR}/versions", mode: 'copy' , overwrite: 'true', pattern: '*versions.yml'
	tag "$id"

	when:
		params.antype == "wgs"

	input:
		set group, id, file(bam), file(bai) from smncnc_bam.mix(bam_SMN_choice)

	output:
		file("*.tsv") into smn_tsv
		set file("*.pdf"), file("*.json")
		set group, file("${group}_smn.INFO") into smn_INFO
		path "*versions.yml"

	script:
		"""
		samtools view -H $bam | \\
			sed -e 's/SN:1/SN:chr1/' | sed -e 's/SN:2/SN:chr2/' |  \\
			sed -e 's/SN:3/SN:chr3/' | sed -e 's/SN:4/SN:chr4/' |  \\
			sed -e 's/SN:5/SN:chr5/' | sed -e 's/SN:6/SN:chr6/' |  \\
			sed -e 's/SN:7/SN:chr7/' | sed -e 's/SN:8/SN:chr8/' |  \\
			sed -e 's/SN:9/SN:chr9/' | sed -e 's/SN:10/SN:chr10/' | \\
			sed -e 's/SN:11/SN:chr11/' | sed -e 's/SN:12/SN:chr12/' |  \\
			sed -e 's/SN:13/SN:chr13/' | sed -e 's/SN:14/SN:chr14/' |  \\
			sed -e 's/SN:15/SN:chr15/' | sed -e 's/SN:16/SN:chr16/' |  \\
			sed -e 's/SN:17/SN:chr17/' | sed -e 's/SN:18/SN:chr18/' |  \\
			sed -e 's/SN:19/SN:chr19/' | sed -e 's/SN:20/SN:chr20/' |  \\
			sed -e 's/SN:21/SN:chr21/' | sed -e 's/SN:22/SN:chr22/' |  \\
			sed -e 's/SN:X/SN:chrX/' | sed -e 's/SN:Y/SN:chrY/' |   \\
			sed -e 's/SN:MT/SN:chrM/' | \\
			samtools reheader - $bam > ${id}.bam
		samtools index -b ${id}.bam -@ ${task.cpus}
		echo ${id}.bam > manifest.txt
		smn_caller.py --manifest manifest.txt --genome 38 --prefix ${id} --outDir . --threads ${task.cpus}
		rm ${id}.bam
		source activate py3-env
		python /SMNCopyNumberCaller/smn_charts.py -s ${id}.json -o .
		mv ${id}.tsv ${group}_SMN.tsv
		echo "SMN ${params.accessdir}/smn/${group}_SMN.tsv" > ${group}_smn.INFO

		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 SAMtools: 
		  version: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
		  container: ${task.container}
		 SMNCopyNumberCaller: 
		  version: 1.1.2
		  container: ${task.container}
		END_VERSIONS
		"""

	stub:
		"""
		touch ${id}.bam
		touch ${id}.tsv
		touch ${id}.pdf
		touch ${id}.json
		touch ${group}_SMN.tsv
		touch ${group}_smn.INFO
	
		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 SAMtools: 
		  version: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
		  container: ${task.container}
		 SMNCopyNumberCaller: 
		  version: 1.1.2
		  container: ${task.container}
		END_VERSIONS
		"""
}
// collects each individual's SMNCNC-tsv and creates one tsv-file
smn_tsv
	.collectFile(keepHeader: true, storeDir: "${OUTDIR}/smn/")



////////////////////////////////////////////////////////////////////////
////////////////////////// EXPANSION HUNTER ////////////////////////////
////////////////////////////////////////////////////////////////////////

// call STRs using ExpansionHunter, and plot alignments with GraphAlignmentViewer
process expansionhunter {
	publishDir "/${OUTDIR}/versions", mode: 'copy' , overwrite: 'true', pattern: '*versions.yml'
	tag "$group"
	cpus 2
	time '10h'
	memory '40 GB'
	scratch true
	stageInMode 'copy'
	stageOutMode 'copy'
	//publishDir "/${OUTDIR}/plots/GAV/${group}", mode: 'copy' , overwrite: 'true', pattern: '*.png'

	when:
		params.str
		
	input:
		set group, id, file(bam), file(bai), sex, type \
			from expansionhunter_bam.mix(expansionhunter_bam_choice).join(meta_exp, by: [0,1]).filter { item -> item[5] == 'proband' }

	output:
		set group, id, file("${group}.eh.vcf") into expansionhunter_vcf
		set group, id, file("${group}.eh_realigned.sort.bam"), file("${group}.eh_realigned.sort.bam.bai"), file("${group}.eh.vcf") into reviewer
		path "*versions.yml"

	script:
		"""
		source activate htslib10
		ExpansionHunter \
			--reads $bam \
			--reference $genome_file \
			--variant-catalog $params.expansionhunter_catalog \
			--output-prefix ${group}.eh
		samtools sort ${group}.eh_realigned.bam -o ${group}.eh_realigned.sort.bam
		samtools index ${group}.eh_realigned.sort.bam

		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 expansionhunter: 
		  version: \$(echo \$(ExpansionHunter --version 2>&1) | sed 's/.*ExpansionHunter v// ; s/]//')
		  container: ${task.container}
		 SAMtools: 
		  version: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
		  container: ${task.container}
		END_VERSIONS
		"""

	stub:
		"""
		source activate htslib10
		touch ${group}.eh.vcf
		touch ${group}.eh_realigned.sort.bam
		touch ${group}.eh_realigned.sort.bam.bai

		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 expansionhunter: 
		  version: \$(echo \$(ExpansionHunter --version 2>&1) | sed 's/.*ExpansionHunter v// ; s/]//')
		  container: ${task.container}
		 SAMtools: 
		  version: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
		  container: ${task.container}
		END_VERSIONS
		"""
}

// annotate expansionhunter vcf
process stranger {
	publishDir "/${OUTDIR}/versions", mode: 'copy' , overwrite: 'true', pattern: '*versions.yml'
	tag "$group"
	memory '1 GB'
	time '10m'
	scratch true
	stageInMode 'copy'
	stageOutMode 'copy'
	container = "/fs1/resources/containers/stranger_0.8.sif"

	input:
		set group, id, file(eh_vcf) from expansionhunter_vcf

	output:
		set group, id, file("${group}.fixinfo.eh.stranger.vcf") into expansionhunter_vcf_anno
		path "*versions.yml"

	script:
		"""
		stranger ${eh_vcf} -f $params.expansionhunter_catalog > ${group}.eh.stranger.vcf
		grep ^# ${group}.eh.stranger.vcf > ${group}.fixinfo.eh.stranger.vcf
		grep -v ^# ${group}.eh.stranger.vcf | sed 's/ /_/g' >> ${group}.fixinfo.eh.stranger.vcf
		
		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 stranger: 
		  version: \$( stranger --version )
		  container: ${task.container}
		END_VERSIONS
		"""

	stub:
		"""
		touch ${group}.fixinfo.eh.stranger.vcf

		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 stranger: 
		  version: \$( stranger --version )
		  container: ${task.container}
		END_VERSIONS
		"""
}
//for i in $( ls *.svg | cut -f 2 -d "." ); do echo "STR_IMG $i /access/!{params.subdir}/plots/reviewer/!{group}/!{group}.${i}.svg" >> !{group}_rev.INFO; done
process reviewer {
	tag "$group"
	cpus 1
	time '10m'
	memory '1 GB'
	scratch true
	stageInMode 'copy'
	stageOutMode 'copy'
	errorStrategy 'ignore'
	container = "/fs1/resources/containers/REViewer_2021-06-07.sif"
	publishDir "/${OUTDIR}/plots/reviewer/${group}", mode: 'copy' , overwrite: 'true', pattern: '*.svg'
	publishDir "/${OUTDIR}/versions", mode: 'copy' , overwrite: 'true', pattern: '*versions.yml'
	
	input:
		set group, id, file(bam), file(bai), file(vcf) from reviewer

	output:
		file("*svg")
		//set group, file("${group}_rev.INFO") into reviewer_INFO
		path "*versions.yml"

	shell:
		'''
		grep LocusId !{params.expansionhunter_catalog} | sed 's/[",^ ]//g' | cut -d':' -f2 | perl -na -e 'chomp; \
		system("REViewer --reads !{bam} \
		--vcf !{vcf} \
		--reference !{genome_file} \
		--catalog !{params.expansionhunter_catalog} \
		--locus $_ \
		--output-prefix !{id}");'

		cat <<-END_VERSIONS > !{task.process}_versions.yml
		!{task.process}:
		 REViewer: 
		  version: $(echo $(REViewer --version 2>&1) | sed 's/^.*REViewer v//')
		  container: !{task.container}
		END_VERSIONS
		'''

	stub:
		"""
		touch ${id}.svg

		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 REViewer: 
		  version: \$(echo \$(REViewer --version 2>&1) | sed 's/^.*REViewer v//')
		  container: ${task.container}
		END_VERSIONS
		"""
}

// split multiallelic sites in expansionhunter vcf
// FIXME: Use env variable for picard path...
process vcfbreakmulti_expansionhunter {
	publishDir "/${OUTDIR}/vcf", mode: 'copy' , overwrite: 'true', pattern: '*.expansionhunter.vcf.gz'
	publishDir "/${OUTDIR}/versions", mode: 'copy' , overwrite: 'true', pattern: '*versions.yml'
	tag "$group"
	time '10m'
	memory '40 GB'
	scratch true
	stageInMode 'copy'
	stageOutMode 'copy'

	input:
		set group, id, file(eh_vcf_anno) from expansionhunter_vcf_anno
		set group, id, sex, mother, father, phenotype, diagnosis, type, assay, clarity_sample_id, ffpe, analysis from meta_str.filter{ item -> item[7] == 'proband' }

	output:
		file("${group}.expansionhunter.vcf.gz") into expansionhunter_scout
		set group, file("${group}_str.INFO") into str_INFO
		path "*versions.yml"

	script:
		if (father == "") { father = "null" }
		if (mother == "") { mother = "null" }
		if (mode == "family") {
			"""
			java -jar /opt/conda/envs/CMD-WGS/share/picard-2.21.2-1/picard.jar RenameSampleInVcf INPUT=${eh_vcf_anno} OUTPUT=${eh_vcf_anno}.rename.vcf NEW_SAMPLE_NAME=${id}
			vcfbreakmulti ${eh_vcf_anno}.rename.vcf > ${group}.expansionhunter.vcf.tmp
			familyfy_str.pl --vcf ${group}.expansionhunter.vcf.tmp --mother $mother --father $father --out ${group}.expansionhunter.vcf
			bgzip ${group}.expansionhunter.vcf
			tabix ${group}.expansionhunter.vcf.gz
			echo "STR	${params.accessdir}/vcf/${group}.expansionhunter.vcf.gz" > ${group}_str.INFO

			cat <<-END_VERSIONS > ${task.process}_versions.yml
			${task.process}:
			 vcflib: 
			  version: 1.0.9
			  container: ${task.container}
			 RenameSampleInVcf: 
			  version: \$(echo \$(java -jar /opt/conda/envs/CMD-WGS/share/picard-2.21.2-1/picard.jar RenameSampleInVcf --version 2>&1) | sed 's/-SNAPSHOT//')
			  container: ${task.container}
			 tabix: 
			  version: \$(echo \$(tabix --version 2>&1) | sed 's/^.*(htslib) // ; s/ Copyright.*//')
			  container: ${task.container}
			 bgzip: 
			  version: \$(echo \$(bgzip --version 2>&1) | sed 's/^.*(htslib) // ; s/ Copyright.*//')
			  container: ${task.container}
			END_VERSIONS
			"""
		}
		else {
			"""
			java -jar /opt/conda/envs/CMD-WGS/share/picard-2.21.2-1/picard.jar RenameSampleInVcf INPUT=${eh_vcf_anno} OUTPUT=${eh_vcf_anno}.rename.vcf NEW_SAMPLE_NAME=${id}
			vcfbreakmulti ${eh_vcf_anno}.rename.vcf > ${group}.expansionhunter.vcf
			bgzip ${group}.expansionhunter.vcf
			tabix ${group}.expansionhunter.vcf.gz
			echo "STR	${params.accessdir}/vcf/${group}.expansionhunter.vcf.gz" > ${group}_str.INFO

			cat <<-END_VERSIONS > ${task.process}_versions.yml
			${task.process}:
			 vcflib: 
			  version: 1.0.9
			  container: ${task.container}
			 RenameSampleInVcf: 
			  version: \$(echo \$(java -jar /opt/conda/envs/CMD-WGS/share/picard-2.21.2-1/picard.jar RenameSampleInVcf --version 2>&1) | sed 's/-SNAPSHOT//')
			  container: ${task.container}
			 tabix: 
			  version: \$(echo \$(tabix --version 2>&1) | sed 's/^.*(htslib) // ; s/ Copyright.*//')
			  container: ${task.container}
			 bgzip: 
			  version: \$(echo \$(bgzip --version 2>&1) | sed 's/^.*(htslib) // ; s/ Copyright.*//')
			  container: ${task.container}
			END_VERSIONS
			"""
		}

	stub:
		"""
		touch ${group}.expansionhunter.vcf.gz
		touch "${group}_str.INFO"

		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 vcflib: 
		  version: 1.0.9
		  container: ${task.container}
		 RenameSampleInVcf: 
		  version: \$(echo \$(java -jar /opt/conda/envs/CMD-WGS/share/picard-2.21.2-1/picard.jar RenameSampleInVcf --version 2>&1) | sed 's/-SNAPSHOT//')
		  container: ${task.container}
		 tabix: 
		  version: \$(echo \$(tabix --version 2>&1) | sed 's/^.*(htslib) // ; s/ Copyright.*//')
		  container: ${task.container}
		 bgzip: 
		  version: \$(echo \$(bgzip --version 2>&1) | sed 's/^.*(htslib) // ; s/ Copyright.*//')
		  container: ${task.container}
		END_VERSIONS
		"""
}
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
process melt_qc_val {
	tag "$id"
	time '2m'
	memory '50 MB'

	when:
		params.onco

	input:
		set group, id, qc from qc_melt

	output:
		set id, val(INS_SIZE), val(MEAN_DEPTH), val(COV_DEV) into qc_melt_val
		set group, id, val(INS_SIZE), val(MEAN_DEPTH), val(COV_DEV) into qc_cnvkit_val
	
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

	stub:
		INS_SIZE = 0
		MEAN_DEPTH = 0
		COV_DEV = 0
		"""
		echo $INS_SIZE $MEAN_DEPTH $COV_DEV > qc.val
		"""
}

// MELT always give VCFs for each type of element defined in mei_list
// If none found -> 0 byte vcf. merge_melt.pl merges the three, if all empty
// it creates a vcf with only header from params.meltheader
// merge_melt.pl gives output ${id}.melt.merged.vcf
process melt {
	publishDir "/${OUTDIR}/versions", mode: 'copy' , overwrite: 'true', pattern: '*versions.yml'
	cpus 3
	errorStrategy 'retry'
	container = '/fs1/resources/containers/melt_2.2.2.sif'
	tag "$id"
	memory '40 GB'
	time '3h'
	publishDir "${OUTDIR}/vcf", mode: 'copy' , overwrite: 'true'

	input:
		set id, group, file(bam), file(bai), val(INS_SIZE), val(MEAN_DEPTH), val(COV_DEV) from bam_melt.mix(bam_melt_choice).join(qc_melt_val)

	when:
		params.onco

	output:
		set group, id, file("${id}.melt.merged.vcf") into melt_vcf_nonfiltered
		path "*versions.yml"

	script:
		"""
		set group, id, file("${id}.melt.merged.vcf") into melt_vcf_nonfiltered

		java -jar /opt/MELTv2.2.2/MELT.jar Single \\
			-bamfile $bam \\
			-r 150 \\
			-h $genome_file \\
			-n /opt/MELTv2.2.2/add_bed_files/Hg38/Hg38.genes.bed \\
			-z 500000 \\
			-d 50 -t /opt/mei_list \\
			-w . \\
			-c $MEAN_DEPTH \\
			-cov $COV_DEV \\
			-e $INS_SIZE
		merge_melt.pl $params.meltheader $id

		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 MELT: 
		  version: \$(echo \$(java -jar /opt/MELTv2.2.2/MELT.jar Single -h | sed "s/.*MELTv// ; s/ -.*//")
		  container: ${task.container}
		END_VERSIONS
		"""

	stub:
		"""
		touch ${group}.expansionhunter.vcf.gz
		touch "${group}_str.INFO"

		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 MELT: 
		  version: \$(echo \$(java -jar /opt/MELTv2.2.2/MELT.jar Single -h | sed "s/.*MELTv// ; s/ -.*//")
		  container: ${task.container}
		END_VERSIONS
		"""
}

process intersect_melt {
	cpus 2
	tag "$id"
	memory '2 GB'
	time '1h'
	publishDir "${OUTDIR}/vcf", mode: 'copy' , overwrite: 'true'

	input:
		set group, id, file(vcf) from melt_vcf_nonfiltered

	when:
		params.onco

	output:
		set group, id, file("${id}.melt.merged.intersected.vcf") into melt_vcf

	"""
	bedtools intersect -a $vcf -b $params.intersect_bed -header > ${id}.melt.merged.intersected.vcf
	"""

}

// When rerunning sample from bam, dnascope has to be run unsharded. this is mixed together with all other vcfs in a trio //
process dnascope {
	publishDir "/${OUTDIR}/versions", mode: 'copy' , overwrite: 'true', pattern: '*versions.yml'
	cpus 54
	memory '40 GB'
	// 12 GB peak giab //
	time '4h'
	tag "$id"
	container = "/fs1/resources/containers/sentieon_202112.sif"

	when:
		params.varcall

	input:
		set group, id, bam, bai, bqsr from complete_bam.mix(dnascope_bam_choice).join(dnascope_bqsr, by: [0,1] )

	output:
		set group, id, file("${id}.dnascope.gvcf.gz"), file("${id}.dnascope.gvcf.gz.tbi") into complete_vcf_choice
		set group, id, file("${id}.dnascope.gvcf.gz") into gvcf_gens_choice
		path "*versions.yml"

	script:
		"""
		sentieon driver \\
			-t ${task.cpus} \\
			-r $genome_file \\
			-q $bqsr \\
			-i $bam --shard 1:1-248956422 --shard 2:1-242193529 --shard 3:1-198295559 --shard 4:1-190214555 --shard 5:1-120339935 --shard 5:120339936-181538259 --shard 6:1-170805979 --shard 7:1-159345973 --shard 8:1-145138636 --shard 9:1-138394717 --shard 10:1-133797422 --shard 11:1-135086622 --shard 12:1-56232327 --shard 12:56232328-133275309 --shard 13:1-114364328 --shard 14:1-107043718 --shard 15:1-101991189 --shard 16:1-90338345 --shard 17:1-83257441 --shard 18:1-80373285 --shard 19:1-58617616 --shard 20:1-64444167 --shard 21:1-46709983 --shard 22:1-50818468 --shard X:1-124998478 --shard X:124998479-156040895 --shard Y:1-57227415 --shard M:1-16569 \\
			--algo DNAscope --emit_mode GVCF ${id}.dnascope.gvcf.gz

		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 Sentieon DRIVER: 
		  version: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
		  container: ${task.container}
		END_VERSIONS
		"""

	stub:
		"""
		touch ${id}.dnascope.gvcf.gz
		touch ${id}.dnascope.gvcf.gz.tbi

		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 Sentieon DRIVER: 
		  version: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
		  container: ${task.container}
		END_VERSIONS
		"""
}

process bamtoyaml {
	cpus 1
	time "5m"
	memory "2MB"

	input:
		set group, id, bam, bai from bampath_start
	
	output:
		set group, file("${group}_bamstart.INFO") into bamchoice_INFO

	script:
		"""
		echo "BAM	$id	/access/${params.subdir}/bam/${bam.getName()}" > ${group}_bamstart.INFO
		"""

	stub:
		"""
		touch ${group}_bamstart.INFO
		"""
}


process gvcf_combine {
	publishDir "/${OUTDIR}/versions", mode: 'copy' , overwrite: 'true', pattern: '*versions.yml'
	cpus 16
	tag "$group"
	memory '5 GB'
	time '5h'
	container = "/fs1/resources/containers/sentieon_202112.sif"

	input:
		set group, id, file(vcf), file(idx) from complete_vcf_choice.groupTuple()

	output:
		set group, id, file("${group}.combined.vcf"), file("${group}.combined.vcf.idx") into combined_vcf
		path "*versions.yml"

	script:
		all_gvcfs = vcf.join(' -v ')

		"""
		sentieon driver \\
			-t ${task.cpus} \\
			-r $genome_file \\
			--algo GVCFtyper \\
			-v $all_gvcfs ${group}.combined.vcf

		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 Sentieon DRIVER: 
		  version: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
		  container: ${task.container}
		END_VERSIONS
		"""

	stub:
		"""
		touch ${group}.combined.vcf
		touch ${group}.combined.vcf.idx

		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 Sentieon DRIVER: 
		  version: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
		  container: ${task.container}
		END_VERSIONS
		"""
}

// Create ped from input variables //
process create_ped {
	tag "$group"
	time '20m'
	publishDir "/${OUTDIR}/ped", mode: 'copy' , overwrite: 'true'	

	input:
		set group, id, sex, mother, father, phenotype, diagnosis, type, assay, clarity_sample_id, ffpe, analysis from ped.filter { item -> item[7] == 'proband' }
		
	output:
		set group, type, file("${group}_base.ped") into ped_mad, ped_peddy, ped_inher, ped_scout, ped_loqus, ped_prescore, ped_compound, ped_pod
		set group, type_ma, file("${group}_ma.ped") optional true into ped_inher_ma, ped_prescore_ma, ped_compound_ma, ped_mad_ma
		set group, type_fa, file("${group}_fa.ped") optional true into ped_inher_fa, ped_prescore_fa, ped_compound_fa, ped_mad_fa


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
		touch ${group}_base.ped
		touch ${group}_ma.ped
		touch ${group}_fa.ped

        echo $type_fa $type_ma > type.val
		"""
}

//madeline ped, run if family mode
process madeline {
	publishDir "/${OUTDIR}/ped", mode: 'copy' , overwrite: 'true', pattern: '*.madeline.xml'
	publishDir "/${OUTDIR}/versions", mode: 'copy' , overwrite: 'true', pattern: '*versions.yml'
	memory '1 GB'
	time '30m'
	container '/fs1/resources/containers/madeline.sif'

	input:
		set group, type, file(ped) from ped_mad.mix(ped_mad_ma,ped_mad_fa)

	output:
		file("${ped}.madeline.xml") into madeline_ped
		set group, file("${group}_madde.INFO") into madde_INFO
		path "*versions.yml"

	when:
		mode == "family" && params.assay == "wgs"

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

		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 ped_parser: 
		  version: \$(echo \$(ped_parser --version 2>&1) | sed -e "s/^.*ped_parser version: //")
		  container: ${task.container}
		 madeline: 
		  version: \$(echo \$(madeline --version 2>&1) | grep : | sed -e"s/^.*Madeline //; s/PDE : 1.*//")
		  container: ${task.container}
		END_VERSIONS
		"""

	stub:
		"""
		source activate tools
		touch ${group}_madde.INFO
		touch ${ped}.madeline.xml

		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 ped_parser: 
		  version: \$(echo \$(ped_parser --version 2>&1) | sed -e "s/^.*ped_parser version: //")
		  container: ${task.container}
		 madeline: 
		  version: \$(echo \$(madeline --version 2>&1) | grep : | sed -e"s/^.*Madeline //; s/PDE : 1.*//")
		  container: ${task.container}
		END_VERSIONS
		"""
}

process freebayes {
	publishDir "/${OUTDIR}/versions", mode: 'copy' , overwrite: 'true', pattern: '*versions.yml'
	cpus 1
	time '2h'
	container '/fs1/resources/containers/twistmyeloid_2020-06-17.sif'
	scratch true
	stageInMode 'copy'
	stageOutMode 'copy'

	when: 
		params.onco || params.assay == "exome"

	input:
		set group, id, file(bam), file(bai) from bam_freebayes.mix(bam_freebayes_choice)

	output:
		set group, file("${id}.pathfreebayes.lines") into freebayes_concat
		path "*versions.yml"

	script:
		if (params.onco) {
			"""
			freebayes -f $genome_file --pooled-continuous --pooled-discrete -t $params.intersect_bed --min-repeat-entropy 1 -F 0.03 $bam > ${id}.freebayes.vcf
			vcfbreakmulti ${id}.freebayes.vcf > ${id}.freebayes.multibreak.vcf
			bcftools norm -m-both -c w -O v -f $genome_file -o ${id}.freebayes.multibreak.norm.vcf ${id}.freebayes.multibreak.vcf
			vcfanno_linux64 -lua /fs1/resources/ref/hg19/bed/scout/sv_tracks/silly.lua $params.vcfanno ${id}.freebayes.multibreak.norm.vcf > ${id}.freebayes.multibreak.norm.anno.vcf
			grep ^# ${id}.freebayes.multibreak.norm.anno.vcf > ${id}.freebayes.multibreak.norm.anno.path.vcf
			grep -v ^# ${id}.freebayes.multibreak.norm.anno.vcf | grep -i pathogenic > ${id}.freebayes.multibreak.norm.anno.path.vcf2
			cat ${id}.freebayes.multibreak.norm.anno.path.vcf ${id}.freebayes.multibreak.norm.anno.path.vcf2 > ${id}.freebayes.multibreak.norm.anno.path.vcf3
			filter_freebayes.pl ${id}.freebayes.multibreak.norm.anno.path.vcf3 > ${id}.pathfreebayes.lines

			cat <<-END_VERSIONS > ${task.process}_versions.yml
			${task.process}:
			 freebayes: 
			  version: \$(echo \$(freebayes --version 2>&1) | sed 's/version:\s*v//g')
			  container: ${task.container}
			 vcflib: 
			  version: 1.0.9
			  container: ${task.container}
			 bcftools: 
			  version: \$(echo \$(bcftools --version 2>&1) | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
			  container: ${task.container}
			 vcfanno: 
			  version: \$(echo \$(vcfanno_linux64 2>&1) | grep version | cut -f3 -d' ' )
			  container: ${task.container}
			END_VERSIONS
			"""
		}
		else {
			"""
			touch ${id}.pathfreebayes.lines

			cat <<-END_VERSIONS > ${task.process}_versions.yml
			${task.process}:
			 freebayes: 
			  version: \$(echo \$(freebayes --version 2>&1) | sed 's/version:\s*v//g')
			  container: ${task.container}
			 vcflib: 
			  version: 1.0.9
			  container: ${task.container}
			 bcftools: 
			  version: \$(echo \$(bcftools --version 2>&1) | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
			  container: ${task.container}
			 vcfanno: 
			  version: \$(echo \$(vcfanno_linux64 2>&1) | grep version | cut -f3 -d' ' )
			  container: ${task.container}
			END_VERSIONS
			"""
		}

	stub:
		"""
		touch ${id}.pathfreebayes.lines

		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 freebayes: 
		  version: \$(echo \$(freebayes --version 2>&1) | sed 's/version:\s*v//g')
		  container: ${task.container}
		 vcflib: 
		  version: 1.0.9
		  container: ${task.container}
		 bcftools: 
		  version: \$(echo \$(bcftools --version 2>&1) | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
		  container: ${task.container}
		 vcfanno: 
		  version: \$(echo \$(vcfanno_linux64 2>&1) | grep version | cut -f3 -d' ' )
		  container: ${task.container}
		END_VERSIONS
		"""
}

/////////////// MITOCHONDRIA SNV CALLING ///////////////
///////////////                          ///////////////

// create an MT BAM file
process fetch_MTseqs {
	cpus 2
	memory '10GB'
	time '30m'
	tag "$id"
	publishDir "/${OUTDIR}/bam", mode: 'copy', overwrite: 'true', pattern: '*.bam*'
	publishDir "/${OUTDIR}/versions", mode: 'copy' , overwrite: 'true', pattern: '*versions.yml'

	when:
		params.antype == "wgs"

	input:
		set group, id, file(bam), file(bai) from bam_mito.mix(bam_mito_choice)

    output:
        set group, id, file ("${id}_mito.bam"), file("${id}_mito.bam.bai") into mutserve_bam, eklipse_bam, qc_mito_bam
		set group, file("${group}_mtbam.INFO") into mtBAM_INFO
		path "*versions.yml"

	script:
		"""
		sambamba view -f bam $bam M > ${id}_mito.bam
		samtools index -b ${id}_mito.bam
		echo "mtBAM	$id	/access/${params.subdir}/bam/${id}_mito.bam" > ${group}_mtbam.INFO

		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 Sambamba: 
		  version: \$(echo \$(sambamba --version 2>&1) | awk '{print \$2}' )
		  container: ${task.container}
		 SAMtools: 
		  version: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
		  container: ${task.container}
		END_VERSIONS
		"""

	stub:
		"""
		touch ${id}_mito.bam
		touch ${id}_mito.bam.bai
		touch ${group}_mtbam.INFO

		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 Sambamba: 
		  version: \$(echo \$(sambamba --version 2>&1) | awk '{print \$2}' )
		  container: ${task.container}
		 SAMtools: 
		  version: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
		  container: ${task.container}
		END_VERSIONS
		"""
}

process sentieon_mitochondrial_qc {

    // Fetch mitochondrial coverage statistics
    // Calculate mean_coverage and pct_above_500x
    
    cpus 52
    memory '30 GB'
	tag "$id"
	time '2h'
	scratch true
	stageInMode 'copy'
	stageOutMode 'copy'
	container = "/fs1/resources/containers/sentieon_202112.sif"

	when:
	    params.antype == "wgs"
    
	input:
        set group, id, file(bam), file(bai) from qc_mito_bam

	output:
    	set group, id, file("${id}_mito_coverage.tsv") into qc_mito

	script:	
		"""
		sentieon driver \\
			-r $genome_file \\
			-t ${task.cpus} \\
			-i $bam \\
			--algo CoverageMetrics \\
			--cov_thresh 500 \\
			mt_cov_metrics.txt

		head -1 mt_cov_metrics.txt.sample_interval_summary > "${id}_mito_coverage.tsv"
		grep "^M" mt_cov_metrics.txt.sample_interval_summary >> "${id}_mito_coverage.tsv"
		"""

	stub:
		"""
		touch ${id}_mito_coverage.tsv
		"""
}

process build_mitochondrial_qc_json {
   
    cpus 4
    tag "$id"
    time "10m"

    input:
        set group, id, file(mito_qc_file) from qc_mito

    output:
        set group, id, file("${id}_mito_qc.json") into qc_mito_json
    
	script:
		"""
		mito_tsv_to_json.py ${mito_qc_file} ${id} > "${id}_mito_qc.json"
		"""
	stub:
		"""
		touch ${id}_mito_qc.json
		"""
}

    
// gatk FilterMutectCalls in future if FPs overwhelms tord/sofie/carro
process run_mutect2 {
	cpus 4
	memory '50 GB'
	time '1h'
	tag "$group"
	publishDir "/${OUTDIR}/vcf", mode: 'copy', overwrite: 'true', pattern: '*.vcf'
	publishDir "/${OUTDIR}/versions", mode: 'copy' , overwrite: 'true', pattern: '*versions.yml'

	when:
		!params.onco
	
	input:
		set group, id, file(bam), file(bai) from mutserve_bam.groupTuple()

	output:
		set group, id, file("${group}.mutect2.vcf") into ms_vcfs_1, ms_vcfs_2
		path "*versions.yml"

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

		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 GATK: 
		  version: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$// ; s/-SNAPSHOT//')
		  container: ${task.container}
		END_VERSIONS
		"""

	stub:
		"""
		source activate gatk4-env
		touch ${group}.mutect2.vcf

		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 GATK: 
		  version: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$// ; s/-SNAPSHOT//')
		  container: ${task.container}
		END_VERSIONS
		"""
}

// split and left-align variants
process split_normalize_mito {
	publishDir "/${OUTDIR}/versions", mode: 'copy' , overwrite: 'true', pattern: '*versions.yml'
	cpus 1
	memory '1GB'
	time '10m'

	input:
		set group, id, file(ms_vcf) from ms_vcfs_1
		set g2, id2, sex, type from meta_mutect2.groupTuple()

	output:
		set group, file("${group}.mutect2.breakmulti.filtered5p.0genotyped.proband.vcf") into adj_vcfs
		path "*versions.yml"

	script:
		proband_idx = type.findIndexOf{ it == "proband" }

		"""
		grep -vP "^M\\s+955" $ms_vcf > ${ms_vcf}.fix
		bcftools norm -m-both -o ${ms_vcf}.breakmulti ${ms_vcf}.fix
		bcftools sort ${ms_vcf}.breakmulti | bgzip > ${ms_vcf}.breakmulti.fix
		tabix -p vcf ${ms_vcf}.breakmulti.fix
		bcftools norm -f $params.rCRS_fasta -o ${ms_vcf.baseName}.adjusted.vcf ${ms_vcf}.breakmulti.fix
		bcftools view -i 'FMT/AF[*]>0.05' ${ms_vcf.baseName}.adjusted.vcf -o ${group}.mutect2.breakmulti.filtered5p.vcf
		bcftools filter -S 0 --exclude 'FMT/AF[*]<0.05' ${group}.mutect2.breakmulti.filtered5p.vcf -o ${group}.mutect2.breakmulti.filtered5p.0genotyped.vcf
		filter_mutect2_mito.pl ${group}.mutect2.breakmulti.filtered5p.0genotyped.vcf ${id2[proband_idx]} > ${group}.mutect2.breakmulti.filtered5p.0genotyped.proband.vcf

		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 vcflib: 
		  version: 1.0.9
		  container: ${task.container}
		 bcftools: 
		  version: \$(echo \$(bcftools --version 2>&1) | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
		  container: ${task.container}
		 tabix: 
		  version: \$(echo \$(tabix --version 2>&1) | sed 's/^.*(htslib) // ; s/ Copyright.*//')
		  container: ${task.container}
		END_VERSIONS
		"""

	stub:
		"""
		touch ${group}.mutect2.breakmulti.filtered5p.0genotyped.proband.vcf

		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 vcflib: 
		  version: 1.0.9
		  container: ${task.container}
		 bcftools: 
		  version: \$(echo \$(bcftools --version 2>&1) | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
		  container: ${task.container}
		 tabix: 
		  version: \$(echo \$(tabix --version 2>&1) | sed 's/^.*(htslib) // ; s/ Copyright.*//')
		  container: ${task.container}
		END_VERSIONS
		"""
}

// use python tool HmtNote for annotating vcf
// future merging with diploid genome does not approve spaces in info-string
process run_hmtnote {
	publishDir "/${OUTDIR}/versions", mode: 'copy' , overwrite: 'true', pattern: '*versions.yml'
	cpus 1
	memory '5GB'
	time '15m'


	input:
		set group, file(adj_vcf) from adj_vcfs

	output:
		set group, file("${group}.fixinfo.vcf") into mito_diplod_vep
		path "*versions.yml"

	script:
		"""
		source activate tools
		hmtnote annotate ${adj_vcf} ${group}.hmtnote --offline
		grep ^# ${group}.hmtnote > ${group}.fixinfo.vcf
		grep -v ^# ${group}.hmtnote | sed 's/ /_/g' >> ${group}.fixinfo.vcf

		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 hmtnote: 
		  version: \$(echo \$(hmtnote --version 2>&1) | sed 's/^.*hmtnote, version //; s/Using.*\$//' )
		  container: ${task.container}
		END_VERSIONS
		"""

	stub:
		"""
		source activate tools
		touch ${group}.fixinfo.vcf

		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 hmtnote: 
		  version: \$(echo \$(hmtnote --version 2>&1) | sed 's/^.*hmtnote, version //; s/Using.*\$//' )
		  container: ${task.container}
		END_VERSIONS
		"""
}

// run haplogrep 2 on resulting vcf
process run_haplogrep {
	time '10m'
	memory '30 GB'
	cpus '2'
	publishDir "/${OUTDIR}/plots/mito", mode: 'copy', overwrite: 'true', pattern: "*.png"
	publishDir "/${OUTDIR}/versions", mode: 'copy' , overwrite: 'true', pattern: '*versions.yml'

	input:
		set group, id, file(ms_vcf) from ms_vcfs_2

	output:
		file("${group}.haplogrep.png")
		set group, file("${group}_haplo.INFO") into haplogrep_INFO
		path "*versions.yml"

	shell:
		'''
		for sample in `bcftools query -l !{ms_vcf}`; do 
			bcftools view -c1 -Oz -s $sample -o $sample.vcf.gz !{ms_vcf}
			java  -Xmx16G -Xms16G -jar /opt/bin/haplogrep.jar classify \
			--in $sample.vcf.gz\
			--out $sample.hg2.vcf \
			--format vcf \
			--lineage 1
			dot $sample.hg2.vcf.dot -Tps2 > $sample.hg2.vcf.ps
			gs -dSAFER -dBATCH -dNOPAUSE -sDEVICE=png16m -dGraphicsAlphaBits=4 -r1200 -dDownScaleFactor=3 -sOutputFile=${sample}.hg2.vcf.png ${sample}.hg2.vcf.ps
		done
		montage -mode concatenate -tile 3x1 *.png !{group}.haplogrep.png
		echo "IMG haplogrep !{params.accessdir}/plots/mito/!{group}.haplogrep.png" > !{group}_haplo.INFO
		
		cat <<-END_VERSIONS > !{task.process}_versions.yml
		!{task.process}:
		 haplogrep: 
		  version: $(echo $(java -jar /opt/bin/haplogrep.jar classify 2>&1) | sed "s/.*Classify v// ; s/ .*//")
		  container: !{task.container}
		 montage: 
		  version: $(echo $(gm -version 2>&1) | head -1 | sed -e "s/GraphicsMagick //" | cut -d" " -f1 )
		  container: !{task.container}
		END_VERSIONS
		'''

	stub:
		"""
		touch ${group}.haplogrep.png
		touch ${group}_haplo.INFO

		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 haplogrep: 
		  version: \$(echo \$(java -jar /opt/bin/haplogrep.jar classify 2>&1) | sed "s/htt.*Classify v// ; s/ .*//")
		  container: ${task.container}
		 montage: 
		  version: \$(echo \$(gm -version 2>&1) | head -1 | sed -e "s/GraphicsMagick //" | cut -d" " -f1 )
		  container: ${task.container}
		END_VERSIONS
		"""
}

// use eKLIPse for detecting mitochondrial deletions
process run_eklipse {
	cpus 2
	memory '10GB'
	time '60m'
	publishDir "/${OUTDIR}/plots/mito", mode: 'copy', overwrite: 'true', pattern: '*.png'
	publishDir "/${OUTDIR}/plots/mito", mode: 'copy', overwrite: 'true', pattern: '*.hetplasmid_frequency.txt'
	publishDir "/${OUTDIR}/versions", mode: 'copy' , overwrite: 'true', pattern: '*versions.yml'

	input:
		set group, id, file(bam), file(bai) from eklipse_bam

	output:
		set file("*.png"), file("${id}.hetplasmid_frequency.txt")
		set group, file("${id}_eklipse.INFO") into eklipse_INFO
		path "*versions.yml"

	script:
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
		echo "IMG eklipse ${params.accessdir}/plots/mito/${id}_eklipse.png" > ${id}_eklipse.INFO

		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 eKLIPse: 
		  version: 1-8
		  container: ${task.container}
		END_VERSIONS
		"""

	stub:
		"""
		source activate htslib10
		touch ${id}.hetplasmid_frequency.txt
		touch ${id}_eklipse.INFO
		touch ${id}.png

		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 eKLIPse: 
		  version: 1-8
		  container: ${task.container}
		END_VERSIONS
		"""
}

//eklipseM_INFO.collectFile(name: "eklipse.INFO").set{ eklipse_INFO }

// Splitting & normalizing variants, merging with Freebayes/Mutect2, intersecting against exome/clinvar introns
process split_normalize {
	cpus 1
	publishDir "/${OUTDIR}/vcf", mode: 'copy', overwrite: 'true', pattern: '*.vcf'
	publishDir "/${OUTDIR}/versions", mode: 'copy' , overwrite: 'true', pattern: '*versions.yml'
	tag "$group"
	memory '10 GB'
	time '1h'
	scratch true
	stageInMode 'copy'
	stageOutMode 'copy'

	when:
		params.annotate

	input:
		set group, id, file(vcf), file(idx), file(vcfconcat) from combined_vcf.join(mito_diplod_vep.mix(freebayes_concat))

	output:
		set group, file("${group}.norm.uniq.DPAF.vcf") into split_norm, vcf_gnomad
		set group, id, file("${group}.intersected.vcf"), file("${group}.multibreak.vcf") into split_vep, split_cadd, vcf_cnvkit
		path "*versions.yml"

	script:
	id = id[0]
	// rename M to MT because genmod does not recognize M
	if(params.onco) {
		"""
		cat $vcf $vcfconcat > ${id}.concat.freebayes.vcf
		vcfbreakmulti ${id}.concat.freebayes.vcf > ${group}.multibreak.vcf
		bcftools norm -m-both -c w -O v -f $genome_file -o ${group}.norm.vcf ${group}.multibreak.vcf
		bcftools sort ${group}.norm.vcf | vcfuniq > ${group}.norm.uniq.vcf
		wgs_DPAF_filter.pl ${group}.norm.uniq.vcf > ${group}.norm.uniq.DPAF.vcf
		bedtools intersect \
			-a ${group}.norm.uniq.DPAF.vcf \\
			-b $params.intersect_bed \\
			-u -header > ${group}.intersected.vcf
		
		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 vcflib: 
		  version: 1.0.9
		  container: ${task.container}
		 bcftools: 
		  version: \$(echo \$(bcftools --version 2>&1) | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
		  container: ${task.container}
		 bedtools: 
		  version: \$(echo \$(bedtools --version 2>&1) | sed -e "s/^.*bedtools v//" )
		  container: ${task.container}
		 MergeVcfs: 
		  version: \$(echo \$(java -jar /opt/conda/envs/CMD-WGS/share/picard-2.21.2-1/picard.jar MergeVcfs --version 2>&1))
		  container: ${task.container}
		END_VERSIONS
		"""
	}

	else {
		"""
		vcfbreakmulti ${vcf} > ${group}.multibreak.vcf
		bcftools norm -m-both -c w -O v -f $genome_file -o ${group}.norm.vcf ${group}.multibreak.vcf
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

		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 vcflib: 
		  version: 1.0.9
		  container: ${task.container}
		 bcftools: 
		  version: \$(echo \$(bcftools --version 2>&1) | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
		  container: ${task.container}
		 bedtools: 
		  version: \$(echo \$(bedtools --version 2>&1) | sed -e "s/^.*bedtools v//" )
		  container: ${task.container}
		 MergeVcfs: 
		  version: \$(echo \$(java -jar /opt/conda/envs/CMD-WGS/share/picard-2.21.2-1/picard.jar MergeVcfs --version 2>&1))
		  container: ${task.container}
		END_VERSIONS
		"""
	}

	stub:
		"""
		touch ${group}.norm.uniq.DPAF.vcf
		touch ${group}.intersected.vcf
		touch ${group}.multibreak.vcf

		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 vcflib: 
		  version: 1.0.9
		  container: ${task.container}
		 bcftools: 
		  version: \$(echo \$(bcftools --version 2>&1) | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
		  container: ${task.container}
		 bedtools: 
		  version: \$(echo \$(bedtools --version 2>&1) | sed -e "s/^.*bedtools v//" )
		  container: ${task.container}
		 MergeVcfs: 
		  version: \$(echo \$(java -jar /opt/conda/envs/CMD-WGS/share/picard-2.21.2-1/picard.jar MergeVcfs --version 2>&1))
		  container: ${task.container}
		END_VERSIONS
		"""

}

/////////////// Collect QC into single file ///////////////

process merge_qc_json {   
    cpus 1
    errorStrategy 'retry'
    maxErrors 5
    publishDir "${OUTDIR}/qc", mode: 'copy' , overwrite: 'true', pattern: '*.QC'
    tag "$id"
    time '10m'

    input:
        set group, id, file(qc) from qc_cdm.mix(qc_mito_json).groupTuple(by: [0,1])

    output:
        set id, file("${id}.QC") into qc_cdm_merged

    script:
        qc_json_files = qc.join(' ')

    """
    merge_json_files.py ${qc_json_files} > ${id}.QC
    """
}   
    
// Load QC data into CDM (via middleman)
process qc_to_cdm {
	cpus 1
	errorStrategy 'retry'
	maxErrors 5
	publishDir "${CRONDIR}/qc", mode: 'copy' , overwrite: 'true'
	tag "$id"
	time '10m'

	when:
		!params.noupload
	
	input:
		set id, file(qc), diagnosis, r1, r2 from qc_cdm_merged.join(qc_extra)

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
    
    
process annotate_vep {
	publishDir "/${OUTDIR}/versions", mode: 'copy' , overwrite: 'true', pattern: '*versions.yml'
	container = '/fs1/resources/containers/ensembl-vep_release_103.sif'
	cpus 54
	tag "$group"
	memory '150 GB'
	time '5h'
	scratch true
	stageInMode 'copy'
	stageOutMode 'copy'

	input:
		set group, id, file(vcf), idx from split_vep.mix(annotate_only_vep)

	output:
		set group, file("${group}.vep.vcf") into vep
		path "*versions.yml"

	script:
		"""
		vep \\
			-i ${vcf} \\
			-o ${group}.vep.vcf \\
			--offline \\
			--everything \\
			--merged \\
			--vcf \\
			--no_stats \\
			--synonyms $params.SYNONYMS \\
			--fork ${task.cpus} \\
			--force_overwrite \\
			--fasta $params.VEP_FASTA \\
			--dir_cache $params.VEP_CACHE \\
			--dir_plugins $params.VEP_CACHE/Plugins \\
			--distance 200 \\
			-cache \\
			--plugin CADD,$params.CADD \\
			--plugin LoFtool \\
			--plugin MaxEntScan,$params.MAXENTSCAN,SWA,NCSS \\
			--plugin dbNSFP,/fs1/resources/ref/hg38/annotation_dbs/dbnsfp/dbNSFP4.3a_grch38.gz,transcript_match=1,REVEL_score,REVEL_rankscore \\
			-custom $params.GNOMAD_EXOMES,gnomADe,vcf,exact,0,AF_popmax,AF,popmax \\
			-custom $params.GNOMAD_GENOMES,gnomADg,vcf,exact,0,AF_popmax,AF,popmax \\
			-custom $params.GNOMAD_MT,gnomAD_mt,vcf,exact,0,AF_hom,AF_het \\
			-custom $params.PHYLOP \\
			-custom $params.PHASTCONS

		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 VEP: 
		  version: \$( echo \$(vep --help 2>&1) | sed 's/^.*Versions:.*ensembl-vep : // ; s/ .*\$//')
		  container: ${task.container}
		END_VERSIONS
		"""

	stub:
		"""
		touch ${group}.vep.vcf

		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 VEP: 
		  version: \$( echo \$(vep --help 2>&1) | sed 's/^.*Versions:.*ensembl-vep : // ; s/ .*\$//')
		  container: ${task.container}
		END_VERSIONS
		"""

}


// gene, clinvar, loqusdb, enigma(onco)
process vcfanno {
	publishDir "/${OUTDIR}/versions", mode: 'copy' , overwrite: 'true', pattern: '*versions.yml'
	cpus params.cpu_some
	memory '32GB'
	time '20m'
	errorStrategy 'retry'
	maxErrors 5

	input:
		set group, file(vcf) from vep

	output:
		set group, file("${group}.clinvar.loqusdb.gene.vcf") into vcfanno_vcf
		path "*versions.yml"

	script:
		"""
		vcfanno_linux64 -lua /fs1/resources/ref/hg19/bed/scout/sv_tracks/silly.lua $params.vcfanno $vcf > ${group}.clinvar.loqusdb.gene.vcf

		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 vcfanno: 
		  version: \$(echo \$(vcfanno 2>&1) | grep version | cut -f3 -d' ' )
		  container: ${task.container}
		END_VERSIONS
		"""

	stub:
		"""
		touch ${group}.clinvar.loqusdb.gene.vcf

		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 vcfanno: 
		  version: \$(echo \$(vcfanno 2>&1) | grep version | cut -f3 -d' ' )
		  container: ${task.container}
		END_VERSIONS
		"""
}

// Extracting most severe consequence: 
// Modifying annotations by VEP-plugins, and adding to info-field: 
// Modifying CLNSIG field to allow it to be used by genmod score properly:
process modify_vcf {
	cpus 1
	tag "$group"
	memory '1 GB'
	time '10m'

	input:
		set group, file(vcf) from vcfanno_vcf

	output:
		set group, file("${group}.mod.vcf") into mod_vcf

	script:
		"""
		modify_vcf_scout.pl $vcf > ${group}.mod.vcf
		"""

	stub:
		"""
		touch ${group}.mod.vcf
		"""
} 


// Marking splice INDELs: 
process mark_splice {
	cpus 1
	tag "$group"
	memory '1 GB'
	time '20m'

	input:
		set group, file(vcf) from mod_vcf

	output:
		set group, file("${group}.marksplice.vcf") into splice_marked

	script:
		"""
		/opt/bin/mark_spliceindels.pl $vcf > ${group}.marksplice.vcf
		"""

	stub:
		"""
		touch ${group}.marksplice.vcf
		"""
}

// Extract all INDELs from VCF for CADD annotation
process extract_indels_for_cadd {
	cpus 1
	tag "$group"
	memory '1 GB'
	time '5m'
	publishDir "/${OUTDIR}/versions", mode: 'copy' , overwrite: 'true', pattern: '*versions.yml'

	input:
		set group, id, file(vcf), idx from split_cadd.mix(annotate_only_cadd)
	
	output:
		set group, file("${group}.only_indels.vcf") into indel_cadd_vep
		path "*versions.yml"

	script:
		"""
		bcftools view $vcf -V snps -o ${group}.only_indels.vcf

		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 bcftools: 
		  version: \$(echo \$(bcftools --version 2>&1) | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
		  container: ${task.container}
		END_VERSIONS
		"""

	stub:
		"""
		touch ${group}.only_indels.vcf

		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 bcftools: 
		  version: \$(echo \$(bcftools --version 2>&1) | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
		  container: ${task.container}
		END_VERSIONS
		""" 
}

// Annotate Indels with VEP+Gnomad genomes. Filter variants below threshold
process indel_vep {
	publishDir "/${OUTDIR}/versions", mode: 'copy' , overwrite: 'true', pattern: '*versions.yml'
	cpus 5
	container = '/fs1/resources/containers/ensembl-vep_release_103.sif'
	tag "$group"
	memory '10 GB'
	time '3h'

	input:
		set group, file(vcf) from indel_cadd_vep

	output:
		set group, file("${group}.only_indels.vep.filtered.vcf") into indel_cadd_vcf
		path "*versions.yml"

	script:
		"""
		vep \\
			-i $vcf \\
			-o ${group}.only_indels.vep.vcf \\
			--offline \\
			--cache \\
			--merged \\
			--vcf \\
			--synonyms $params.SYNONYMS \\
			--fasta $params.VEP_FASTA \\
			-custom $params.GNOMAD_GENOMES,gnomADg,vcf,exact,0,AF \\
			-custom $params.GNOMAD_MT,gnomAD_mt,vcf,exact,0,AF_hom,AF_het \\
			--dir_cache $params.VEP_CACHE \\
			--force_overwrite \\
			--no_stats \\
			--fork ${task.cpus}
		filter_indels.pl ${group}.only_indels.vep.vcf > ${group}.only_indels.vep.filtered.vcf

		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 VEP: 
		  version: \$( echo \$(vep --help 2>&1) | sed 's/^.*Versions:.*ensembl-vep : // ; s/ .*\$//')
		  container: ${task.container}
		END_VERSIONS
		"""

	stub:
		"""
		touch ${group}.only_indels.vep.filtered.vcf

		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 VEP: 
		  version: \$( echo \$(vep --help 2>&1) | sed 's/^.*Versions:.*ensembl-vep : // ; s/ .*\$//')
		  container: ${task.container}
		END_VERSIONS
		"""
}

// Calculate CADD scores for all indels
process calculate_indel_cadd {
	publishDir "/${OUTDIR}/versions", mode: 'copy' , overwrite: 'true', pattern: '*versions.yml'
	cpus 2
	container = '/fs1/resources/containers/cadd_v1.6.sif'
	scratch true
	stageInMode 'copy'
	stageOutMode 'copy'
	tag "$group"
	memory '15 GB'
	time '3h'

	input:
		set group, file(vcf) from indel_cadd_vcf

	output:
		set group, file("${group}.indel_cadd.gz") into indel_cadd
		path "*versions.yml"

	script:
		"""
		/CADD-scripts/CADD.sh -c ${task.cpus} -g GRCh38 -o ${group}.indel_cadd.gz $vcf

		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 CADD: 
		  version: \$(echo \$(/CADD-scripts/CADD.sh -v 2>&1) | sed -e "s/^.*CADD-v// ; s/ (c).*//")
		  container: ${task.container}
		END_VERSIONS
		"""

	stub:
		"""
		touch ${group}.indel_cadd.gz

		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 CADD: 
		  version: \$(echo \$(/CADD-scripts/CADD.sh -v 2>&1) | sed -e "s/^.*CADD-v// ; s/ (c).*//")
		  container: ${task.container}
		END_VERSIONS
		"""
}

// Add the calculated indel CADDs to the vcf
process add_cadd_scores_to_vcf {
	publishDir "/${OUTDIR}/versions", mode: 'copy' , overwrite: 'true', pattern: '*versions.yml'
	cpus 4
	tag "$group"
	memory '1 GB'
	time '5m'
	container = '/fs1/resources/containers/genmod.sif'

	input: 
		set group, file(vcf), file(cadd_scores) from splice_marked.join(indel_cadd)

	output:
		set group, file("${group}.cadd.vcf") into ma_vcf, fa_vcf, base_vcf
		path "*versions.yml"

	script:
		"""
		gunzip -c $cadd_scores > cadd
		bgzip -@ ${task.cpus} cadd
		tabix -p vcf cadd.gz
		genmod annotate --cadd-file cadd.gz $vcf > ${group}.cadd.vcf

		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 gunzip: 
		  version: \$(echo \$(gunzip --version 2>&1) | sed -e "s/^.*(gzip) // ; s/Copyright.*//")
		  container: ${task.container}
		 bgzip: 
		  version: \$(echo \$(bgzip --version 2>&1) | sed 's/^.*(htslib) // ; s/ Copyright.*//')
		  container: ${task.container}
		 tabix: 
		  version: \$(echo \$(tabix --version 2>&1) | sed 's/^.*(htslib) // ; s/ Copyright.*//')
		  container: ${task.container}
		 genmod: 
		  version: \$(echo \$(genmod --version 2>&1) | sed -e "s/^.*genmod version: //")
		  container: ${task.container}
		END_VERSIONS
		"""

	stub:
		"""
		touch ${group}.cadd.vcf

		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 gunzip: 
		  version: \$(echo \$(gunzip --version 2>&1) | sed -e "s/^.*(gzip) // ; s/Copyright.*//")
		  container: ${task.container}
		 bgzip: 
		  version: \$(echo \$(bgzip --version 2>&1) | sed 's/^.*(htslib) // ; s/ Copyright.*//')
		  container: ${task.container}
		 tabix: 
		  version: \$(echo \$(tabix --version 2>&1) | sed 's/^.*(htslib) // ; s/ Copyright.*//')
		  container: ${task.container}
		 genmod: 
		  version: \$(echo \$(genmod --version 2>&1) | sed -e "s/^.*genmod version: //")
		  container: ${task.container}
		END_VERSIONS
		"""
}


// # Annotating variant inheritance models:
process inher_models {
	publishDir "/${OUTDIR}/versions", mode: 'copy' , overwrite: 'true', pattern: '*versions.yml'
	cpus 6
	memory '64 GB'
	tag "$group"
	time '10m'
	scratch true
	stageInMode 'copy'
	stageOutMode 'copy'
	container = '/fs1/resources/containers/genmod.sif'

	input:
		set group, file(vcf), type, file(ped) from base_vcf.mix(ma_vcf, fa_vcf).join(ped_inher.mix(ped_inher_ma,ped_inher_fa)).view()

	output:
		set group, type, file("${group}.models.vcf") into inhermod
		path "*versions.yml"

	script:
		"""
		genmod models $vcf -p ${task.cpus} -f $ped > ${group}.models.vcf

		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 genmod: 
		  version: \$(echo \$(genmod --version 2>&1) | sed -e "s/^.*genmod version: //")
		  container: ${task.container}
		END_VERSIONS
		"""

	stub:
		"""
		touch ${group}.models.vcf

		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 genmod: 
		  version: \$(echo \$(genmod --version 2>&1) | sed -e "s/^.*genmod version: //")
		  container: ${task.container}
		END_VERSIONS
		"""
}


// Scoring variants: 
// Adjusting compound scores: 
// Sorting VCF according to score: 
process genmodscore {
	publishDir "/${OUTDIR}/versions", mode: 'copy' , overwrite: 'true', pattern: '*versions.yml'
	cpus 2
	tag "$group"
	memory '10 GB'
	time '30m'
	container = '/fs1/resources/containers/genmod.sif'

	input:
		set group, type, file(vcf) from inhermod

	output:
		set group, type, file("${group_score}.scored.vcf") into scored_vcf
		path "*versions.yml"

	script:
		group_score = group
		if ( type == "ma" || type == "fa") {
			group_score = group + "_" + type
		}

		if ( mode == "family" && params.antype == "wgs" ) {
			"""
			genmod score -i $group_score -c $params.rank_model -r $vcf -o ${group_score}.score1.vcf
			genmod compound ${group_score}.score1.vcf > ${group_score}.score2.vcf
			sed 's/RankScore=${group}:/RankScore=${group_score}:/g' -i ${group_score}.score2.vcf
			genmod sort -p -f $group_score ${group_score}.score2.vcf -o ${group_score}.scored.vcf

			cat <<-END_VERSIONS > ${task.process}_versions.yml
			${task.process}:
			 genmod: 
			  version: \$(echo \$(genmod --version 2>&1) | sed -e "s/^.*genmod version: //")
			  container: ${task.container}
			END_VERSIONS
			"""
		}
		else {
			"""
			genmod score -i $group_score -c $params.rank_model_s -r $vcf -o ${group_score}.score1.vcf
			genmod sort -p -f $group_score ${group_score}.score1.vcf -o ${group_score}.scored.vcf

			cat <<-END_VERSIONS > ${task.process}_versions.yml
			${task.process}:
			 genmod: 
			  version: \$(echo \$(genmod --version 2>&1) | sed -e "s/^.*genmod version: //")
			  container: ${task.container}
			END_VERSIONS
			"""
		}

	stub:
		group_score = group
		"""
		touch ${group_score}.scored.vcf

		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 genmod: 
		  version: \$(echo \$(genmod --version 2>&1) | sed -e "s/^.*genmod version: //")
		  container: ${task.container}
		END_VERSIONS
		"""
}

// Bgzipping and indexing VCF: 
process vcf_completion {
	cpus 16
	publishDir "/${OUTDIR}/vcf", mode: 'copy', overwrite: 'true', pattern: '*.vcf.gz*'
	publishDir "/${OUTDIR}/versions", mode: 'copy' , overwrite: 'true', pattern: '*versions.yml'
	tag "$group"
	time '1h'

	input:
		set group, type, file(vcf) from scored_vcf

	output:
		set group, type, file("${group_score}.scored.vcf.gz"), file("${group_score}.scored.vcf.gz.tbi") into vcf_peddy, snv_sv_vcf,snv_sv_vcf_ma,snv_sv_vcf_fa, vcf_loqus
		set group, file("${group}_snv.INFO") into snv_INFO
		path "*versions.yml"

	script:
		group_score = group
		if ( type == "ma" || type == "fa") {
			group_score = group + "_" + type
		}

		"""
		sed 's/^MT/M/' -i $vcf
		sed 's/ID=MT,length/ID=M,length/' -i $vcf
		bgzip -@ ${task.cpus} $vcf -f
		tabix ${vcf}.gz -f
		echo "SNV	$type	${params.accessdir}/vcf/${group_score}.scored.vcf.gz" > ${group}_snv.INFO

		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 bgzip: 
		  version: \$(echo \$(bgzip --version 2>&1) | sed 's/^.*(htslib) // ; s/ Copyright.*//')
		  container: ${task.container}
		 tabix: 
		  version: \$(echo \$(tabix --version 2>&1) | sed 's/^.*(htslib) // ; s/ Copyright.*//')
		  container: ${task.container}
		END_VERSIONS
		"""

	stub:
		group_score = group
		"""
		touch ${group_score}.scored.vcf.gz
		touch ${group_score}.scored.vcf.gz.tbi
		touch ${group}_snv.INFO

		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 bgzip: 
		  version: \$(echo \$(bgzip --version 2>&1) | sed 's/^.*(htslib) // ; s/ Copyright.*//')
		  container: ${task.container}
		 tabix: 
		  version: \$(echo \$(tabix --version 2>&1) | sed 's/^.*(htslib) // ; s/ Copyright.*//')
		  container: ${task.container}
		END_VERSIONS
		"""
}


// Running PEDDY: 
process peddy {
	publishDir "/${OUTDIR}/ped", mode: 'copy' , overwrite: 'true', pattern: '*.csv'
	publishDir "/${OUTDIR}/ped", mode: 'copy' , overwrite: 'true', pattern: '*.ped'
	publishDir "/${OUTDIR}/versions", mode: 'copy' , overwrite: 'true', pattern: '*versions.yml'
	//container = '/fs1/resources/containers/wgs_20200115.sif'
	cpus 6
	tag "$group"
	time '1h'

	when:
		!params.annotate_only

	input:
		set group, type, file(vcf), file(idx), type, file(ped) from vcf_peddy.join(ped_peddy)

	output:
		set file("${group}.ped_check.csv"),file("${group}.peddy.ped"), file("${group}.sex_check.csv") into peddy_files
		set group, file("${group}_peddy.INFO") into peddy_INFO
		path "*versions.yml"

	script:
		"""
		source activate py3-env
		python -m peddy --sites hg38 -p ${task.cpus} $vcf $ped --prefix $group
		echo "PEDDY	${params.accessdir}/ped/${group}.ped_check.csv,${params.accessdir}/ped/${group}.peddy.ped,${params.accessdir}/ped/${group}.sex_check.csv" > ${group}_peddy.INFO

		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 peddy: 
		  version: \$(echo \$(python -m peddy --version 2>&1) | sed 's/^.*peddy, version //')
		  container: ${task.container}
		END_VERSIONS
		"""

	stub:
		"""
		source activate py3-env
		touch ${group}.ped_check.csv
		touch ${group}.peddy.ped
		touch ${group}.sex_check.csv
		touch ${group}_peddy.INFO

		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 peddy: 
		  version: \$(echo \$(python -m peddy --version 2>&1) | sed 's/^.*peddy, version //')
		  container: ${task.container}
		END_VERSIONS
		"""
}

// Extract all variants (from whole genome) with a gnomAD af > x%
process fastgnomad {
	cpus 2
	memory '32 GB'
	tag "$group"
	publishDir "/${OUTDIR}/vcf", mode: 'copy', overwrite: 'true', pattern: '*.vcf'
	publishDir "/${OUTDIR}/versions", mode: 'copy' , overwrite: 'true', pattern: '*versions.yml'
	time '2h'

	when:
		!params.onco && !params.exome

	input:
		set group, file(vcf) from vcf_gnomad

	output:
		set group, file("${group}.SNPs.vcf") into vcf_upd, vcf_roh, vcf_pod
		path "*versions.yml"

	script:
		"""
		gzip -c $vcf > ${vcf}.gz
		annotate -g $params.FASTGNOMAD_REF -i ${vcf}.gz > ${group}.SNPs.vcf

		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 gzip: 
		  version: \$(echo \$(gzip --version 2>&1) | sed 's/^.*gzip // ; s/Copyright.*//' )
		  container: ${task.container}
		END_VERSIONS
		"""

	stub:
		"""
		touch ${group}.SNPs.vcf

		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 gzip: 
		  version: \$(echo \$(gzip --version 2>&1) | sed 's/^.*gzip // ; s/Copyright.*//' )
		  container: ${task.container}
		END_VERSIONS
		"""
}


// Call UPD regions from SNP vcf
process upd {
	publishDir "/${OUTDIR}/versions", mode: 'copy' , overwrite: 'true', pattern: '*versions.yml'
	tag "$group"
	time '10m'
	memory '1 GB'

	input:
		set gr, file(vcf) from vcf_upd
		set group, id, sex, mother, father, phenotype, diagnosis, type, assay, clarity_sample_id, ffpe, analysis from meta_upd.filter{ item -> item[7] == 'proband' }

	output:
		file("upd.bed") into upd_plot
		set group, file("upd.sites.bed") into upd_table
		path "*versions.yml"

	script:
		if( mode == "family" && trio == true ) {
			"""
			upd --vcf $vcf --proband $id --mother $mother --father $father --af-tag GNOMADAF regions > upd.bed
			upd --vcf $vcf --proband $id --mother $mother --father $father --af-tag GNOMADAF sites > upd.sites.bed

			cat <<-END_VERSIONS > ${task.process}_versions.yml
			${task.process}:
			 upd: 
			  version: \$(echo \$(upd --version 2>&1)
			  container: ${task.container}
			END_VERSIONS
			"""
		}
		else {
			"""
			touch upd.bed
			touch upd.sites.bed

			cat <<-END_VERSIONS > ${task.process}_versions.yml
			${task.process}:
			 upd: 
			  version: \$(echo \$(upd --version 2>&1)
			  container: ${task.container}
			END_VERSIONS
			"""
		}

	stub:
		"""
		touch upd.bed
		touch upd.sites.bed

		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 upd: 
		  version: \$(echo \$(upd --version 2>&1)
		  container: ${task.container}
		END_VERSIONS
		"""
}


process upd_table {
	publishDir "/${OUTDIR}/plots", mode: 'copy' , overwrite: 'true'
	tag "$group"
	time '10m'
	memory '1 GB'

	input:
		set group, file(upd_sites) from upd_table

	output:
		file("${group}.UPDtable.xls")

	when:
		mode == "family" && trio == true

	script:
		"""
		upd_table.pl $upd_sites > ${group}.UPDtable.xls
		"""

	stub:
		"""
		touch ${group}.UPDtable.xls
		"""
}


// Call ROH regions from SNP vcf
process roh {
	publishDir "/${OUTDIR}/versions", mode: 'copy' , overwrite: 'true', pattern: '*versions.yml'
	tag "$group"
	time '10m'
	memory '1 GB'

	input:
		set gr, file(vcf) from vcf_roh

	output:
		set gr, file("roh.txt") into roh_plot
		path "*versions.yml"

	script:
		"""
		bcftools roh --rec-rate 1e-9 --AF-tag GNOMADAF ${vcf} -o roh.txt

		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 bcftools: 
		  version: \$(echo \$(bcftools --version 2>&1) | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
		  container: ${task.container}
		END_VERSIONS
		"""

	stub:
		"""
		touch roh.txt

		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 bcftools: 
		  version: \$(echo \$(bcftools --version 2>&1) | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
		  container: ${task.container}
		END_VERSIONS
		"""
}

// Create coverage profile using GATK
process gatkcov {
	publishDir "/${OUTDIR}/cov", mode: 'copy' , overwrite: 'true', pattern: '*.tsv'
	publishDir "/${OUTDIR}/versions", mode: 'copy' , overwrite: 'true', pattern: '*versions.yml'
	tag "$group"
	cpus 2
	memory '60 GB'
	time '5h'

	input:
		set id, group, file(bam), file(bai), gr, sex, type from cov_bam.mix(cov_bam_choice).join(meta_gatkcov, by:1)

	output:
		set group, id, type, sex, file("${id}.standardizedCR.tsv"), file("${id}.denoisedCR.tsv") into cov_plot, cov_gens
		path "*versions.yml"

	when:
		params.gatkcov

	script:
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
		
		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 GATK: 
		  version: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$// ; s/-SNAPSHOT//')
		  container: ${task.container}
		END_VERSIONS
		"""

	stub:
		"""
		source activate gatk4-env
		touch ${id}.standardizedCR.tsv
		touch ${id}.denoisedCR.tsv

		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 GATK: 
		  version: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$// ; s/-SNAPSHOT//')
		  container: ${task.container}
		END_VERSIONS
		"""
}


// Plot ROH, UPD and coverage in a genomic overview plot
process overview_plot {
	publishDir "/${OUTDIR}/plots", mode: 'copy' , overwrite: 'true', pattern: "*.png"
	tag "$group"
	time '20m'
	memory '5 GB'

	input:
		file(upd) from upd_plot
		set gr, file(roh) from roh_plot
		set group, id, type, sex, file(cov_stand), file(cov_denoised) from cov_plot.groupTuple()


	output:
		file("${group}.genomic_overview.png")
		set group, file("${group}_oplot.INFO") into oplot_INFO

	script:
		proband_idx = type.findIndexOf{ it == "proband" }

	script:
		"""
		genome_plotter.pl --dict $params.GENOMEDICT \\
			--sample ${id[proband_idx]} \\
			--upd $upd \\
			--roh $roh \\
			--sex ${sex[proband_idx]} \\
			--cov ${cov_denoised[proband_idx]} \\
			--out ${group}.genomic_overview.png
		echo "IMG overviewplot	${params.accessdir}/plots/${group}.genomic_overview.png" > ${group}_oplot.INFO 
		"""

	stub:
		"""
		touch ${group}.genomic_overview.png
		touch ${group}_oplot.INFO
		"""
}

process generate_gens_data {
	publishDir "/${OUTDIR}/plot_data", mode: 'copy' , overwrite: 'true', pattern: "*.gz*"
	publishDir "${CRONDIR}/gens", mode: 'copy', overwrite: 'true', pattern: "*.gens"
	tag "$group"
	cpus 1
	time '3h'
	memory '5 GB'

	when:
		!params.onco && !params.exome

	input:
		set id, group, file(gvcf), g, type, sex, file(cov_stand), file(cov_denoise) from gvcf_gens_choice.join(cov_gens, by:[1])

	output:
		set file("${id}.cov.bed.gz"), file("${id}.baf.bed.gz"), file("${id}.cov.bed.gz.tbi"), file("${id}.baf.bed.gz.tbi"), file("${id}.overview.json.gz")
		file("${id}.gens") into gens_middleman

	script:
		"""
		generate_gens_data.pl $cov_stand $gvcf $id $params.GENS_GNOMAD
		echo "gens load sample --sample-id $id --genome-build 38 --baf ${params.gens_accessdir}/${id}.baf.bed.gz --coverage ${params.gens_accessdir}/${id}.cov.bed.gz --overview-json ${params.gens_accessdir}/${id}.overview.json.gz" > ${id}.gens
		"""

	stub:
		"""
		touch ${id}.cov.bed.gz
		touch ${id}.baf.bed.gz
		touch ${id}.cov.bed.gz.tbi
		touch ${id}.baf.bed.gz.tbi
		touch ${id}.overview.json.gz
		touch ${id}.gens
		"""
}

// SV-calling //

// GATK panel+wgs //

process gatk_coverage {
	publishDir "/${OUTDIR}/versions", mode: 'copy' , overwrite: 'true', pattern: '*versions.yml'
	cpus 10
	memory '50GB'
	time '2h'
	container = '/fs1/resources/containers/gatk_4.1.9.0.sif'
	scratch true
	stageInMode 'copy'
	stageOutMode 'copy'
	tag "$id"   

	when:
		params.sv && params.gatkcnv

	input:
		set group, id, file(bam), file(bai) from bam_gatk.mix(bam_gatk_choice)

	output:
		set group, id, file("${id}.tsv") into call_ploidy, call_cnv
		path "*versions.yml"

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

		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 GATK: 
		  version: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$// ; s/-SNAPSHOT//')
		  container: ${task.container}
		END_VERSIONS
		"""

	stub:
		"""
		export THEANO_FLAGS="base_compiledir=."
		export MKL_NUM_THREADS=${task.cpus}
		export OMP_NUM_THREADS=${task.cpus}
		set +u
		source activate gatk
		touch ${id}.tsv

		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 GATK: 
		  version: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$// ; s/-SNAPSHOT//')
		  container: ${task.container}
		END_VERSIONS
		"""
}

process gatk_call_ploidy {
	publishDir "/${OUTDIR}/versions", mode: 'copy' , overwrite: 'true', pattern: '*versions.yml'
	cpus 10
	memory '40GB'
	time '2h'
	container = '/fs1/resources/containers/gatk_4.1.9.0.sif'
	scratch true
	stageInMode 'copy'
	stageOutMode 'copy'
	tag "$id"

	input:
		set group, id, file(tsv) from call_ploidy

	output:
		set group, id, file("ploidy.tar") into ploidy_to_cnvcall, ploidy_to_post
		path "*versions.yml"

	script:
		"""
		export THEANO_FLAGS="base_compiledir=."
		export MKL_NUM_THREADS=${task.cpus}
		export OMP_NUM_THREADS=${task.cpus}
		set +u
		source activate gatk
		gatk --java-options "-Xmx20g" DetermineGermlineContigPloidy \\
			--model $params.ploidymodel \\
			-I $tsv \\
			-O ploidy/ \\
			--output-prefix $group
		tar -cvf ploidy.tar ploidy/

		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 GATK: 
		  version: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$// ; s/-SNAPSHOT//')
		  container: ${task.container}
		END_VERSIONS
		"""

	stub:
		"""
		export THEANO_FLAGS="base_compiledir=."
		export MKL_NUM_THREADS=${task.cpus}
		export OMP_NUM_THREADS=${task.cpus}
		set +u
		source activate gatk
		touch ploidy.tar

		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 GATK: 
		  version: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$// ; s/-SNAPSHOT//')
		  container: ${task.container}
		END_VERSIONS
		"""
}

process gatk_call_cnv {
	publishDir "/${OUTDIR}/versions", mode: 'copy' , overwrite: 'true', pattern: '*versions.yml'
	cpus 8
	memory '45GB'
	time '3h'
	container = '/fs1/resources/containers/gatk_4.1.9.0.sif'
	scratch true
	stageInMode 'copy'
	stageOutMode 'copy'
	tag "$id"

	input:
		set group, id, file(tsv), file(ploidy), i, refpart \
			from call_cnv.join(ploidy_to_cnvcall, by: [0,1]).combine(gatk_ref)

	output:
		set group, id, i, file("${group}_${i}.tar") into postprocessgatk
		path "*versions.yml"

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
		
		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 GATK: 
		  version: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$// ; s/-SNAPSHOT//')
		  container: ${task.container}
		END_VERSIONS
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
		touch ${group}_${i}.tar

		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 GATK: 
		  version: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$// ; s/-SNAPSHOT//')
		  container: ${task.container}
		END_VERSIONS
		"""
}

process postprocessgatk {
	publishDir "/${OUTDIR}/versions", mode: 'copy' , overwrite: 'true', pattern: '*versions.yml'
	cpus 8
	memory '40GB'
	time '3h'
	container = '/fs1/resources/containers/gatk_4.1.9.0.sif'
	
	//scratch true
	// stageInMode 'copy'
	// stageOutMode 'copy'
	publishDir "/${OUTDIR}/sv_vcf/", mode: 'copy', overwrite: 'true'
	tag "$id"


	input:
		set group, id, i, file(tar), file(ploidy), shard_no, shard \
			from postprocessgatk.groupTuple(by: [0,1]).join(ploidy_to_post, by: [0,1]).combine(gatk_postprocess.groupTuple(by: [3]))


	output:
		set group, id, \
			file("genotyped-intervals-${group}-vs-cohort30.vcf.gz"), \
			file("genotyped-segments-${group}-vs-cohort30.vcf.gz"), \
			file("denoised-${group}-vs-cohort30.vcf.gz") into called_gatk
		path "*versions.yml"

	script:
		modelshards = shard.join(' --model-shard-path ') // join each reference shard
		caseshards = []
		for (n = 1; n <= i.size(); n++) { // join each shard(n) that's been called
			tmp = group+'_'+i[n-1]+'/'+group+'_'+i[n-1]+'-calls' 
			caseshards = caseshards + tmp
		}
		caseshards = caseshards.join( ' --calls-shard-path ')
 	
	shell:
		'''
		THEANO_FLAGS="base_compiledir=/fs1/resources/theano"
		for model in !{tar}; do
		tar -xvf $model
		done
		tar -xvf !{ploidy}
		set +u
		source activate gatk
		export MKL_NUM_THREADS=!{task.cpus}
		export OMP_NUM_THREADS=!{task.cpus}
		gatk --java-options "-Xmx25g" PostprocessGermlineCNVCalls \
			--allosomal-contig X --allosomal-contig Y \
			--contig-ploidy-calls ploidy/!{group}-calls/ \
			--sample-index 0 \\
			--output-genotyped-intervals genotyped-intervals-!{group}-vs-cohort30.vcf.gz \
			--output-genotyped-segments genotyped-segments-!{group}-vs-cohort30.vcf.gz \
			--output-denoised-copy-ratios denoised-!{group}-vs-cohort30.vcf.gz \
			--sequence-dictionary !{params.GENOMEDICT} \
			--calls-shard-path !{caseshards} \
			--model-shard-path !{modelshards}

		cat <<-END_VERSIONS > !{task.process}_versions.yml
		!{task.process}:
		 GATK: 
		  version: $(echo $(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*// ; s/-SNAPSHOT//')
		  container: !{task.container}
		END_VERSIONS
		'''

	stub:
		"""
		THEANO_FLAGS="base_compiledir=/fs1/resources/theano"
		set +u
		source activate gatk
		export MKL_NUM_THREADS=!{task.cpus}
		export OMP_NUM_THREADS=!{task.cpus}
		source activate gatk
		touch genotyped-intervals-${group}-vs-cohort30.vcf.gz
		touch genotyped-segments-${group}-vs-cohort30.vcf.gz
		touch denoised-${group}-vs-cohort30.vcf.gz

		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 GATK: 
		  version: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$// ; s/-SNAPSHOT//')
		  container: ${task.container}
		END_VERSIONS
		"""
}

process filter_merge_gatk {
	cpus 1
	tag "$group"
	time '2h'
	memory '1 GB'
	publishDir "/${OUTDIR}/sv_vcf", mode: 'copy', overwrite: 'true'

	input:
		set group, id, file(inter), file(gatk), file(denoised) from called_gatk

	output:
		set group, id, file("${id}.gatk.filtered.merged.vcf") into merged_gatk,merged_gatk_panel

	script:
		"""
		filter_gatk.pl $gatk > ${id}.gatk.filtered.vcf
		mergeGATK.pl ${id}.gatk.filtered.vcf > ${id}.gatk.filtered.merged.vcf
		"""

	
	stub:
		"""
		touch ${id}.gatk.filtered.merged.vcf
		"""
}


process manta {
	cpus = 56
	publishDir "/${OUTDIR}/sv_vcf/", mode: 'copy', overwrite: 'true', pattern: '*.vcf.gz'
	publishDir "/${OUTDIR}/versions", mode: 'copy' , overwrite: 'true', pattern: '*versions.yml'
	tag "$id"
	time '10h'
	memory '150 GB'
	scratch true
	stageInMode 'copy'
	stageOutMode 'copy'

	when:
		params.sv && !params.onco && !params.exome

	input:
		set group, id, file(bam), file(bai) from bam_manta.mix(bam_manta_choice)

	output:
		set group, id, file("${id}.manta.vcf.gz") into called_manta
		path "*versions.yml"

	script:
		bams = bam.join('--bam ')

		"""
		configManta.py --bam $bam --reference $genome_file --runDir .
		python runWorkflow.py -m local -j ${task.cpus}
		mv results/variants/diploidSV.vcf.gz ${id}.manta.vcf.gz
		mv results/variants/diploidSV.vcf.gz.tbi ${id}.manta.vcf.gz.tbi

		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 Manta: 
		  version: \$( configManta.py --version )
		  container: ${task.container}
		END_VERSIONS
		"""

	stub:
		"""
		touch ${id}.manta.vcf.gz

		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 Manta: 
		  version: \$( configManta.py --version )
		  container: ${task.container}
		END_VERSIONS
		"""
}

process manta_panel {
	cpus = 56
	publishDir "/${OUTDIR}/sv_vcf/", mode: 'copy', overwrite: 'true', pattern: '*.vcf.gz'
	publishDir "/${OUTDIR}/versions", mode: 'copy' , overwrite: 'true', pattern: '*versions.yml'
	tag "$id"
	time '1h'
	memory '150 GB'
	scratch true
	stageInMode 'copy'
	stageOutMode 'copy'

	when:
		params.sv && params.onco

	input:
		set group, id, file(bam), file(bai) from bam_manta_panel.mix(bam_mantapanel_choice)

	output:
		set group, id, file("${id}.manta.vcf.gz") into called_manta_panel
		path "*versions.yml"


	script:
		"""
		configManta.py --bam $bam --reference $genome_file --runDir . --exome --callRegions $params.bedgz --generateEvidenceBam
		python runWorkflow.py -m local -j ${task.cpus}
		mv results/variants/diploidSV.vcf.gz ${id}.manta.vcf.gz
		mv results/variants/diploidSV.vcf.gz.tbi ${id}.manta.vcf.gz.tbi

		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 Manta: 
		  version: \$( configManta.py --version )
		  container: ${task.container}
		END_VERSIONS
		"""

	stub:
		"""
		touch ${id}.manta.vcf.gz

		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 Manta: 
		  version: \$( configManta.py --version )
		  container: ${task.container}
		END_VERSIONS
		"""
}

process delly_panel {
	cpus = 5
	publishDir "/${OUTDIR}/sv_vcf/", mode: 'copy', overwrite: 'true', pattern: '*.vcf.gz'
	publishDir "/${OUTDIR}/versions", mode: 'copy' , overwrite: 'true', pattern: '*versions.yml'
	tag "$id"
	time '3h'
	memory '10 GB'
	scratch true
	stageInMode 'copy'
	stageOutMode 'copy'
	cache 'deep'
	
	when:
		params.sv && params.onco && params.delly

	input:
		set group, id, file(bam), file(bai) from bam_delly_panel.mix(bam_dellypanel_choice)

	output:
		set group, id, file("${id}.delly.vcf.gz") into called_delly_panel
		path "*versions.yml"

	script:
		"""
		delly call -g $genome_file -o ${id}.bcf $bam
		bcftools view ${id}.bcf > ${id}.vcf
		filter_delly.pl --vcf ${id}.vcf --bed $params.intersect_bed > ${id}.delly.vcf
		bgzip -c ${id}.delly.vcf > ${id}.delly.vcf.gz

		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 DELLY: 
		  version: \$( echo \$(delly --version 2>&1) | sed 's/^.*Delly version: v//; s/ using.*\$//')
		  container: ${task.container}
		END_VERSIONS
		"""

	stub:
		"""
		touch ${id}.delly.vcf.gz

		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 DELLY: 
		  version: \$( echo \$(delly --version 2>&1) | sed 's/^.*Delly version: v//; s/ using.*\$//')
		  container: ${task.container}
		END_VERSIONS
		"""
}

process cnvkit_panel {
	cpus = 5
	container = '/fs1/resources/containers/twistmyeloid_active.sif'
	publishDir "/${OUTDIR}/sv_vcf/", mode: 'copy', overwrite: 'true', pattern: '*.vcf'
	publishDir "/${OUTDIR}/plots/", mode: 'copy', overwrite: 'true', pattern: '*.png'
	publishDir "/${OUTDIR}/versions", mode: 'copy' , overwrite: 'true', pattern: '*versions.yml'
	tag "$id"
	time '20m'
	memory '20 GB'
	scratch true
	stageInMode 'copy'
	stageOutMode 'copy'

	when:
		params.sv && params.onco

	input:
		set group, id, file(bam), file(bai), file(vcf), file(multi), val(INS_SIZE), val(MEAN_DEPTH), val(COV_DEV) from bam_cnvkit_panel.mix(bam_cnvkitpanel_choice).join(vcf_cnvkit, by:[0,1]).join(qc_cnvkit_val, by:[0,1]).view()
		//set id, val(INS_SIZE), val(MEAN_DEPTH), val(COV_DEV) from qc_cnvkit_val.view()
		//set group, id, file(vcf) from vcf_cnvkit.view()
	
	output:
		set group, id, file("${id}.cnvkit_filtered.vcf") into called_cnvkit_panel
		file("${id}.call.cns") into unfiltered_cns
		file("${group}.genomic_overview.png")
		set group, file("${group}_oplot.INFO") into cnvkit_INFO
		path "*versions.yml"

	script:
		"""
		cnvkit.py batch $bam -r $params.cnvkit_reference -p 5 -d results/
		cnvkit.py call results/*.cns -v $vcf -o ${id}.call.cns
		filter_cnvkit.pl ${id}.call.cns $MEAN_DEPTH > ${id}.filtered
		cnvkit.py export vcf ${id}.filtered -i "$id" > ${id}.cnvkit_filtered.vcf
		cnvkit.py scatter -s results/*dedup.cn{s,r} -o ${group}.genomic_overview.png -v $vcf -i $id
		echo "IMG overviewplot	${params.accessdir}/plots/${group}.genomic_overview.png" > ${group}_oplot.INFO
		
		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 CNVkit: 
		  version: \$(cnvkit.py version | sed -e 's/cnvkit v//g')
		  container: ${task.container}
		END_VERSIONS
		"""

	stub:
		"""
		touch ${id}.cnvkit_filtered.vcf
		touch ${id}.call.cns
		touch ${group}.genomic_overview.png
		touch ${group}_oplot.INFO

		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 CNVkit: 
		  version: \$(cnvkit.py version | sed -e 's/cnvkit v//g')
		  container: ${task.container}
		END_VERSIONS
		"""
}

process svdb_merge_panel {
	cpus 1
	cache 'deep'
	tag "$group"
	publishDir "/${OUTDIR}/sv_vcf/merged/", mode: 'copy', overwrite: 'true', pattern: '*.vcf'
	publishDir "/${OUTDIR}/versions", mode: 'copy' , overwrite: 'true', pattern: '*versions.yml'
	time '10m'
	memory '1 GB'
	scratch true
	stageInMode 'copy'
	stageOutMode 'copy'

	input:
		//set group, id, file(mantaV), file(dellyV), file(melt), file(cnvkitV) \
		//	from called_manta_panel.join(called_delly_panel, by:[0,1]).join(melt_vcf, by:[0,1]).join(called_cnvkit_panel, by:[0,1])
		set group, id, file(vcfs), id, file(melt) from called_manta_panel.mix(called_delly_panel,called_cnvkit_panel,merged_gatk_panel).groupTuple().join(melt_vcf).view()
				
	output:
		set group, id, file("${group}.merged.filtered.melt.vcf") into vep_sv_panel, annotsv_panel 
		//set group, id, file("${group}.merged.filtered.vcf") into annotsv_panel
		set group, file("${group}.merged.filtered.melt.vcf") into loqusdb_sv_panel
		path "*versions.yml"

	script:
		//tmp = mantaV.collect {it + ':manta ' } + dellyV.collect {it + ':delly ' } + cnvkitV.collect {it + ':cnvkit ' }
		//vcfs = tmp.join(' ')
		if (vcfs.size() > 1) {
			// for each sv-caller add idx, find vcf and find priority, add in priority order! //
			// index of vcfs added from mix //
			manta_idx = vcfs.findIndexOf{ it =~ 'manta' }
			delly_idx = vcfs.findIndexOf{ it =~ 'delly' }
			cnvkit_idx = vcfs.findIndexOf{ it =~ 'cnvkit' }
			gatk_idx = vcfs.findIndexOf{ it =~ 'gatk' }

			// find vcfs //
			manta = manta_idx >= 0 ? vcfs[manta_idx].collect {it + ':manta ' } : null
			delly = delly_idx >= 0 ? vcfs[delly_idx].collect {it + ':delly ' } : null
			cnvkit = cnvkit_idx >= 0 ? vcfs[cnvkit_idx].collect {it + ':cnvkit ' } : null
			gatk = gatk_idx >= 0 ? vcfs[gatk_idx].collect {it + ':gatk ' } : null
			tmp = manta + delly + gatk + cnvkit
			tmp = tmp - null
			vcfs_svdb = tmp.join(' ')

			// find priorities //
			mantap = manta_idx >= 0 ? 'manta' : null
			dellyp = delly_idx >= 0 ? 'delly' : null
			gatkp = gatk_idx >= 0 ? 'gatk' : null
			cnvkitp = cnvkit_idx >= 0 ? 'cnvkit' : null
			tmpp = [mantap, dellyp, gatkp, cnvkitp]
			tmpp = tmpp - null
			priority = tmpp.join(',')
			
			"""
			source activate py3-env
			svdb --merge --vcf $vcfs_svdb --no_intra --pass_only --bnd_distance 2500 --overlap 0.7 --priority $priority > ${group}.merged.vcf
			filter_panel_cnv.pl --mergedvcf ${group}.merged.vcf --callers $priority > ${group}.merged.filtered.vcf
			vcf-concat ${group}.merged.filtered.vcf $melt | vcf-sort -c > ${group}.merged.filtered.melt.vcf
			
			cat <<-END_VERSIONS > ${task.process}_versions.yml
			${task.process}:
			 SVDB: 
			  version: \$( echo \$(svdb) | head -1 | sed 's/usage: SVDB-\\([0-9]\\.[0-9]\\.[0-9]\\).*/\\1/' )
			  container: ${task.container}
			 SAMtools: 
			  version: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
			  container: ${task.container}
			END_VERSIONS
			"""
		}
		else {
			"""
			mv $vcf ${group}.merged.filtered.melt.vcf

			cat <<-END_VERSIONS > ${task.process}_versions.yml
			${task.process}:
			 SVDB: 
			  version: \$( echo \$(svdb) | head -1 | sed 's/usage: SVDB-\\([0-9]\\.[0-9]\\.[0-9]\\).*/\\1/' )
			  container: ${task.container}
			 SAMtools: 
			  version: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
			  container: ${task.container}
			END_VERSIONS
			"""
		}
//			vcf-concat ${group}.merged.filtered.vcf $melt | vcf-sort -c > ${group}.merged.filtered.melt.vcf // need to add melt somewhere!
	
	stub:
		"""
		touch ${group}.merged.filtered.melt.vcf

		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 SVDB: 
		  version: \$( echo \$(svdb) | head -1 | sed 's/usage: SVDB-\\([0-9]\\.[0-9]\\.[0-9]\\).*/\\1/' )
		  container: ${task.container}
		 SAMtools: 
		  version: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
		  container: ${task.container}
		END_VERSIONS
		"""
}

process tiddit {
	cpus = 2
	publishDir "/${OUTDIR}/sv_vcf/", mode: 'copy', overwrite: 'true', pattern: '*.vcf'
	publishDir "/${OUTDIR}/versions", mode: 'copy' , overwrite: 'true', pattern: '*versions.yml'
	time '10h'
	tag "$id"
	memory '10 GB'
	scratch true
	stageInMode 'copy'
	stageOutMode 'copy'

	when:
		params.sv && !params.onco &&  !params.exome


	input:
		set group, id, file(bam), file(bai) from bam_tiddit.mix(bam_tiddit_choice)

	output:
		set group, id, file("${id}.tiddit.filtered.vcf") into called_tiddit
		path "*versions.yml"

	script:
		"""
		TIDDIT.py --sv -o ${id}.tiddit --bam $bam
		grep -E \"#|PASS\" ${id}.tiddit.vcf > ${id}.tiddit.filtered.vcf

		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 tiddit: 
		  version: \$(echo \$(TIDDIT.py 2>&1) | sed 's/^.*TIDDIT-//; s/ .*\$//')
		  container: ${task.container}
		END_VERSIONS
		"""

	stub:
		"""
		touch ${id}.tiddit.filtered.vcf

		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 tiddit: 
		  version: \$(echo \$(TIDDIT.py 2>&1) | sed 's/^.*TIDDIT-//; s/ .*\$//')
		  container: ${task.container}
		END_VERSIONS
		"""
}

process svdb_merge {
	cpus 1
	tag "$group"
	publishDir "/${OUTDIR}/sv_vcf/merged/", mode: 'copy', overwrite: 'true', pattern: '*.vcf'
	publishDir "/${OUTDIR}/versions", mode: 'copy' , overwrite: 'true', pattern: '*versions.yml'
	time '2h'
	memory '1 GB'
	scratch true
	stageInMode 'copy'
	stageOutMode 'copy'

	input:
		set group, id, file(mantaV) from called_manta.groupTuple()
		set group, id, file(tidditV) from called_tiddit.groupTuple()
		set group, id, file(gatkV) from merged_gatk.groupTuple()
		
	output:
		set group, id, file("${group}.merged.bndless.vcf") into vcf_vep, annotsv_vcf
		set group, file("${group}.merged.vcf") into loqusdb_sv
		path "*versions.yml"

	script:
		if (mode == "family") {
			vcfs = []
			manta = []
			tiddit = []
			gatk = []
			for (i = 1; i <= mantaV.size(); i++) {
				tmp = mantaV[i-1] + ':manta' + "${i}"
				tmp1 = tidditV[i-1] + ':tiddit' + "${i}"
				tmp2 = gatkV[i-1] + ':gatk' + "${i}"
				vcfs = vcfs + tmp + tmp1 + tmp2
				mt = 'manta' + "${i}"
				tt = 'tiddit' + "${i}"
				ct = 'gatk' + "${i}"
				manta = manta + mt
				tiddit = tiddit + tt
				gatk = gatk + ct
			}
			prio = manta + tiddit + gatk
			prio = prio.join(',')
			vcfs = vcfs.join(' ')
			"""
			source activate py3-env
			svdb --merge --vcf $vcfs --no_intra --pass_only --bnd_distance 2500 --overlap 0.7 --priority $prio > ${group}.merged_tmp.vcf
			merge_callsets.pl ${group}.merged_tmp.vcf > ${group}.merged.vcf
			grep -v BND ${group}.merged.vcf > ${group}.merged.bndless.vcf
			
			cat <<-END_VERSIONS > ${task.process}_versions.yml
			${task.process}:
			 SVDB: 
			  version: \$( echo \$(svdb) | head -1 | sed 's/usage: SVDB-\\([0-9]\\.[0-9]\\.[0-9]\\).*/\\1/' )
			  container: ${task.container}
			 SAMtools: 
			  version: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
			  container: ${task.container}
			END_VERSIONS
			"""
		}

		else {
			tmp = mantaV.collect {it + ':manta ' } + tidditV.collect {it + ':tiddit ' } + gatkV.collect {it + ':gatk ' }
			vcfs = tmp.join(' ')
			"""
			source activate py3-env
			svdb --merge --vcf $vcfs --no_intra --pass_only --bnd_distance 2500 --overlap 0.7 --priority manta,tiddit,gatk > ${group}.merged_tmp.vcf
			merge_callsets.pl ${group}.merged_tmp.vcf > ${group}.merged.vcf
			grep -v BND ${group}.merged.vcf > ${group}.merged.bndless.vcf
			
			cat <<-END_VERSIONS > ${task.process}_versions.yml
			${task.process}:
			 SVDB: 
			  version: \$( echo \$(svdb) | head -1 | sed 's/usage: SVDB-\\([0-9]\\.[0-9]\\.[0-9]\\).*/\\1/' )
			  container: ${task.container}
			 SAMtools: 
			  version: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
			  container: ${task.container}
			END_VERSIONS
			"""
		}

	stub:
		"""
		touch ${group}.merged.vcf
		touch ${group}.merged.bndless.vcf

		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 SVDB: 
		  version: \$( echo \$(svdb) | head -1 | sed 's/usage: SVDB-\\([0-9]\\.[0-9]\\.[0-9]\\).*/\\1/' )
		  container: ${task.container}
		 SAMtools: 
		  version: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
		  container: ${task.container}
		END_VERSIONS
		"""
}

process add_to_loqusdb {
	cpus 1
	publishDir "${CRONDIR}/loqus", mode: 'copy' , overwrite: 'true'
	tag "$group"
	memory '1 MB'
	time '5m'

	when:
		!params.noupload

	input:
		set group, type, file(vcf), file(tbi), type, file(ped), file(svvcf) from vcf_loqus.join(ped_loqus).join(loqusdb_sv.mix(loqusdb_sv_panel)).view()

	output:
		file("${group}*.loqus") into loqusdb_done

	script:
		if (params.assay == "wgs") {
			"""
			echo "-db $params.loqusdb load -f ${params.accessdir}/ped/${ped} --variant-file ${params.accessdir}/vcf/${vcf} --sv-variants ${params.accessdir}/sv_vcf/merged/${svvcf}" > ${group}.loqus
			"""
		}
		else {
			"""
			echo "-db $params.loqusdb load -f ${params.accessdir}/ped/${ped} --variant-file ${params.accessdir}/vcf/${vcf} --sv-variants ${params.accessdir}/sv_vcf/merged/${svvcf}" > ${group}.loqus
			"""
		}

	stub:
		"""
		touch ${group}.loqus
		"""
}

//create AnnotSV tsv file
process annotsv {
	container = '/fs1/resources/containers/annotsv.v2.3.sif'
	cpus 2
	tag "$group"
	publishDir "/${OUTDIR}/annotsv/", mode: 'copy', overwrite: 'true', pattern: '*.tsv'
	publishDir "/${OUTDIR}/versions", mode: 'copy' , overwrite: 'true', pattern: '*versions.yml'
	time '5h'
	memory '20 GB'

	input:
		set group, id, file(sv) from annotsv_vcf.mix(annotsv_panel)
			
	output:
		set group, file("${group}_annotsv.tsv") into annotsv, annotsv_ma, annotsv_fa
		path "*versions.yml"

	script:
		"""
		export ANNOTSV="/AnnotSV"
		/AnnotSV/bin/AnnotSV -SvinputFile $sv \\
			-typeOfAnnotation full \\
			-outputDir $group \\
			-genomeBuild GRCh38
		mv $group/*.annotated.tsv ${group}_annotsv.tsv

		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 AnnotSV: 
		  version: \$( echo \$(/AnnotSV/bin/AnnotSV --version) | sed -e "s/AnnotSV //g ; s/Copyright.*//" )
		  container: ${task.container}
		END_VERSIONS
		"""

	stub:
		"""
		export ANNOTSV="/AnnotSV"
		touch ${group}_annotsv.tsv

		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 AnnotSV: 
		  version: \$( echo \$(/AnnotSV/bin/AnnotSV --version) | sed -e "s/AnnotSV //g ; s/Copyright.*//" )
		  container: ${task.container}
		END_VERSIONS
		"""
}

process vep_sv {
	publishDir "/${OUTDIR}/versions", mode: 'copy' , overwrite: 'true', pattern: '*versions.yml'
	cpus 56
	container = '/fs1/resources/containers/ensembl-vep_release_103.sif'
	tag "$group"
	memory '150 GB'
	time '1h'
	
	input:
		set group, id, file(vcf) from vcf_vep.mix(vep_sv_panel)

	output:
		set group, id, file("${group}.vep.vcf") into vep_vcf
		path "*versions.yml"

	script:
		"""
		vep \\
			-i $vcf \\
			-o ${group}.vep.vcf \\
			--offline \\
			--merged \\
			--everything \\
			--synonyms $params.SYNONYMS \\
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

		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 VEP: 
		  version: \$( echo \$(vep --help 2>&1) | sed 's/^.*Versions:.*ensembl-vep : // ; s/ .*\$//')
		  container: ${task.container}
		END_VERSIONS
		"""

	stub:
		"""
		touch ${group}.vep.vcf

		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 VEP: 
		  version: \$( echo \$(vep --help 2>&1) | sed 's/^.*Versions:.*ensembl-vep : // ; s/ .*\$//')
		  container: ${task.container}
		END_VERSIONS
		"""
}

process postprocess_vep {
	publishDir "/${OUTDIR}/versions", mode: 'copy' , overwrite: 'true', pattern: '*versions.yml'
	cpus = 1
	tag "$group"

	input:
		set group, id, file(vcf) from vep_vcf

	output:
		set group, file("${group}.vep.clean.merge.omim.vcf") into artefact_vcf
		path "*versions.yml"
	
	script:
		"""
		cleanVCF.py --vcf $vcf > ${group}.vep.clean.vcf
		svdb --merge --overlap 0.9 --notag --vcf ${group}.vep.clean.vcf > ${group}.vep.clean.merge.vcf
		sed -i '3 i ##INFO=<ID=set,Number=1,Type=String,Description="Source VCF for the merged record in SVDB">' ${group}.vep.clean.merge.vcf
		sed -i '3 i ##INFO=<ID=VARID,Number=1,Type=String,Description="The variant ID of merged samples">' ${group}.vep.clean.merge.vcf
		add_omim.pl ${group}.vep.clean.merge.vcf > ${group}.vep.clean.merge.omim.vcf

		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 SVDB: 
		  version: \$( echo \$(svdb) | head -1 | sed 's/usage: SVDB-\\([0-9]\\.[0-9]\\.[0-9]\\).*/\\1/' )
		  container: ${task.container}
		 SAMtools: 
		  version: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
		  container: ${task.container}
		END_VERSIONS
		"""

	stub:
		"""
		touch ${group}.vep.clean.merge.omim.vcf

		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 SVDB: 
		  version: \$( echo \$(svdb) | head -1 | sed 's/usage: SVDB-\\([0-9]\\.[0-9]\\.[0-9]\\).*/\\1/' )
		  container: ${task.container}
		 SAMtools: 
		  version: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
		  container: ${task.container}
		END_VERSIONS
		"""
}

// Query artefact db
process artefact {
	publishDir "/${OUTDIR}/versions", mode: 'copy' , overwrite: 'true', pattern: '*versions.yml'
	cpus 1
	tag "$group"
	time '10h'
	memory '10 GB'
	scratch true
	stageInMode 'copy'
	stageOutMode 'copy'

	input:
		set group, file(sv) from artefact_vcf

	output:
		set group, file("${group}.artefact.vcf") into manip_vcf,manip_vcf_ma,manip_vcf_fa
		path "*versions.yml"

	script:
		// use loqusdb dump not svdb database //
		if (params.gatkcnv) {
			"""
			source activate py3-env
			svdb \\
			--query --bnd_distance 25000 --overlap 0.7 --in_occ Obs --out_occ ACOUNT --in_frq Frq --out_frq AFRQ  \\
			--db $params.svdb \\
			--query_vcf $sv > ${group}.artefact.vcf

			cat <<-END_VERSIONS > ${task.process}_versions.yml
			${task.process}:
			 SVDB: 
			  version: \$( echo \$(svdb) | head -1 | sed 's/usage: SVDB-\\([0-9]\\.[0-9]\\.[0-9]\\).*/\\1/' )
			  container: ${task.container}
			 SAMtools: 
			  version: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
			  container: ${task.container}
			END_VERSIONS
			"""
		}
		// for oncov1-0 still use svdb database remove in future//
		else {
			"""
			source activate py3-env
			svdb \\
			--sqdb $params.svdb --query \\
			--query_vcf $sv --out_occ ACOUNT --out_frq AFRQ > ${group}.artefact.vcf

			cat <<-END_VERSIONS > ${task.process}_versions.yml
			${task.process}:
			 SVDB: 
			  version: \$( echo \$(svdb) | head -1 | sed 's/usage: SVDB-\\([0-9]\\.[0-9]\\.[0-9]\\).*/\\1/' )
			  container: ${task.container}
			 SAMtools: 
			  version: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
			  container: ${task.container}
			END_VERSIONS
			"""
		}

	stub:
		"""
		touch ${group}.artefact.vcf

		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 SVDB: 
		  version: \$( echo \$(svdb) | head -1 | sed 's/usage: SVDB-\\([0-9]\\.[0-9]\\.[0-9]\\).*/\\1/' )
		  container: ${task.container}
		 SAMtools: 
		  version: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
		  container: ${task.container}
		END_VERSIONS
		"""
}


process prescore {
	cpus 1
	tag "$group"
	memory '10 GB'
	time '30m'

	input:
		set group, file(sv_artefact), type, file(ped), file(annotsv) from manip_vcf.mix(manip_vcf_ma,manip_vcf_fa).join(ped_prescore.mix(ped_prescore_ma,ped_prescore_fa)).join(annotsv.mix(annotsv_ma,annotsv_fa))

	output:
		set group, type, file("${group}.annotatedSV.vcf") into annotatedSV

	script:
		"""
		prescore_sv.pl \\
		--sv $sv_artefact --ped $ped --annotsv $annotsv --osv ${group}.annotatedSV.vcf
		"""

	stub:
		"""
		touch ${group}.annotatedSV.vcf
		"""
}

process score_sv {
	cpus 5
	tag "$group $mode"
	publishDir "/${OUTDIR}/vcf", mode: 'copy', overwrite: 'true', pattern: '*.vcf.gz*'
	publishDir "/${OUTDIR}/versions", mode: 'copy' , overwrite: 'true', pattern: '*versions.yml'
	memory '10 GB'
	time '2h'
	container = '/fs1/resources/containers/genmod.sif'

	input:
		set group, type, file(vcf) from annotatedSV

	output:
		set group, type, file("${group_score}.sv.scored.sorted.vcf.gz"), file("${group_score}.sv.scored.sorted.vcf.gz.tbi") into sv_rescore,sv_rescore_ma,sv_rescore_fa
		set group, file("${group}_sv.INFO") into sv_INFO
		set group, file("${group_score}.sv.scored.sorted.vcf.gz") into svvcf_bed, svvcf_pod
		path "*versions.yml"

	script:
		group_score = group
		if ( type == "ma" || type == "fa") {
			group_score = group + "_" + type
		}
	
		if (mode == "family" && params.antype == "wgs") {
			"""
			genmod score -i $group_score -c $params.svrank_model -r $vcf -o ${group_score}.sv.scored_tmp.vcf
			bcftools sort -O v -o ${group_score}.sv.scored.sorted.vcf ${group_score}.sv.scored_tmp.vcf 
			bgzip -@ ${task.cpus} ${group_score}.sv.scored.sorted.vcf -f
			tabix ${group_score}.sv.scored.sorted.vcf.gz -f
			echo "SV	$type	${params.accessdir}/vcf/${group_score}.sv.scored.sorted.vcf.gz" > ${group}_sv.INFO

			cat <<-END_VERSIONS > ${task.process}_versions.yml
			${task.process}:
			 gunzip: 
			  version: \$(echo \$(gunzip --version 2>&1) | sed -e "s/^.*(gzip) // ; s/Copyright.*//")
			  container: ${task.container}
			 bgzip: 
			  version: \$(echo \$(bgzip --version 2>&1) | sed 's/^.*(htslib) // ; s/ Copyright.*//')
			  container: ${task.container}
			 tabix: 
			  version: \$(echo \$(tabix --version 2>&1) | sed 's/^.*(htslib) // ; s/ Copyright.*//')
			  container: ${task.container}
			 genmod: 
			  version: \$(echo \$(genmod --version 2>&1) | sed -e "s/^.*genmod version: //")
			  container: ${task.container}
			 bcftools: 
			  version: \$(echo \$(bcftools --version 2>&1) | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
			  container: ${task.container}
			END_VERSIONS
			"""
		}
		else {
			"""
			genmod score -i $group_score -c $params.svrank_model_s -r $vcf -o ${group_score}.sv.scored.vcf
			bcftools sort -O v -o ${group_score}.sv.scored.sorted.vcf ${group}.sv.scored.vcf
			bgzip -@ ${task.cpus} ${group_score}.sv.scored.sorted.vcf -f
			tabix ${group_score}.sv.scored.sorted.vcf.gz -f
			echo "SV	$type	${params.accessdir}/vcf/${group_score}.sv.scored.sorted.vcf.gz" > ${group}_sv.INFO

			cat <<-END_VERSIONS > ${task.process}_versions.yml
			${task.process}:
			 gunzip: 
			  version: \$(echo \$(gunzip --version 2>&1) | sed -e "s/^.*(gzip) // ; s/Copyright.*//")
			  container: ${task.container}
			 bgzip: 
			  version: \$(echo \$(bgzip --version 2>&1) | sed 's/^.*(htslib) // ; s/ Copyright.*//')
			  container: ${task.container}
			 tabix: 
			  version: \$(echo \$(tabix --version 2>&1) | sed 's/^.*(htslib) // ; s/ Copyright.*//')
			  container: ${task.container}
			 genmod: 
			  version: \$(echo \$(genmod --version 2>&1) | sed -e "s/^.*genmod version: //")
			  container: ${task.container}
			 bcftools: 
			  version: \$(echo \$(bcftools --version 2>&1) | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
			  container: ${task.container}
			END_VERSIONS
			"""
		}

	stub:
		group_score = group
		"""
		touch ${group_score}.sv.scored.sorted.vcf.gz
		touch ${group_score}.sv.scored.sorted.vcf.gz.tbi
		touch ${group}_sv.INFO

		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 gunzip: 
		  version: \$(echo \$(gunzip --version 2>&1) | sed -e "s/^.*(gzip) // ; s/Copyright.*//")
		  container: ${task.container}
		 bgzip: 
		  version: \$(echo \$(bgzip --version 2>&1) | sed 's/^.*(htslib) // ; s/ Copyright.*//')
		  container: ${task.container}
		 tabix: 
		  version: \$(echo \$(tabix --version 2>&1) | sed 's/^.*(htslib) // ; s/ Copyright.*//')
		  container: ${task.container}
		 genmod: 
		  version: \$(echo \$(genmod --version 2>&1) | sed -e "s/^.*genmod version: //")
		  container: ${task.container}
		 bcftools: 
		  version: \$(echo \$(bcftools --version 2>&1) | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
		  container: ${task.container}
		END_VERSIONS
		"""
}

process compound_finder {
	cpus 5
	tag "$group $mode"
	publishDir "/${OUTDIR}/vcf", mode: 'copy', overwrite: 'true', pattern: '*.vcf.gz*'
	publishDir "/${OUTDIR}/versions", mode: 'copy' , overwrite: 'true', pattern: '*versions.yml'
	memory '10 GB'
	time '2h'

	when:
		mode == "family" && params.assay == "wgs"

	input:
		set group, type, file(vcf), file(tbi), file(ped), file(snv), file(tbi) from sv_rescore.mix(sv_rescore_ma,sv_rescore_fa).join(ped_compound.mix(ped_compound_ma,ped_compound_fa), by: [0,1]).join(snv_sv_vcf.mix(snv_sv_vcf_ma,snv_sv_vcf_fa),by: [0,1])
		//set group, file(snv), file(tbi) from snv_sv_vcf

	output:
		set group, file("${group_score}.snv.rescored.sorted.vcf.gz"), file("${group_score}.snv.rescored.sorted.vcf.gz.tbi") into vcf_yaml
		set group, file("${group}_svp.INFO") into svcompound_INFO
		path "*versions.yml"

	script:
		group_score = group
		if ( type == "ma" || type == "fa") {
			group_score = group + "_" + type
		}

		"""
		compound_finder.pl \\
			--sv $vcf --ped $ped --snv $snv \\
			--osv ${group_score}.sv.rescored.sorted.vcf \\
			--osnv ${group_score}.snv.rescored.sorted.vcf \\
			--skipsv
		bgzip -@ ${task.cpus} ${group_score}.snv.rescored.sorted.vcf -f
		tabix ${group_score}.snv.rescored.sorted.vcf.gz -f
		echo "SVc	$type	${params.accessdir}/vcf/${group_score}.sv.scored.sorted.vcf.gz,${params.accessdir}/vcf/${group_score}.snv.rescored.sorted.vcf.gz" > ${group}_svp.INFO
		
		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 bgzip: 
		  version: \$(echo \$(bgzip --version 2>&1) | sed 's/^.*(htslib) // ; s/ Copyright.*//')
		  container: ${task.container}
		 tabix: 
		  version: \$(echo \$(tabix --version 2>&1) | sed 's/^.*(htslib) // ; s/ Copyright.*//')
		  container: ${task.container}
		END_VERSIONS
		"""

	stub:
		group_score = group
		"""
		touch ${group_score}.snv.rescored.sorted.vcf.gz
		touch ${group_score}.snv.rescored.sorted.vcf.gz.tbi
		touch ${group}_svp.INFO

		cat <<-END_VERSIONS > ${task.process}_versions.yml
		${task.process}:
		 bgzip: 
		  version: \$(echo \$(bgzip --version 2>&1) | sed 's/^.*(htslib) // ; s/ Copyright.*//')
		  container: ${task.container}
		 tabix: 
		  version: \$(echo \$(tabix --version 2>&1) | sed 's/^.*(htslib) // ; s/ Copyright.*//')
		  container: ${task.container}
		END_VERSIONS
		"""
}


process ouput_files {
	cpus 1
	memory '1MB'
	time '2m'

	input:
		set group, files from bam_INFO.mix(snv_INFO,sv_INFO,str_INFO,peddy_INFO,madde_INFO,svcompound_INFO,smn_INFO,bamchoice_INFO,mtBAM_INFO,oplot_INFO,haplogrep_INFO,eklipse_INFO,cnvkit_INFO).groupTuple()

	output:
		set group, file("${group}.INFO") into yaml_INFO

	script:
		files = files.join( ' ' )

		"""
		cat $files > ${group}.INFO
		"""

	stub:
		"""
		touch ${group}.INFO
		"""
}


process svvcf_to_bed {
	publishDir "/${OUTDIR}/bed", mode: 'copy' , overwrite: 'true'
	tag "group"
	memory '1 GB'
	time '10m'

	when:
		!params.onco && !params.exome

	input:
		set group, file(vcf) from svvcf_bed
		set group, id, sex, type from meta_svbed.filter { item -> item[3] == 'proband' }

	output:
		file("${group}.sv.bed")


	script:
		"""
		cnv2bed.pl --cnv $vcf --pb $id > ${group}.sv.bed
		"""

	stub:
		"""
		touch ${group}.sv.bed
		"""
}

process plot_pod {
	container = '/fs1/resources/containers/POD_2020-05-19.sif'
	publishDir "/${OUTDIR}/pod", mode: 'copy' , overwrite: 'true', pattern: '*.pdf'
	publishDir "/${OUTDIR}/pod", mode: 'copy' , overwrite: 'true', pattern: '*.html'
	tag "$group"
	time '20m'
	memory '1 GB'

	input:
		set group, file(snv) from vcf_pod
		set group, file(cnv), type, file(ped) from svvcf_pod.join(ped_pod)
		set group, id, sex, type from meta_pod.filter { item -> item[3] == 'proband' }		

	output:
		set file("${id}_POD_karyotype.pdf"), file("${id}_POD_results.html")

	when:
		mode == "family" && trio == true

	script:
		"""
		parental_origin_of_duplication.pl --snv $snv --cnv $cnv --proband $id --ped $ped
		"""

	stub:
		"""
		touch ${id}_POD_karyotype.pdf
		touch ${id}_POD_results.html
		"""
}

process create_yaml {
	publishDir "/${OUTDIR}/yaml", mode: 'copy' , overwrite: 'true', pattern: '*.yaml'
	publishDir "/${OUTDIR}/yaml/alt_affect", mode: 'copy' , overwrite: 'true', pattern: '*.yaml.*a'
	publishDir "${CRONDIR}/scout", mode: 'copy' , overwrite: 'true', pattern: '*.yaml'
	errorStrategy 'retry'
	maxErrors 5
	tag "$group"
	time '5m'
	memory '1 GB'

	input:
		set group, id, sex, mother, father, phenotype, diagnosis, type, assay, clarity_sample_id, ffpe, analysis, type, file(ped), file(INFO) from yml_diag.join(ped_scout).join(yaml_INFO)

	output:
		set group, file("${group}.yaml*") into yaml

	script:
		"""
		create_yml.pl \\
			--g $group,$clarity_sample_id --d $diagnosis --panelsdef $params.panelsdef --out ${group}.yaml --ped $ped --files $INFO --assay $assay,$analysis --antype $params.antype
		"""

	stub:
		"""
		touch ${group}.yaml
		"""
}
