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
	bam_gatk_choice }

vcf_choice.into{
	split_cadd_choice;
	split_vep_choice;
}
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
	memory '120 GB'
	tag "$id $shard"
	time '5h'

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
	time '1h'
	memory '120 GB'

	input:
		set val(id), group, file(shard), file(shard_bai) from bwa_shards_ch.groupTuple(by: [0,1])

	output:
		set id, group, file("${id}_merged.bam"), file("${id}_merged.bam.bai") into merged_bam_locusc
		set id, file("${id}_merged.bam"), file("${id}_merged.bam.bai") into merged_bam_dedup, merged_bam_bqsr

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
	memory '80 GB'
	// 64 GB peak giab //
	scratch true
	stageInMode 'copy'
	stageOutMode 'copy'
	tag "$id"

	input:
		set val(group), val(id), file(r1), file(r2) from fastq.mix(fastq_trimmed)

	output:
		set id, group, file("${id}_merged.bam"), file("${id}_merged.bam.bai") into bam_locusc, bam_markdup, bam_bqsr
		// set id, file("${id}_merged.bam"), file("${id}_merged.bam.bai") into bam_dedup remnant of distri

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
	publishDir "${OUTDIR}/bam", mode: 'copy' , overwrite: 'true', pattern: '*_dedup.bam*'

	input:
		set id, group, file(bam), file(bai) from bam_markdup.mix(merged_bam_dedup)

	output:
		set group, id, file("${id}_dedup.bam"), file("${id}_dedup.bam.bai") into complete_bam, chanjo_bam, expansionhunter_bam, yaml_bam, cov_bam, bam_manta, bam_nator, bam_tiddit, bam_manta_panel, bam_delly_panel, bam_cnvkit_panel, bam_freebayes, bam_mito, smncnc_bam, bam_gatk 
		set id, group, file("${id}_dedup.bam"), file("${id}_dedup.bam.bai") into qc_bam, bam_melt
		set val(id), file("dedup_metrics.txt") into dedupmet_sentieonqc
		set group, file("${group}_bam.INFO") into bam_INFO

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
		-i $bam --shard 1:1-248956422 --shard 2:1-242193529 --shard 3:1-198295559 --shard 4:1-190214555 --shard 5:1-120339935 --shard 5:120339936-181538259 --shard 6:1-170805979 --shard 7:1-159345973 --shard 8:1-145138636 --shard 9:1-138394717 --shard 10:1-133797422 --shard 11:1-135086622 --shard 12:1-56232327 --shard 12:56232328-133275309 --shard 13:1-114364328 --shard 14:1-107043718 --shard 15:1-101991189 --shard 16:1-90338345 --shard 17:1-83257441 --shard 18:1-80373285 --shard 19:1-58617616 --shard 20:1-64444167 --shard 21:1-46709983 --shard 22:1-50818468 --shard X:1-124998478 --shard X:124998479-156040895 --shard Y:1-57227415 --shard M:1-16569 \\
		--algo Dedup --score_info ${id}.score \\
		--metrics dedup_metrics.txt \\
		--rmdup ${id}_dedup.bam
	echo "BAM	$id	/access/${params.subdir}/bam/${id}_dedup.bam" > ${group}_bam.INFO
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
	publishDir "${OUTDIR}/bqsr", mode: 'copy' , overwrite: 'true', pattern: '*table'

	input:
		set id, group, file(bam), file(bai) from bam_bqsr.mix(merged_bam_bqsr)

	output:
		set group, id, file("${id}.bqsr.table") into dnascope_bqsr

	"""
	sentieon driver -t ${task.cpus} \\
		-r $genome_file -i $bam \\
		--algo QualCal ${id}.bqsr.table \\
		-k $params.KNOWN
	"""	
}

//Collect various QC data: 
process sentieon_qc {
	cpus 52
	memory '30 GB'
	publishDir "${OUTDIR}/qc", mode: 'copy' , overwrite: 'true', pattern: '*.QC'
	tag "$id"
	cache 'deep'
	time '2h'
	scratch true
	stageInMode 'copy'
	stageOutMode 'copy'

	input:
		set id, group, file(bam), file(bai), file(dedup) from qc_bam.mix(bam_qc_choice).join(dedupmet_sentieonqc.mix(dedup_dummy))

	output:
		set group, id, file("${id}.QC") into qc_cdm, qc_melt
		file("*.txt")

	script:
		target = ""
		panel = ""
		cov = "WgsMetricsAlgo wgs_metrics.txt"
		assay = "wgs"
		if( params.onco || params.exome) {
			target = "--interval $params.intervals"
			panel = params.panelhs + "${bam.toRealPath()}" + params.panelhs2 
			cov = "CoverageMetrics --cov_thresh 1 --cov_thresh 10 --cov_thresh 30 --cov_thresh 100 --cov_thresh 250 --cov_thresh 500 cov_metrics.txt"
			assay = "panel"
		}
		
	"""
	sentieon driver \\
		-r $genome_file $target \\
		-t ${task.cpus} \\
		-i ${bam.toRealPath()} \\
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
	time '10m'

	when:
		!params.noupload
	
	input:
		set group, id, file(qc), diagnosis, r1, r2 from qc_cdm.join(qc_extra)

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
	memory '10 GB'
	publishDir "${OUTDIR}/cov", mode: 'copy', overwrite: 'true'
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

	"""
	sambamba depth region -t ${task.cpus} -L $params.scoutbed -T 10 -T 15 -T 20 -T 50 -T 100 ${bam.toRealPath()} > ${id}.bwa.chanjo.cov
	"""
}

process SMNCopyNumberCaller {
	cpus 10
	memory '25GB'
	time '2h'
	publishDir "${OUTDIR}/plots/SMNcnc", mode: 'copy' , overwrite: 'true', pattern: '*.pdf*'
	tag "$id"

	when:
		params.antype == "wgs"

	input:
        set group, id, file(bam), file(bai) from smncnc_bam.mix(bam_SMN_choice)

	output:
		file("*.tsv") into smn_tsv
		set file("*.pdf"), file("*.json")
		set group, file("${group}_smn.INFO") into smn_INFO

	"""
	samtools view -H ${bam.toRealPath()} | \\
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
		samtools reheader - ${bam.toRealPath()} > ${id}.bam
	samtools index -b ${id}.bam -@ ${task.cpus}
	echo ${id}.bam > manifest.txt
	smn_caller.py --manifest manifest.txt --genome 38 --prefix ${id} --outDir . --threads ${task.cpus}
	rm ${id}.bam
	source activate py3-env
	python /SMNCopyNumberCaller/smn_charts.py -s ${id}.json -o .
	mv ${id}.tsv ${group}_SMN.tsv
	echo "SMN ${OUTDIR}/smn/${group}_SMN.tsv" > ${group}_smn.INFO
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
	tag "$group"
	cpus 2
	time '10h'
	memory '40 GB'
	scratch true
	stageInMode 'copy'
	stageOutMode 'copy'
	//publishDir "${OUTDIR}/plots/GAV/${group}", mode: 'copy' , overwrite: 'true', pattern: '*.png'

	when:
		params.str
		
	input:
		set group, id, file(bam), file(bai), sex, type \
			from expansionhunter_bam.mix(expansionhunter_bam_choice).join(meta_exp, by: [0,1]).filter { item -> item[5] == 'proband' }

	output:
		set group, id, file("${group}.eh.vcf") into expansionhunter_vcf
		set group, id, file("${group}.eh_realigned.sort.bam"), file("${group}.eh_realigned.sort.bam.bai"), file("${group}.eh.vcf") into reviewer

	"""
	source activate htslib10
	ExpansionHunter \
		--reads ${bam.toRealPath()} \
		--reference $genome_file \
		--variant-catalog $params.expansionhunter_catalog \
		--output-prefix ${group}.eh
	samtools sort ${group}.eh_realigned.bam -o ${group}.eh_realigned.sort.bam
	samtools index ${group}.eh_realigned.sort.bam
	"""
}

// annotate expansionhunter vcf
process stranger {
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

	"""
	stranger ${eh_vcf} -f $params.expansionhunter_catalog > ${group}.eh.stranger.vcf
	grep ^# ${group}.eh.stranger.vcf > ${group}.fixinfo.eh.stranger.vcf
    grep -v ^# ${group}.eh.stranger.vcf | sed 's/ /_/g' >> ${group}.fixinfo.eh.stranger.vcf
	"""
}

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
	publishDir "${OUTDIR}/plots/reviewer/${group}", mode: 'copy' , overwrite: 'true', pattern: '*.svg'
	
	input:
		set group, id, file(bam), file(bai), file(vcf) from reviewer

	output:
		file("*svg")

	shell:
	'''
	for locus in $(grep LocusId !{params.expansionhunter_catalog} | cut -f 2 -d ":" | sed 's/\"//g' | sed 's/,//g'); do 
		REViewer \
		--reads !{bam} \
		--vcf !{vcf} \
		--reference !{genome_file} \
		--catalog !{params.expansionhunter_catalog} \
		--locus $locus \
		--output-prefix 7156-15.${locus}
	done
    '''
}

// split multiallelic sites in expansionhunter vcf
// FIXME: Use env variable for picard path...
process vcfbreakmulti_expansionhunter {
	publishDir "${OUTDIR}/vcf", mode: 'copy' , overwrite: 'true'
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
			echo "STR	${OUTDIR}/vcf/${group}.expansionhunter.vcf.gz" > ${group}_str.INFO
			"""
		}
		else {
			"""
			java -jar /opt/conda/envs/CMD-WGS/share/picard-2.21.2-1/picard.jar RenameSampleInVcf INPUT=${eh_vcf_anno} OUTPUT=${eh_vcf_anno}.rename.vcf NEW_SAMPLE_NAME=${id}
			vcfbreakmulti ${eh_vcf_anno}.rename.vcf > ${group}.expansionhunter.vcf
			bgzip ${group}.expansionhunter.vcf
			tabix ${group}.expansionhunter.vcf.gz
			echo "STR	${OUTDIR}/vcf/${group}.expansionhunter.vcf.gz" > ${group}_str.INFO
			"""
		}
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
	memory '40 GB'
	time '3h'
	scratch true
	stageInMode 'copy'
	stageOutMode 'copy'

	input:
		set id, group, file(bam), file(bai), val(INS_SIZE), val(MEAN_DEPTH), val(COV_DEV) from bam_melt.mix(bam_melt_choice).join(qc_melt_val)

	when:
		params.onco

	output:
		set group, id, file("${id}.melt.merged.vcf") into melt_vcf

	"""
	java -jar  /opt/MELT.jar Single \\
		-bamfile ${bam.toRealPath()} \\
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
process dnascope {
	cpus 54
	memory '40 GB'
	// 12 GB peak giab //
	time '4h'
	tag "$id"

	when:
		params.varcall

	input:
		set group, id, bam, bai, bqsr from complete_bam.mix(dnascope_bam_choice).join(dnascope_bqsr, by: [0,1] )

	output:
		set group, id, file("${id}.dnascope.gvcf.gz"), file("${id}.dnascope.gvcf.gz.tbi") into complete_vcf_choice
		set group, id, file("${id}.dnascope.gvcf.gz") into gvcf_gens_choice
		//set group, file("${group}_bamstart.INFO") into bamchoice_INFO

	"""
	sentieon driver \\
		-t ${task.cpus} \\
		-r $genome_file \\
		-q $bqsr \\
		-i ${bam.toRealPath()} --shard 1:1-248956422 --shard 2:1-242193529 --shard 3:1-198295559 --shard 4:1-190214555 --shard 5:1-120339935 --shard 5:120339936-181538259 --shard 6:1-170805979 --shard 7:1-159345973 --shard 8:1-145138636 --shard 9:1-138394717 --shard 10:1-133797422 --shard 11:1-135086622 --shard 12:1-56232327 --shard 12:56232328-133275309 --shard 13:1-114364328 --shard 14:1-107043718 --shard 15:1-101991189 --shard 16:1-90338345 --shard 17:1-83257441 --shard 18:1-80373285 --shard 19:1-58617616 --shard 20:1-64444167 --shard 21:1-46709983 --shard 22:1-50818468 --shard X:1-124998478 --shard X:124998479-156040895 --shard Y:1-57227415 --shard M:1-16569 \\
		--algo DNAscope --emit_mode GVCF ${id}.dnascope.gvcf.gz
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

	"""
	echo "BAM	$id	/access/${params.subdir}/bam/${bam.getName()}" > ${group}_bamstart.INFO
	"""
}


process gvcf_combine {
	cpus 16
	tag "$group"
	memory '5 GB'
	time '5h'

	input:
		set group, id, file(vcf), file(idx) from complete_vcf_choice.mix(gvcf_choice).groupTuple()

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
	time '5m'

	input:
		set group, id, sex, mother, father, phenotype, diagnosis, type, assay, clarity_sample_id, ffpe, analysis from ped
		

	output:
		set group, id, file("${id}_slice.ped") into ped_ch
		set id, val(group) into madde_group
		//file("${group}.INFO") into tissue_INFO

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
	echo "${group}\t${id}\t${father}\t${mother}\t${sex}\t${phenotype}" > ${id}_slice.ped
	"""
}

// collects each individual's ped-line and creates one ped-file
// ped_ch
// 	.collectFile(sort: true, storeDir: "${OUTDIR}/ped/")
// 	.into{ ped_mad; ped_peddy; ped_inher; ped_scout; ped_loqus; ped_prescore; ped_compound; ped_pod }
process join_pedigree {
	cpus 1
	memory '1MB'
	time '1m'
	publishDir "${OUTDIR}/ped", mode: 'copy' , overwrite: 'true'

	input:
		set group, id, pedslice from ped_ch.groupTuple()

	output:
		set group, file("${group}.ped") into ped_mad, ped_peddy, ped_inher, ped_scout, ped_loqus, ped_prescore, ped_compound, ped_pod

	script:
	ped = pedslice.join( ' ' )

	"""
	cat $ped > ${group}.ped
	"""
}

//madeline ped, run if family mode
process madeline {
	publishDir "${OUTDIR}/ped", mode: 'copy' , overwrite: 'true'
	memory '1 GB'
	time '30m'
	container '/fs1/resources/containers/madeline.sif'

	input:
		set group, file(ped) from ped_mad
		set id, val(group) from madde_group

	output:
		file("${ped}.madeline.xml") into madeline_ped
		set group, file("${group}_madde.INFO") into madde_INFO

	when:
		mode == "family" && params.assay == "wgs"

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
	echo "MADDE ${OUTDIR}/ped/${ped}.madeline.xml" > ${group}_madde.INFO
	"""
}

process freebayes {
    cpus 1
    time '2h'
    container '/fs1/resources/containers/twistmyeloid_2020-06-17.sif'
	scratch true
	stageInMode 'copy'
	stageOutMode 'copy'
	cache 'deep'

	when: 
		params.onco

    input:
        set group, id, file(bam), file(bai) from bam_freebayes.mix(bam_freebayes_choice)

    output:
        set group, file("${id}.pathfreebayes.lines") into freebayes_concat

	script:
		if (params.onco) {
			"""
			freebayes -f $genome_file --pooled-continuous --pooled-discrete -t $params.intersect_bed --min-repeat-entropy 1 -F 0.03 ${bam.toRealPath()} > ${id}.freebayes.vcf
			vcfbreakmulti ${id}.freebayes.vcf > ${id}.freebayes.multibreak.vcf
			bcftools norm -m-both -c w -O v -f $genome_file -o ${id}.freebayes.multibreak.norm.vcf ${id}.freebayes.multibreak.vcf
			vcfanno_linux64 -lua /fs1/resources/ref/hg19/bed/scout/sv_tracks/silly.lua $params.vcfanno ${id}.freebayes.multibreak.norm.vcf > ${id}.freebayes.multibreak.norm.anno.vcf
			grep ^# ${id}.freebayes.multibreak.norm.anno.vcf > ${id}.freebayes.multibreak.norm.anno.path.vcf
			grep -v ^# ${id}.freebayes.multibreak.norm.anno.vcf | grep -i pathogenic > ${id}.freebayes.multibreak.norm.anno.path.vcf2
			cat ${id}.freebayes.multibreak.norm.anno.path.vcf ${id}.freebayes.multibreak.norm.anno.path.vcf2 > ${id}.freebayes.multibreak.norm.anno.path.vcf3
			filter_freebayes.pl ${id}.freebayes.multibreak.norm.anno.path.vcf3 > ${id}.pathfreebayes.lines
			"""
		}
		else {
			"""
			touch ${id}.pathfreebayes.lines
			"""
		}

}

/////////////// MITOCHONDRIA SNV CALLING ///////////////
///////////////                          ///////////////

// create an MT BAM file
process fetch_MTseqs {
	cpus 2
	memory '10GB'
	time '30m'
	tag "$id"
	publishDir "${OUTDIR}/bam", mode: 'copy', overwrite: 'true', pattern: '*.bam*'

	when:
		params.antype == "wgs"

    input:
        set group, id, file(bam), file(bai) from bam_mito.mix(bam_mito_choice)

    output:
        set group, id, file ("${id}_mito.bam"), file("${id}_mito.bam.bai") into mutserve_bam, eklipse_bam
		set group, file("${group}_mtbam.INFO") into mtBAM_INFO

    """
    sambamba view -f bam ${bam.toRealPath()} M > ${id}_mito.bam
    samtools index -b ${id}_mito.bam
	echo "mtBAM	$id	/access/${params.subdir}/bam/${id}_mito.bam" > ${group}_mtbam.INFO
    """

}

// gatk FilterMutectCalls in future if FPs overwhelms tord/sofie/carro
process run_mutect2 {
    cpus 4
    memory '50 GB'
    time '30m'
	tag "$group"
	publishDir "${OUTDIR}/vcf", mode: 'copy', overwrite: 'true'

	when:
		!params.onco
    
    input:
        set group, id, file(bam), file(bai) from mutserve_bam.groupTuple()

    output:
        set group, id, file("${group}.mutect2.vcf") into ms_vcfs_1, ms_vcfs_2

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
    """

}

// split and left-align variants
process split_normalize_mito {
    cpus 1
    memory '1GB'
    time '5m'

    input:
        set group, id, file(ms_vcf) from ms_vcfs_1
		set g2, id2, sex, type from meta_mutect2.groupTuple()

    output:
        set group, file("${group}.mutect2.breakmulti.filtered5p.0genotyped.proband.vcf") into adj_vcfs

	script:
		proband_idx = type.findIndexOf{ it == "proband" }

    """
    vcfbreakmulti $ms_vcf > ${ms_vcf}.breakmulti
    bcftools sort ${ms_vcf}.breakmulti | bgzip > ${ms_vcf}.breakmulti.fix
    tabix -p vcf ${ms_vcf}.breakmulti.fix
    bcftools norm -f $params.rCRS_fasta -o ${ms_vcf.baseName}.adjusted.vcf ${ms_vcf}.breakmulti.fix
	bcftools view -i 'FMT/AF[*]>0.05' ${ms_vcf.baseName}.adjusted.vcf -o ${group}.mutect2.breakmulti.filtered5p.vcf
	bcftools filter -S 0 --exclude 'FMT/AF[*]<0.05' ${group}.mutect2.breakmulti.filtered5p.vcf -o ${group}.mutect2.breakmulti.filtered5p.0genotyped.vcf
	filter_mutect2_mito.pl ${group}.mutect2.breakmulti.filtered5p.0genotyped.vcf ${id2[proband_idx]} > ${group}.mutect2.breakmulti.filtered5p.0genotyped.proband.vcf
    """

}

// use python tool HmtNote for annotating vcf
// future merging with diploid genome does not approve spaces in info-string
process run_hmtnote {
    cpus 1
    memory '5GB'
    time '15m'


    input:
        set group, file(adj_vcf) from adj_vcfs

    output:
        set group, file("${group}.fixinfo.vcf") into mito_diplod_vep

    """
    source activate tools
    hmtnote annotate ${adj_vcf} ${group}.hmtnote --offline
    grep ^# ${group}.hmtnote > ${group}.fixinfo.vcf
    grep -v ^# ${group}.hmtnote | sed 's/ /_/g' >> ${group}.fixinfo.vcf
    """
    
}

// run haplogrep 2 on resulting vcf
process run_haplogrep {
    time '10m'
    memory '30 GB'
    cpus '2'
	publishDir "${OUTDIR}/plots/mito", mode: 'copy', overwrite: 'true'

    input:
        set group, id, file(ms_vcf) from ms_vcfs_2

    output:
       file("${group}.haplogrep.png")

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
    '''

}

// use eKLIPse for detecting mitochondrial deletions
process run_eklipse {
    cpus 2
    memory '10GB'
    time '60m'
	publishDir "${OUTDIR}/plots/mito", mode: 'copy', overwrite: 'true'

    input:
        set group, id, file(bam), file(bai) from eklipse_bam

	output:
		set file("*.png"), file("${id}.hetplasmid_frequency.txt")
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
    """

}



// Splitting & normalizing variants, merging with Freebayes/Mutect2, intersecting against exome/clinvar introns
process split_normalize {
	cpus 1
	publishDir "${OUTDIR}/vcf", mode: 'copy', overwrite: 'true'
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
		set group, id, file("${group}.intersected.vcf"), file("${group}.multibreak.vcf") into split_vep, split_cadd, vcf_loqus, vcf_cnvkit
		
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
		"""

	}



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
		set group, id, file(vcf), file(multi), file(ped) from vcf_loqus.join(ped_loqus)
		//file(ped) from ped_loqus

	output:
		file("${group}.loqus") into loqusdb_done

	"""
	echo "-db $params.loqusdb load -f ${OUTDIR}/ped/${ped} --variant-file ${vcf.toRealPath()}" > ${group}.loqus
	"""
}

process annotate_vep {
	container = '/fs1/resources/containers/ensembl-vep_release_103.sif'
	cpus 54
	tag "$group"
	memory '150 GB'
	time '5h'
	scratch true
	stageInMode 'copy'
	stageOutMode 'copy'

	input:
		set group, id, file(vcf), idx from split_vep.mix(split_vep_choice)

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
		--synonyms $params.SYNONYMS \\
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
		-custom $params.GNOMAD_MT,gnomAD_mt,vcf,exact,0,AF_hom,AF_het \\
		-custom $params.PHYLOP \\
		-custom $params.PHASTCONS
	"""
}

// gene, clinvar, loqusdb, enigma(onco)
process vcfanno {
	cpus params.cpu_some
	memory '32GB'
	time '20m'
	scratch true
	stageInMode 'copy'
	stageOutMode 'copy'

	input:
		set group, file(vcf) from vep

	output:
		set group, file("${group}.clinvar.loqusdb.gene.vcf") into vcfanno_vcf

	"""
	vcfanno_linux64 -lua /fs1/resources/ref/hg19/bed/scout/sv_tracks/silly.lua $params.vcfanno $vcf > ${group}.clinvar.loqusdb.gene.vcf
	"""
}


// # Annotating variant inheritance models:
process inher_models {
	cpus 6
	memory '64 GB'
	tag "$group"
	time '10m'
	scratch true
	stageInMode 'copy'
	stageOutMode 'copy'

	input:
		set group, file(vcf), file(ped) from vcfanno_vcf.join(ped_inher)
		//file(ped) from ped_inher

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
	memory '1 GB'
	time '10m'

	input:
		set group, file(vcf) from inhermod

	output:
		set group, file("${group}.mod.vcf") into mod_vcf

	"""
	modify_vcf_scout.pl $vcf > ${group}.mod.vcf
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

	"""
	/opt/bin/mark_spliceindels.pl $vcf > ${group}.marksplice.vcf
	"""
}

// Extract all INDELs from VCF for CADD annotation
process extract_indels_for_cadd {
	cpus 1
	tag "$group"
	memory '1 GB'
	time '5m'

	input:
		set group, id, file(vcf), idx from split_cadd.mix(split_cadd_choice)
	
	output:
		set group, file("${group}.only_indels.vcf") into indel_cadd_vep

	"""
	bcftools view $vcf -V snps -o ${group}.only_indels.vcf 
	"""    
}

// Annotate Indels with VEP+Gnomad genomes. Filter variants below threshold
process indel_vep {
	cpus 5
	container = '/fs1/resources/containers/ensembl-vep_release_103.sif'
	tag "$group"
	memory '10 GB'
	time '3h'

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
		--synonyms $params.SYNONYMS \\
		--fasta $params.VEP_FASTA \\
		-custom $params.GNOMAD_GENOMES,gnomADg,vcf,exact,0,AF \\
		-custom $params.GNOMAD_MT,gnomAD_mt,vcf,exact,0,AF_hom,AF_het \\
		--dir_cache $params.VEP_CACHE \\
		--force_overwrite \\
		--no_stats \\
		--fork ${task.cpus}
	filter_indels.pl ${group}.only_indels.vep.vcf > ${group}.only_indels.vep.filtered.vcf
	"""
}

// Calculate CADD scores for all indels
process calculate_indel_cadd {
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

	"""
	/CADD-scripts/CADD.sh -c ${task.cpus} -g GRCh38 -o ${group}.indel_cadd.gz $vcf
	"""
}

// Add the calculated indel CADDs to the vcf
process add_cadd_scores_to_vcf {
	cpus 4
	tag "$group"
	memory '1 GB'
	time '5m'

	input: 
		set group, file(vcf), file(cadd_scores) from splice_marked.join(indel_cadd)

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
	memory '10 GB'
	time '30m'

	input:
		set group, file(vcf) from indel_cadd_added

	output:
		set group, file("${group}.scored.vcf") into scored_vcf

	script:
		if ( mode == "family" && params.antype == "wgs" ) {
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
		set group, file("${group}_snv.INFO") into snv_INFO

	"""
	sed 's/^MT/M/' -i $vcf
	sed 's/ID=MT,length/ID=M,length/' -i $vcf
	bgzip -@ ${task.cpus} $vcf -f
	tabix ${vcf}.gz -f
	echo "SNV	${OUTDIR}/vcf/${group}.scored.vcf.gz" > ${group}_snv.INFO
	"""
}


// Running PEDDY: 
process peddy {
	publishDir "${OUTDIR}/ped", mode: 'copy' , overwrite: 'true'
	//container = '/fs1/resources/containers/wgs_20200115.sif'
	cpus 6
	tag "$group"

	input:
		//file(ped) from ped_peddy
		set group, file(vcf), file(idx), file(ped) from vcf_peddy.join(ped_peddy)

	output:
		set file("${group}.ped_check.csv"),file("${group}.peddy.ped"), file("${group}.sex_check.csv") into peddy_files
		set group, file("${group}_peddy.INFO") into peddy_INFO

	"""
	source activate py3-env
	python -m peddy --sites hg38 -p ${task.cpus} $vcf $ped --prefix $group
	echo "PEDDY	${OUTDIR}/ped/${group}.ped_check.csv,${OUTDIR}/ped/${group}.peddy.ped,${OUTDIR}/ped/${group}.sex_check.csv" > ${group}_peddy.INFO
	"""
}

// Extract all variants (from whole genome) with a gnomAD af > x%
process fastgnomad {
	cpus 2
	memory '32 GB'
	tag "$group"
	publishDir "${OUTDIR}/vcf", mode: 'copy', overwrite: 'true'
	time '2h'

	when:
		!params.onco && !params.exome

	input:
		set group, file(vcf) from vcf_gnomad

	output:
		set group, file("${group}.SNPs.vcf") into vcf_upd, vcf_roh, vcf_pod

	"""
	gzip -c $vcf > ${vcf}.gz
	annotate -g $params.FASTGNOMAD_REF -i ${vcf}.gz > ${group}.SNPs.vcf
	"""
	
}


// Call UPD regions from SNP vcf
process upd {
	tag "$group"
	time '10m'
	memory '1 GB'

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
	time '10m'
	memory '1 GB'

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
	time '10m'
	memory '1 GB'

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
	memory '60 GB'
	time '5h'

	input:
		set id, group, file(bam), file(bai), gr, sex, type from cov_bam.mix(cov_bam_choice).join(meta_gatkcov, by:1)

	output:
		set group, id, type, sex, file("${id}.standardizedCR.tsv"), file("${id}.denoisedCR.tsv") into cov_plot, cov_gens

	when:
		params.gatkcov

	"""
	source activate gatk4-env

	gatk CollectReadCounts \\
		-I ${bam.toRealPath()} -L $params.COV_INTERVAL_LIST \\
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
	time '20m'
	memory '5 GB'

	input:
		file(upd) from upd_plot
		set gr, file(roh) from roh_plot
		set group, id, type, sex, file(cov_stand), file(cov_denoised) from cov_plot.groupTuple()


	output:
		file("${group}.genomic_overview.png")

	script:
		proband_idx = type.findIndexOf{ it == "proband" }

	"""
	genome_plotter.pl --dict $params.GENOMEDICT \\
		 --sample ${id[proband_idx]} \\
		 --upd $upd \\
		 --roh $roh \\
		 --sex ${sex[proband_idx]} \\
		 --cov ${cov_denoised[proband_idx]} \\
		 --out ${group}.genomic_overview.png
	"""
}

process generate_gens_data {
	publishDir "${OUTDIR}/plot_data", mode: 'copy' , overwrite: 'true'
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

	"""
	generate_gens_data.pl $cov_stand $gvcf $id $params.GENS_GNOMAD
	"""
}

process manta {
	cpus = 56
	publishDir "${OUTDIR}/sv_vcf/", mode: 'copy', overwrite: 'true'
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

	script:
		bams = bam.join('--bam ')

	"""
	configManta.py --bam ${bam.toRealPath()} --reference $genome_file --runDir .
	python runWorkflow.py -m local -j ${task.cpus}
	mv results/variants/diploidSV.vcf.gz ${id}.manta.vcf.gz
	mv results/variants/diploidSV.vcf.gz.tbi ${id}.manta.vcf.gz.tbi
	"""
}

process manta_panel {
	cpus = 56
	publishDir "${OUTDIR}/sv_vcf/", mode: 'copy', overwrite: 'true'
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


	"""
	configManta.py --bam ${bam.toRealPath()} --reference $genome_file --runDir . --exome --callRegions $params.bedgz --generateEvidenceBam
	python runWorkflow.py -m local -j ${task.cpus}
	mv results/variants/diploidSV.vcf.gz ${id}.manta.vcf.gz
	mv results/variants/diploidSV.vcf.gz.tbi ${id}.manta.vcf.gz.tbi
	"""
}

process delly_panel {
	cpus = 5
	publishDir "${OUTDIR}/sv_vcf/", mode: 'copy', overwrite: 'true'
	tag "$id"
	time '3h'
	memory '10 GB'
	scratch true
	stageInMode 'copy'
	stageOutMode 'copy'
	cache 'deep'
	
	when:
		params.sv && params.onco

	input:
		set group, id, file(bam), file(bai) from bam_delly_panel.mix(bam_dellypanel_choice)

	output:
		set group, id, file("${id}.delly.vcf.gz") into called_delly_panel


	"""
	delly call -g $genome_file -o ${id}.bcf ${bam.toRealPath()}
	bcftools view ${id}.bcf > ${id}.vcf
	filter_delly.pl --vcf ${id}.vcf --bed $params.intersect_bed > ${id}.delly.vcf
	bgzip -c ${id}.delly.vcf > ${id}.delly.vcf.gz
	"""
}

process cnvkit_panel {
	cpus = 5
	container = '/fs1/resources/containers/twistmyeloid_active.sif'
	publishDir "${OUTDIR}/sv_vcf/", mode: 'copy', overwrite: 'true'
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

	"""
	cnvkit.py batch ${bam.toRealPath()} -r $params.cnvkit_reference -p ${task.cpus} -d results/
	cnvkit.py call results/*.cns -v $vcf -o ${id}.call.cns
	filter_cnvkit.pl ${id}.call.cns $MEAN_DEPTH > ${id}.filtered
	cnvkit.py export vcf ${id}.filtered -i "$id" > ${id}.cnvkit_filtered.vcf
	"""

}

process svdb_merge_panel {
	cpus 1
	cache 'deep'
	tag "$group"
	publishDir "${OUTDIR}/sv_vcf/merged/", mode: 'copy', overwrite: 'true'
	time '10m'
	memory '1 GB'
	scratch true
	stageInMode 'copy'
	stageOutMode 'copy'

	input:
		set group, id, file(mantaV), file(dellyV), file(melt), file(cnvkitV) from called_manta_panel.join(called_delly_panel, by:[0,1]).join(melt_vcf, by:[0,1]).join(called_cnvkit_panel, by:[0,1])
		// set group, id, file(dellyV) from called_delly_panel.groupTuple()
		// set group, id, file(melt) from melt_vcf.groupTuple()
		// set group, id, file(cnvkitV) from called_cnvkit_panel.groupTuple()
				
	output:
		set group, id, file("${group}.merged.filtered.melt.vcf") into vep_sv_panel, annotsv_panel 
		//set group, id, file("${group}.merged.filtered.vcf") into annotsv_panel

	script:
		tmp = mantaV.collect {it + ':manta ' } + dellyV.collect {it + ':delly ' } + cnvkitV.collect {it + ':cnvkit ' }
		vcfs = tmp.join(' ')

		"""
		source activate py3-env
		svdb --merge --vcf $vcfs --no_intra --pass_only --bnd_distance 2500 --overlap 0.7 --priority manta,delly,cnvkit > ${group}.merged.vcf
		filter_panel_cnv.pl ${group}.merged.vcf > ${group}.merged.filtered.vcf
		vcf-concat ${group}.merged.filtered.vcf $melt | vcf-sort -c > ${group}.merged.filtered.melt.vcf
		"""


}

process tiddit {
	cpus = 2
	publishDir "${OUTDIR}/sv_vcf/", mode: 'copy', overwrite: 'true'   
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

	"""
	TIDDIT.py --sv -o ${id}.tiddit --bam ${bam.toRealPath()}
	grep -E \"#|PASS\" ${id}.tiddit.vcf > ${id}.tiddit.filtered.vcf
	"""
}

process gatk_coverage {
    cpus 10
    memory '40GB'
    time '2h'
    container = '/fs1/resources/containers/gatk_4.1.9.0.sif'
    scratch true
	stageInMode 'copy'
	stageOutMode 'copy'
    tag "$id"   

	when:
		params.sv && !params.onco  &&  !params.exome

    input:
        set group, id, file(bam), file(bai) from bam_gatk.mix(bam_gatk_choice)

    output:
        set group, id, file("${id}.tsv") into call_ploidy, call_cnv

    """
	THEANO_FLAGS="base_compiledir=/fs1/resources/theano"
    export MKL_NUM_THREADS=${task.cpus}
    export OMP_NUM_THREADS=${task.cpus}
	set +u
	source activate gatk
    gatk --java-options "-Xmx20g" CollectReadCounts \\
        -L $params.gatk_intervals \\
        -R $params.genome_file \\
        -imr OVERLAPPING_ONLY \\
        -I ${bam.toRealPath()} \\
        --format TSV -O ${id}.tsv
    """
}

process gatk_call_ploidy {
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

    """
	THEANO_FLAGS="base_compiledir=/fs1/resources/theano"
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
    """
}

process gatk_call_cnv {
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
    """
}

process postprocessgatk {
    cpus 8
    memory '40GB'
    time '3h'
    container = '/fs1/resources/containers/gatk_4.1.9.0.sif'
	
    //scratch true
	// stageInMode 'copy'
	// stageOutMode 'copy'
    publishDir "${OUTDIR}/sv_vcf/", mode: 'copy', overwrite: 'true'
    tag "$id"


    input:
        set group, id, i, file(tar), file(ploidy), shard_no, shard \
			from postprocessgatk.groupTuple(by: [0,1]).join(ploidy_to_post, by: [0,1]).combine(gatk_postprocess.groupTuple(by: [3]))


    output:
        set group, id, \
            file("genotyped-intervals-${group}-vs-cohort30.vcf.gz"), \
            file("genotyped-segments-${group}-vs-cohort30.vcf.gz"), \
            file("denoised-${group}-vs-cohort30.vcf.gz") into called_gatk

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
	'''

}

process filter_merge_gatk {
	cpus 1
	tag "$group"
	time '2h'
	memory '1 GB'
	publishDir "${OUTDIR}/sv_vcf", mode: 'copy', overwrite: 'true'

	input:
		set group, id, file(inter), file(gatk), file(denoised) from called_gatk

	output:
		set group, id, file("${id}.gatk.filtred.merged.vcf") into merged_gatk

	"""
	filter_gatk.pl $gatk > ${id}.gatk.filtered.vcf
	mergeGATK.pl ${id}.gatk.filtered.vcf > ${id}.gatk.filtred.merged.vcf
	"""
}

process svdb_merge {
	cpus 1
	tag "$group"
	publishDir "${OUTDIR}/sv_vcf/merged/", mode: 'copy', overwrite: 'true'
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
			"""
		}

}

//create AnnotSV tsv file
process annotsv {
	container = '/fs1/resources/containers/annotsv.v2.3.sif'
	cpus 2
	tag "$group"
	publishDir "${OUTDIR}/annotsv/", mode: 'copy', overwrite: 'true'
	time '5h'
	memory '20 GB'

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
	container = '/fs1/resources/containers/ensembl-vep_release_103.sif'
	tag "$group"
	memory '150 GB'
	time '1h'
	
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
	cleanVCF.py --vcf $vcf > ${group}.vep.clean.vcf
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
	time '10h'
	memory '10 GB'
	scratch true
	stageInMode 'copy'
	stageOutMode 'copy'

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
	memory '10 GB'
	time '30m'

	input:
		set group, file(sv_artefact), file(ped), file(annotsv) from manip_vcf.join(ped_prescore).join(annotsv)
		//set group,  from annotsv

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
	memory '10 GB'
	time '2h'

	input:
		set group, file(vcf) from annotatedSV

	output:
		set group, file("${group}.sv.scored.sorted.vcf.gz"), file("${group}.sv.scored.sorted.vcf.gz.tbi") into sv_rescore
		set group, file("${group}_sv.INFO") into sv_INFO
		set group, file("${group}.sv.scored.sorted.vcf.gz") into svvcf_bed, svvcf_pod
				
	script:
	
		if (mode == "family" && params.antype == "wgs") {
			"""
			genmod score -i $group -c $params.svrank_model -r $vcf -o ${group}.sv.scored_tmp.vcf
			bcftools sort -O v -o ${group}.sv.scored.sorted.vcf ${group}.sv.scored_tmp.vcf 
			bgzip -@ ${task.cpus} ${group}.sv.scored.sorted.vcf -f
			tabix ${group}.sv.scored.sorted.vcf.gz -f
			echo "SV	${OUTDIR}/vcf/${group}.sv.scored.sorted.vcf.gz" > ${group}_sv.INFO
			"""
		}
		else {
			"""
			genmod score -i $group -c $params.svrank_model_s -r $vcf -o ${group}.sv.scored.vcf
			bcftools sort -O v -o ${group}.sv.scored.sorted.vcf ${group}.sv.scored.vcf
			bgzip -@ ${task.cpus} ${group}.sv.scored.sorted.vcf -f
			tabix ${group}.sv.scored.sorted.vcf.gz -f
			echo "SV	${OUTDIR}/vcf/${group}.sv.scored.sorted.vcf.gz" > ${group}_sv.INFO
			"""
		}
}

process compound_finder {
	cpus 5
	tag "$group $mode"
	publishDir "${OUTDIR}/vcf", mode: 'copy', overwrite: 'true', pattern: '*.vcf.gz*'
	memory '10 GB'
	time '2h'

	when:
		mode == "family" && params.assay == "wgs"

	input:
		set group, file(vcf), file(tbi), file(ped) from sv_rescore.join(ped_compound)
		set group, file(snv), file(tbi) from snv_sv_vcf

	output:
		set group, file("${group}.snv.rescored.sorted.vcf.gz"), file("${group}.snv.rescored.sorted.vcf.gz.tbi") into vcf_yaml
			//file("${group}.sv.rescored.sorted.vcf.gz"), file("${group}.sv.rescored.sorted.vcf.gz.tbi") 
		set group, file("${group}_svp.INFO") into svcompound_INFO
				
	// OBS, if flag --skipsv is chosen modify what sv.vcf is presented to scout in yaml (last line)
	// and change output file.
	script:
		"""
		compound_finder.pl \\
			--sv $vcf --ped $ped --snv $snv \\
			--osv ${group}.sv.rescored.sorted.vcf \\
			--osnv ${group}.snv.rescored.sorted.vcf \\
			--skipsv
		bgzip -@ ${task.cpus} ${group}.snv.rescored.sorted.vcf -f
		tabix ${group}.snv.rescored.sorted.vcf.gz -f
		echo "SVc	${OUTDIR}/vcf/${group}.sv.scored.sorted.vcf.gz,${OUTDIR}/vcf/${group}.snv.rescored.sorted.vcf.gz" > ${group}_svp.INFO
		"""

}

// Collects $group.INFO files from each process output that should be included in the yaml for scout loading //
// If a new process needs to be added to yaml. It needs to follow this procedure, as well as be handled in create_yml.pl //
// bam_INFO
// 	.mix(snv_INFO,sv_INFO,str_INFO,peddy_INFO,madde_INFO,svcompound_INFO,smn_INFO,bamchoice_INFO,mtBAM_INFO)
// 	.collectFile()
// 	.set{ yaml_INFO }

process ouput_files {
	cpus 1
	memory '1MB'
	time '2m'

	input:
		set group, files from bam_INFO.mix(snv_INFO,sv_INFO,str_INFO,peddy_INFO,madde_INFO,svcompound_INFO,smn_INFO,bamchoice_INFO,mtBAM_INFO).groupTuple()

	output:
		set group, file("${group}.INFO") into yaml_INFO

	script:
		files = files.join( ' ' )

	"""
	cat $files > ${group}.INFO
	"""
}


process svvcf_to_bed {
	publishDir "${OUTDIR}/bed", mode: 'copy' , overwrite: 'true'
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


	"""
	cnv2bed.pl --cnv $vcf --pb $id > ${group}.sv.bed
	"""
}

process plot_pod {
	container = '/fs1/resources/containers/POD_2020-05-19.sif'
	publishDir "${OUTDIR}/pod", mode: 'copy' , overwrite: 'true'
	tag "$group"
	time '20m'
	memory '1 GB'

	input:
		set group, file(snv) from vcf_pod
		set group, file(cnv), file(ped) from svvcf_pod.join(ped_pod)
		set group, id, sex, type from meta_pod.filter { item -> item[3] == 'proband' }		

	output:
		set file("${id}_POD_karyotype.pdf"), file("${id}_POD_results.html")

	when:
		mode == "family" && trio == true

	"""
	parental_origin_of_duplication.pl --snv $snv --cnv $cnv --proband $id --ped $ped
	"""
}

process create_yaml {
	queue 'bigmem'
	publishDir "${OUTDIR}/yaml", mode: 'copy' , overwrite: 'true'
	publishDir "${CRONDIR}/scout", mode: 'copy' , overwrite: 'true'
	errorStrategy 'retry'
	maxErrors 5
	tag "$group"
	time '5m'
	memory '1 GB'


	when:
		!params.noupload

	input:
		//file(INFO) from yaml_INFO
		//file(ped) from ped_scout
		set group, id, sex, mother, father, phenotype, diagnosis, type, assay, clarity_sample_id, ffpe, analysis, file(ped), file(INFO) from yml_diag.join(ped_scout).join(yaml_INFO)

	output:
		set group, file("${group}.yaml") into yaml

	script:

	"""
	export PORT_CMDSCOUT2_MONGODB=33002 #TA BORT VÄLDIGT FULT
	create_yml.pl \\
		--g $group,$clarity_sample_id --d $diagnosis --panelsdef $params.panelsdef --out ${group}.yaml --ped $ped --files $INFO --assay $assay,$analysis --antype $params.antype
	"""
}
