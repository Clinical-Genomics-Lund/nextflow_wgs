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

// BED FILES //
scoutbed = params.scoutbed

PON = [F: params.GATK_PON_FEMALE, M: params.GATK_PON_MALE]

// Count lines of input csv, if more than 2(header + 1 ind) then mode is set to family //
csv = file(params.csv)
mode = csv.countLines() > 2 ? "family" : "single"
trio = csv.countLines() > 3 ? true : false
println(csv)
println(mode)
println(trio)

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

// If input files are fastq -> normal path. Flags affecting; --shardbwa (sharded bwa) --align(req), --varcall(if variant calling is to be done) and --annotate(if --varcall)
// bam -> skips align and is variant called (if --varcall is present) and annotated (if --annotate is present)
// vcf skips align + varcall and is only annotated (if --annotate is present)

// If .bam -> value 1, else if .vcf -> value 2 else if .gvcf -> value 4 else if none(.fq.gz) if params.shardbwa true -> value 3 otherwise 0
input_files.view().choice(fastq, bam_choice, vcf_choice, fastq_sharded, gvcf_choice) { it[2]  =~ /\.bam/ ? 1 : ( it[2] =~ /\.vcf/ ? 2 : ( it[2] =~ /\.gvcf/ ? 4 : (params.shardbwa ? 3 : 0)))  }

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
    .map{ row-> tuple(row.group, row.id, row.sex, row.mother, row.father, row.phenotype, row.diagnosis) }
    .into { ped; yml_diag; meta_upd }


Channel
    .fromPath(params.csv)
    .splitCsv(header:true)
    .map{ row-> tuple(row.group, row.id, row.sex, row.type) }
    .set { meta_gatkcov }


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
		set val(id), group, file(shard), file(shard_bai) from bwa_shards_ch.groupTuple()

	output:
		set id, group, file("${id}_merged.bam"), file("${id}_merged.bam.bai") into merged_bam, qc_merged_bam

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
		set val(group), val(id), file(r1), file(r2) from fastq

	output:
		set id, group, file("${id}_merged.bam"), file("${id}_merged.bam.bai") into bam, qc_bam

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


// Cannot be put inside process, groupTuple + any other combining operator causes an array of issues //
// All scores are collecting per sampleID, these are joined with one unique instance of their respective merged BAM //
// these are combined with all genomic shards. Each shard is run with merged bam + all scores. //
locus_collector_scores
    .groupTuple(by: [0,1])
    .join(merged_bam_id)
    .combine(dedup_shards)
    .set{ dedup_input }

// Remove duplicate reads
process dedup {
	cpus 16
	cache 'deep'
	tag "$id ($shard_name)"

	input:
		set val(id), group, file(score), file(idx), file(bam), file(bai), val(shard_name), val(shard) from dedup_input
		

	output:
		set val(id), group, file("${shard_name}_${id}.bam"), file("${shard_name}_${id}.bam.bai") into shard_dedup_bam
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
    .groupTuple(by: [0,1])
    .into{ all_dedup_bams_bqsr; all_dedup_bams_dnascope; all_dedup_bams_mergepublish }
//merge genomic shards with neighbouring shard combinations
genomicshards
    .merge(tuple(neighbour_shards))
    .into{ bqsr_shard_shard; varcall_shard_shard }

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

//Collect various QC data: TODO MOVE qc_sentieon to container!
process sentieon_qc {
	cpus 54
	memory '64 GB'
	publishDir "${OUTDIR}/qc", mode: 'copy' , overwrite: 'true'
	tag "$id"

	input:
		set id, group, file(bam), file(bai), file(dedup) from qc_bam.mix(qc_merged_bam).join(merged_dedup_metrics)

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
	echo "--run-folder $rundir --sample-id $id --subassay $diagnosis --assay wgs --qc ${OUTDIR}/qc/${id}.QC" > ${id}.cdm
	"""
}



process bqsr {
	cpus 16
	errorStrategy 'retry'
	maxErrors 5
	tag "$id ($shard_name)"

	input:
		set val(id), group, file(bams), file(bai), val(shard_name), val(shard), val(one), val(two), val(three) from all_dedup_bams_bqsr.combine(bqsr_shard_shard)

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
	publishDir "${OUTDIR}/bam", mode: 'copy', overwrite: 'true'
	tag "$id"

	input:
		set val(id), group, file(bams), file(bais) from all_dedup_bams_mergepublish

	output:
		set bgroup, id, file("${id}_merged_dedup.bam"), file("${id}_merged_dedup.bam.bai") into chanjo_bam, expansionhunter_bam, yaml_bam, cov_bam

	script:
		bams_sorted_str = bams.sort(false) { a, b -> a.getBaseName().tokenize("_")[0] as Integer <=> b.getBaseName().tokenize("_")[0] as Integer } .join(' -i ')
		bgroup = "bams"

	"""
	sentieon util merge -i ${bams_sorted_str} -o ${id}_merged_dedup.bam --mergemode 10
	"""
}


// Calculate coverage for chanjo
process chanjo_sambamba {
	cpus 16
	memory '64 GB'
	publishDir "${OUTDIR}/cov"
	tag "$id"

	input:	
		set group, id, file(bam), file(bai) from chanjo_bam.mix(chanjo_bam_choice)

	output:
		file("${id}.bwa.chanjo.cov") into chanjocov

	"""
	sambamba depth region -t ${task.cpus} -L $scoutbed -T 10 -T 15 -T 20 -T 50 -T 100 ${bam.toRealPath()} > ${id}.bwa.chanjo.cov
	"""
}




////////////////////////////////////////////////////////////////////////
////////////////////////// EXPANSION HUNTER ////////////////////////////
////////////////////////////////////////////////////////////////////////

// call STRs using ExpansionHunter
process expansionhunter {
	tag "$id"

	when:
		params.varcall
		
	input:
		set group, id, file(bam), file(bai) from expansionhunter_bam.mix(expansionhunter_bam_choice)

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

// Collect each bqsr and group by sampleID //
bqsr_merged
    .groupTuple()
    .set{ bqsr_merged_grouped }

// join the bqsr to all genomic sharded bams //
// combine all with each shard and their neighbouring shards //
all_dedup_bams_dnascope
    .join(bqsr_merged_grouped)
    .combine(varcall_shard_shard)
    .set{ varcall_shard_shard_bams }

// Do variant calling using DNAscope, sharded
process dnascope {
	cpus 16
	tag "$id ($shard_name)"

	when:
		params.varcall

	input:
		set id, group, file(bams), file(bai), file(bqsr), val(shard_name), val(shard), val(one), val(two), val(three) from varcall_shard_shard_bams
		
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
		set vgroup, ph, file("${id}.dnascope.gvcf.gz"), file("${id}.dnascope.gvcf.gz.tbi") into complete_vcf
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

complete_vcf
	.mix(complete_vcf_choice)
	.mix(gvcf_choice)
    .groupTuple()
    .set{ gvcfs }

process gvcf_combine {
    cpus 16
	tag "$group"

    input:
	set vgroup, ph, file(vcf), file(idx) from gvcfs
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


//madeline ped, run if family mode
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
	tag "$group"

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
		set group, file("${group}.intersected.vcf") into split_vep, split_cadd, vcf_loqus

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
	echo "loqusdb load -f ${ped.toRealPath()} --variant-file ${vcf.toRealPath()}" > ${group}.loqus
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

// Annotating variants with clinvar
process annotate_clinvar {
    cpus 1
    memory '32GB'
	tag "$group"

	input:
		set group, file(vcf) from vep

	output:
		set group, file("${group}.clinvar.vcf") into snpsift

	"""
	SnpSift -Xmx60g annotate $params.CLINVAR \\
		-info CLNSIG,CLNACC,CLNREVSTAT $vcf > ${group}.clinvar.vcf
	"""

}

// Annotating variants with Genmod
process annotate_genmod {
	cpus 2
	tag "$group"

	input:
		set group, file(vcf) from snpsift

	output:
		set group, file("${group}.genmod.vcf") into genmod

	"""
	genmod annotate --genome-build 38 --annotate_regions $vcf -o ${group}.genmod.vcf
	"""
}

// # Annotating variant inheritance models:
process inher_models {
	cpus 6
	memory '64 GB'
	tag "$group"

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
process loqdb {
	cpus 1
	queue 'bigmem'
	errorStrategy 'retry'
	maxErrors 5
	tag "$group"

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
	tag "$group"

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
	publishDir "${OUTDIR}/vcf", mode: 'copy', overwrite: 'true'
	tag "$group"

	input:
		set group, file(vcf) from scored_vcf

	output:
		set group, file("${group}.scored.vcf.gz"), file("${group}.scored.vcf.gz.tbi") into vcf_peddy, vcf_yaml

	"""
	bgzip -@ ${task.cpus} $vcf -f
	tabix ${vcf}.gz -f
	"""
}


// Running PEDDY: 
process peddy {
	publishDir "${OUTDIR}/ped", mode: 'copy' , overwrite: 'true'
	cpus 6
	tag "$group"

	input:
		file(ped) from ped_peddy
		set group, file(vcf), file(idx) from vcf_peddy

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
	tag "$group"

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
	tag "$group"

	input:
		set gr, file(vcf) from vcf_upd
		set group, id, sex, mother, father, phenotype, diagnosis from meta_upd

	output:
		file("upd.bed") into upd_plot
		set group, file("upd.sites.bed") into upd_table

	when:
		mode == "family" && trio == true

	"""
	upd --vcf $vcf --proband $id --mother $mother --father $father --af-tag GNOMADAF regions > upd.bed
	upd --vcf $vcf --proband $id --mother $mother --father $father --af-tag GNOMADAF sites > upd.sites.bed
	"""
}


process upd_table {
	publishDir "${OUTDIR}/plots", mode: 'copy' , overwrite: 'true'
	tag "$group"

	input:
		set group, file(upd_sites) from upd_table

	output:
		file("${group}.UPDtable.xls")

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
	memory '16 GB'

	input:
		set id, group, file(bam), file(bai), gr, sex, type from cov_bam.join(meta_gatkcov, by:1)

	output:
		set group, id, type, file("${id}.standardizedCR.tsv"), file("${id}.denoisedCR.tsv") into cov_plot, cov_gens

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
	tag "$group"

	input:
		file(upd) from upd_plot
		set gr, file(roh) from roh_plot
		set group, id, type, file(cov_stand), file(cov_denoised) from cov_plot.groupTuple()


	output:
		file("${id[proband_idx]}.genomic_overview.png")

	script:
		proband_idx = type.findIndexOf{ it == "proband" }

	"""
	genome_plotter.pl --dict $params.GENOMEDICT \\
		 --sample ${id[proband_idx]} \\
		 --upd $upd \\
		 --roh $roh \\
		 --cov ${cov_denoised[proband_idx]} \\
		 --out ${id[proband_idx]}.genomic_overview.png
	"""
}

process generate_gens_data {
	publishDir "${OUTDIR}/cov", mode: 'copy' , overwrite: 'true'
	tag "$group"
	cpus 1

	input:
		set group, id, file(gvcf), type, file(cov_stand), file(cov_denoise) from gvcf_gens.join(cov_gens, by:[1]).view()

	output:
		set file("${id}.cov.bed.gz"), file("${id}.baf.bed.gz")

	"""
	generate_gens_data.pl $cov_stand $gvcf $id $params.GENS_GNOMAD
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
		set group, id, file(bam), file(bai) from yaml_bam.groupTuple()
		set group, file(vcf), file(idx) from vcf_yaml
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
	create_yml.pl \\
		$bams \\
		$group \\
		$OUTDIR \\
		$diagnosis \\
		PORT_CMDSCOUT2_MONGODB \\
		> ${group}.yaml
	"""
}
