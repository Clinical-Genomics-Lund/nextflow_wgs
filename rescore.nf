
// Input channels for various meta information //
Channel
	.fromPath(params.csv)
	.splitCsv(header:true)
	.map{ row-> tuple(row.group, file(row.sv_vcf), file(row.snv_vcf) ) }
	.into{ input; madde_group; vcf_peddy; yaml_input }

csv = file(params.csv)
base = csv.getParent()
name = csv.getSimpleName()
ped_file = base/name+'.ped'
ped = file( ped_file )

Channel
	.fromPath(params.csv)
	.splitCsv(header:true)
	.map{ row-> tuple(row.group, file(row.str_vcf) ) }
	.set{ str_input }

Channel
    .fromPath(file(ped))
    .into{ strip_ped; ped_compound; ped_inher; ped_mad; ped_peddy; ped_scout; ped_bam }

OUTDIR = params.outdir+'/'+params.subdir
CRONDIR = params.crondir


mode = ped.countLines() > 1 ? "family" : "single"
println(mode)

process bam_to_info {
	cpus 1
	time '2m'
	memory '2 MB'

	input:
		set ped from ped_bam

	output:
		file("${pgroup}.INFO") into bam_INFO

	script:
		id = 'ph'
		pgroup = 'ph'
		bamrows = ''
		ped.readLines().each{
			if (it =~ /^\S+\t+\S+\t+\S+\t+\S+\t+\S+\t+\S+/) {
				id = it =~ /^(\S+)\t+(\S+)\t+\S+\t+\S+\t+\S+\t+\S+/
				pgroup = id[0][1]
				row = 'BAM	'+id[0][2]+' /access/'+params.subdir+'/bam/'+id[0][2]+'_merged_dedup.bam¤'
				bamrows = bamrows + row
			}
		}
	"""
	echo $bamrows | sed "s/¤/\\n/g" > ${pgroup}.INFO
	"""


}

process str_to_info {
    cpus 1
    time '5m'
    memory '100 MB'

    input:
    	set group, file(str_vcf) from str_input
    
    output:
    	file("${group}.INFO") into str_INFO

    """
	echo "STR	${OUTDIR}/vcf/${str_vcf}" > ${group}.INFO
    """
}

process strip_score {
    cpus 2
    time '30m'
    memory '5 GB'

    input:
        set group, file(sv_vcf), file(snv_vcf) from input
        file(ped) from strip_ped
    
    output:
        set group, file("${group}.sv_stripped.vcf") into stripped_vcf
		set group, file("${group}.snv_stripped.vcf") into stripped_vcf_snv
    """
    rescore_vcf.pl --vcf $sv_vcf --ped $ped > ${group}.sv_stripped.vcf &
	rescore_vcf.pl --vcf $snv_vcf --ped $ped > ${group}.snv_stripped.vcf
    """
}

// # Annotating variant inheritance models:
process inher_models {
	cpus 6
	memory '64 GB'
	tag "$group"
	time '10m'

	input:
		set group, file(vcf) from stripped_vcf_snv
		file(ped) from ped_inher

	output:
		set group, file("${group}.models.vcf") into inhermod

	"""
	genmod models $vcf -p ${task.cpus} -f $ped > ${group}.models.vcf
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
		set group, file(vcf) from inhermod

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
	publishDir "${OUTDIR}/vcf/rescore", mode: 'copy', overwrite: 'true', pattern: '*.vcf.gz*'
	tag "$group"

	input:
		set group, file(vcf) from scored_vcf

	output:
		set group, file("${group}.scored.vcf.gz"), file("${group}.scored.vcf.gz.tbi") into snv_sv_vcf
		file("${group}.INFO") into snv_INFO

	"""
	bgzip -@ ${task.cpus} $vcf -f
	tabix ${vcf}.gz -f
	echo "SNV	${OUTDIR}/vcf/rescore/${group}.scored.vcf.gz" > ${group}.INFO
	"""
}

process score_sv {
	cpus 5
	tag "$group $mode"
	publishDir "${OUTDIR}/vcf/rescore", mode: 'copy', overwrite: 'true', pattern: '*.vcf.gz*'
	memory '10 GB'
	time '2h'

	input:
		set group, file(vcf) from stripped_vcf

	output:
		set group, file("${group}.sv.scored.sorted.vcf.gz"), file("${group}.sv.scored.sorted.vcf.gz.tbi") into sv_rescore
		file("${group}.INFO") into sv_INFO
				
	script:
	
		if (mode == "family") {
			"""
			genmod score -i $group -c $params.svrank_model -r $vcf -o ${group}.sv.scored_tmp.vcf
			bcftools sort -O v -o ${group}.sv.scored.sorted.vcf ${group}.sv.scored_tmp.vcf 
			bgzip -@ ${task.cpus} ${group}.sv.scored.sorted.vcf -f
			tabix ${group}.sv.scored.sorted.vcf.gz -f
			echo "SV	${OUTDIR}/vcf/rescore/${group}.sv.scored.sorted.vcf.gz" > ${group}.INFO
			"""
		}
		else {
			"""
			genmod score -i $group -c $params.svrank_model_s -r $vcf -o ${group}.sv.scored.vcf
			bcftools sort -O v -o ${group}.sv.scored.sorted.vcf ${group}.sv.scored.vcf
			bgzip -@ ${task.cpus} ${group}.sv.scored.sorted.vcf -f
			tabix ${group}.sv.scored.sorted.vcf.gz -f
			echo "SV	${OUTDIR}/vcf/rescore/${group}.sv.scored.sorted.vcf.gz" > ${group}.INFO
			"""
		}
}

process compound_finder {
	cpus 5
	tag "$group $mode"
	publishDir "${OUTDIR}/vcf/rescore", mode: 'copy', overwrite: 'true', pattern: '*.vcf.gz*'
	memory '10 GB'
	time '2h'

	when:
		mode == "family"

	input:
		set group, file(vcf), file(tbi) from sv_rescore
		file(ped) from ped_compound
		set group2, file(snv), file(tbi) from snv_sv_vcf
		

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
		echo "SVc	${OUTDIR}/vcf/rescore/${group}.sv.rescored.sorted.vcf.gz,${OUTDIR}/vcf/rescore/${group}.snv.rescored.sorted.vcf.gz" > ${group}.INFO
		"""

}

//madeline ped, run if family mode
process madeline {
	publishDir "${OUTDIR}/ped/rescore", mode: 'copy' , overwrite: 'true', pattern: '*.xml'
	memory '1 GB'
	time '30m'


	input:
		file(ped) from ped_mad
		set group, file(sv_vcf), file(snv_vcf) from madde_group

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
	echo "MADDE ${OUTDIR}/ped/rescore/${ped}.madeline.xml" > ${group}.INFO
	"""
}

// Running PEDDY: 
process peddy {
	publishDir "${OUTDIR}/ped/rescore", mode: 'copy' , overwrite: 'true'
	container = '/fs1/resources/containers/wgs_20200115.sif'
	cpus 6
	tag "$group"

	input:
		file(ped) from ped_peddy
		set group, file(sv_vcf), file(vcf) from vcf_peddy

	output:
		set file("${group}.ped_check.csv"),file("${group}.peddy.ped"), file("${group}.sex_check.csv") into peddy_files
		file("${group}.INFO") into peddy_INFO

	"""
	source activate py3-env
	python -m peddy --sites hg38 -p ${task.cpus} ${vcf.toRealPath()} $ped --prefix $group
	echo "PEDDY	${OUTDIR}/ped/rescore/${group}.ped_check.csv,${OUTDIR}/ped/rescore/${group}.peddy.ped,${OUTDIR}/ped/rescore/${group}.sex_check.csv" > ${group}.INFO
	"""
}

// Collects $group.INFO files from each process output that should be included in the yaml for scout loading //
// If a new process needs to be added to yaml. It needs to follow this procedure, as well as be handled in create_yml.pl //
bam_INFO
	.mix(snv_INFO,sv_INFO,peddy_INFO,madde_INFO,svcompound_INFO,str_INFO)
	.collectFile()
	.set{ yaml_INFO }

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
		file(INFO) from yaml_INFO
		file(ped) from ped_scout
		set group, file(sv_vcf), file(snv_vcf) from yaml_input

	output:
		set group, file("${group}.yaml") into yaml_output

	script:

	"""
	export PORT_CMDSCOUT2_MONGODB=33002 #TA BORT VÄLDIGT FULT
	create_yml.pl \\
		--g $group --d other --panelsdef $params.panelsdef --out ${group}.yaml --ped $ped --files $INFO --assay wgs-hg38 --antype $params.antype
	"""
}

// export from scout-db needs DIAGNOSIS(genlista), clarity-lims-id, assay, bam-file-path? vcf-path, new groupname 