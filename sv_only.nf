
csv = file(params.csv)
ped = file(params.ped)
mode = "family"
Channel
	.fromPath(params.csv)
	.splitCsv(header:true)
	.map{ row-> tuple(row.group, row.id, file(row.gatk)) }
	.set{ called_gatk }

Channel
	.fromPath(params.csv)
	.splitCsv(header:true)
	.map{ row-> tuple(row.group, row.id, file(row.manta)) }
	.set{ called_manta }

Channel
	.fromPath(params.csv)
	.splitCsv(header:true)
	.map{ row-> tuple(row.group, row.id, file(row.tiddit)) }
	.set{ called_tiddit }

process merge_gatk {
	cpus 1
	tag "$group"
	time '2h'
	memory '1 GB'

	input:
		set group, id, file(gatk) from called_gatk

	output:
		set group, id, file("${id}.gatk.merged.vcf") into merged_gatk

	"""
	mergeGATK.pl $gatk > ${id}.gatk.merged.vcf
	"""
}

process svdb_merge {
	cpus 1
	tag "$group"
	publishDir "/fs1/viktor/gatk_ref/called_vcfs/merged/", mode: 'copy', overwrite: 'true'
	time '2h'
	memory '1 GB'
	// scratch true
	// stageInMode 'copy'
	// stageOutMode 'copy'

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
			tmp = mantaV.collect {it + ':manta ' } + tidditV.collect {it + ':tiddit ' } + gatk.collect {it + ':gatk ' }
			vcfs = tmp.join(' ')
			"""
			source activate py3-env
			svdb --merge --vcf $vcfs --no_intra --pass_only --bnd_distance 2500 --overlap 0.7 --priority manta,tiddit,gatk > ${group}.merged.vcf
			grep -v BND ${group}.merged.vcf > ${group}.merged.bndless.vcf
			"""
		}

}

//create AnnotSV tsv file
process annotsv {
	container = '/fs1/resources/containers/annotsv.v2.3.sif'
	cpus 2
	tag "$group"
	//publishDir "${OUTDIR}/annotsv/", mode: 'copy', overwrite: 'true'
	time '5h'
	memory '20 GB'

	input:
		set group, id, file(sv) from annotsv_vcf
		
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
	memory '150 GB'
	time '1h'
	
	input:
		set group, id, file(vcf) from vcf_vep

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
		set group, file(sv_artefact) from manip_vcf
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
	cache false
	tag "$group"
	publishDir "/fs1/viktor/gatk_ref/called_vcfs/scored", mode: 'copy', overwrite: 'true', pattern: '*.vcf*'
	memory '10 GB'
	time '2h'

	input:
		set group, file(vcf) from annotatedSV

	output:
		set group, file("${group}.sv.scored.vcf"), file("${group}.sv.scored.sorted.vcf.gz.tbi") into sv_rescore
				
	script:
	"""
	genmod score -i $group -c $params.svrank_model_s -r $vcf -o ${group}.sv.scored.vcf
	bcftools sort -O v -o ${group}.sv.scored.sorted.vcf ${group}.sv.scored.vcf
	bgzip -@ ${task.cpus} ${group}.sv.scored.sorted.vcf -f
	tabix ${group}.sv.scored.sorted.vcf.gz -f
	"""

}





