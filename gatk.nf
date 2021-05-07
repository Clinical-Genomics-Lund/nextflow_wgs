#!/usr/bin/env nextflow
// GENERAL PATHS //
OUTDIR = "/fs1/results/wgs/"
OUTDIR_RECALL = "/fs1/viktor/wgs/gatk_recall/sv_vcf"

Channel
	.fromPath(params.gatkreffolders)
	.splitCsv(header:true)
	.map{ row-> tuple(row.i, row.refpart) }
	.into{ gatk_ref; gatk_postprocess }

Channel
	.fromPath(params.csv)
	.splitCsv(header:true)
	.map{ row-> tuple(row.group, row.id, file(row.bam)) }
	.set{ input_files }

// for testing purposes, tsv as input
tsv_choice = Channel.create()
bam_gatk = Channel.create()
input_files.view().choice(bam_gatk, tsv_choice) { it[2]  =~ /\.tsv/ ? 1 : 0  }

tsv_choice.into{ call_ploidy_choice; call_cnv_choice }

mode = "single"

process gatk_coverage {
    cpus 10
    memory '20GB'
    time '2h'
    container = '/fs1/resources/containers/gatk_4.1.9.0.sif'
    scratch true
	stageInMode 'copy'
	stageOutMode 'copy'
    tag "$id"   

	when:
		params.sv && !params.onco  &&  !params.exome

    input:
        set group, id, file(bam), file(bai) from bam_gatk

    output:
        set group, id, file("${id}.tsv") into call_ploidy, call_cnv

    """
    export MKL_NUM_THREADS=${task.cpus}
    export OMP_NUM_THREADS=${task.cpus}
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
    memory '20GB'
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
    export MKL_NUM_THREADS=${task.cpus}
    export OMP_NUM_THREADS=${task.cpus}
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
    memory '25GB'
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
	export HOME=/local/scratch
	echo "[global]" > ~/.theanorc
	echo config.compile.timeout = 1000 >> ~/.theanorc
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
    memory '25GB'
    time '3h'
    container = '/fs1/resources/containers/gatk_4.1.9.0.sif'
    scratch true
	stageInMode 'copy'
	stageOutMode 'copy'
    //publishDir "${OUTDIR}/sv_vcf/", mode: 'copy', overwrite: 'true'
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
    export MKL_NUM_THREADS=!{task.cpus}
    export OMP_NUM_THREADS=!{task.cpus}
	for model in !{tar}; do
	tar -xvf $model
	done
    tar -xvf !{ploidy}
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
    publishDir "${OUTDIR_RECALL}", mode: 'copy', overwrite: 'true'

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
	publishDir "${OUTDIR_RECALL}/merged/", mode: 'copy', overwrite: 'true'
	time '2h'
	memory '1 GB'
	scratch true
	stageInMode 'copy'
	stageOutMode 'copy'

	input:
		// set group, id, file(mantaV) from called_manta.groupTuple()
		// set group, id, file(tidditV) from called_tiddit.groupTuple()
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
            // mantaV = file("${OUTDIR}/sv_vcf/${id[0]}.manta.vcf.gz")
            // tidditV = file("${OUTDIR}/sv_vcf/${id[0]}.tiddit.filtered.vcf")
			// tmp = mantaV.collect {it + ':manta ' } + tidditV.collect {it + ':tiddit ' } + gatkV.collect {it + ':gatk ' }
			// vcfs = tmp.join(' ')
			"""
			source activate py3-env
			svdb --merge --vcf ${OUTDIR}/sv_vcf/${id[0]}.manta.vcf.gz:manta ${OUTDIR}/sv_vcf/${id[0]}.tiddit.filtered.vcf:tiddit $gatkV:gatk --no_intra --pass_only --bnd_distance 2500 --overlap 0.7 --priority manta,tiddit,gatk > ${group}.merged.vcf
			grep -v BND ${group}.merged.vcf > ${group}.merged.bndless.vcf
			"""
		}

}

