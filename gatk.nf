#!/usr/bin/env nextflow

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

process gatk_coverage {
    cpus 10
    memory '20GB'
    time '2h'
    container = '/fs1/resources/containers/gatk_4.1.9.0.sif'
    scratch true
	stageInMode 'copy'
	stageOutMode 'copy'

    input:
        set group, id, file(bam) from bam_gatk

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

    input:
        set group, id, file(tsv) from call_ploidy.mix(call_ploidy_choice)

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

    input:
        set group, id, file(tsv), file(ploidy), i, refpart \
            from call_cnv.mix(call_cnv_choice).join(ploidy_to_cnvcall, by: [0,1]).combine(gatk_ref)

    output:
        set group, id, i, file("${group}_${i}.tar") into postprocessgatk

    """
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


    input:
        set group, id, i, file(tar), file(ploidy) from postprocessgatk.groupTuple(by: [0,1]).join(ploidy_to_post, by: [0,1])
        set shard_no, shard from gatk_postprocess.groupTuple(by: [3])

    output:
        set group, id, \
            file("genotyped-intervals-${group}-vs-cohort30.vcf.gz"), \
            file("genotyped-segments-${group}-vs-cohort30.vcf.gz"), \
            file("denoised-${group}-vs-cohort30.vcf.gz")

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