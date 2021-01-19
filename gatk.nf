#!/usr/bin/env nextflow

Channel
	.fromPath(params.gatkreffolders)
	.splitCsv(header:true)
	.map{ row-> tuple(row.i, row.refpart) }
	.set{ gatk_ref }

Channel
	.fromPath(params.csv)
	.splitCsv(header:true)
	.map{ row-> tuple(row.group, row.id, file(row.tsv)) }
	.into{ call_ploidy; called_cnvkit_panel }

// process gatk_coverage {
//     cpus 10
//     memory '20GB'
//     time '2h'
//     container = '/fs1/resources/containers/gatk_4.1.9.0.sif'

//     input:
//         set group, id, file(bam) from somewhere

//     output:
//         set group, id, file("${id}.tsv") into call_ploidy, called_cnvkit_panel

//     """
//     gatk CollectReadCounts \\
//         -L grch38.preprocessed.blacklisted.gcfiltered.interval_list \\
//         -R /fs1/resources/ref/hg38/fasta/GCA_000001405.15_GRCh38_no_alt_analysis_set_nochr.fna \\
//         -imr OVERLAPPING_ONLY \\
//         -I $bam \\
//         --format TSV -O ${i}.tsv
//     """
// }

process gatk_call_ploidy {
    cpus 10
    memory '20GB'
    time '2h'
    container = '/fs1/resources/containers/gatk_4.1.9.0.sif'

    input:
        set group, id, file(tsv) from call_ploidy

    output:
        set group, id, file("ploidy.tar") into ploidy_to_cnvcall

    """
    gatk DetermineGermlineContigPloidy \\
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
            from called_cnvkit_panel.join(ploidy_to_cnvcall, by: [0,1]).combine(gatk_ref)

    """
    export MKL_NUM_THREADS=${task.cpus}
    export OMP_NUM_THREADS=${task.cpus}
    tar -xvf ploidy.tar
    mkdir $group
    gatk --java-options "-Xmx25g" GermlineCNVCaller \\
        --run-mode CASE \\
        -I $tsv \\
        --contig-ploidy-calls ploidy/${group}-calls/ \\
        --model ${refpart} \\
        --output ${group}/ \\
        --output-prefix ${group}_${i}
    """
}