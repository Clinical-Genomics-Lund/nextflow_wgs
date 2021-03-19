#!/usr/bin/env nextflow
genome_file = params.genome_file
Channel
	.fromPath(params.csv)
	.splitCsv(header:true)
	.map{ row-> tuple(row.group, row.id, file(row.bam), file(row.bai)) }
	.set { bam_ch }

Channel
	.fromPath(params.csv)
	.splitCsv(header:true)
	.map{ row-> tuple(row.group, row.id, file(row.vcf)) }
	.set { diploid_vcf }

// create an MT BAM file
process fetch_MTseqs {

    input:
        set group, id, file(bam), file(bai) from bam_ch

    output:
        set group, id, file ("${id}_mito.bam"), file("${id}_mito.bam.bai") into mutserve_bam, eklipse_bam

    """
    sambamba view -f bam $bam M > ${id}_mito.bam
    samtools index -b ${id}_mito.bam
    """

}

// gatk FilterMutectCalls in future if FPs overwhelms tord/sofie/carro
process run_mutect2 {
    cpus 4
    memory '16 GB'
    time '15m'
    
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
process run_bcftools {
    cpus 1
    memory '1GB'
    time '5m'

    input:
        set group, id, file(ms_vcf) from ms_vcfs_1

    output:
        set group, id, file("${ms_vcf.baseName}.adjusted.vcf") into adj_vcfs


    """
    vcfbreakmulti $ms_vcf > ${ms_vcf}.breakmulti
    bcftools sort ${ms_vcf}.breakmulti | bgzip > ${ms_vcf}.breakmulti.fix
    tabix -p vcf ${ms_vcf}.breakmulti.fix
    bcftools norm -f $params.rCRS_fasta -o ${ms_vcf.baseName}.adjusted.vcf ${ms_vcf}.breakmulti.fix    
    """

}


// use python tool HmtNote for annotating vcf
// future merging with diploid genome does not approve spaces in info-string
process run_hmtnote {
    cpus 1
    memory '5GB'
    time '15m'


    input:
        set group, id, file(adj_vcf) from adj_vcfs
        set group, id, file(vcf) from diploid_vcf.first()

    output:
        set group, file("${group}.concatmito.vcf") into hmtnote_vcfs

    """
    source activate tools
    hmtnote annotate ${adj_vcf} ${group}.hmtnote --offline
    grep ^# ${group}.hmtnote > ${group}.fixinfo.vcf
    grep -v ^# ${group}.hmtnote | sed 's/ /_/g' >> ${group}.fixinfo.vcf
    java -jar /opt/conda/envs/CMD-WGS/share/picard-2.21.2-1/picard.jar MergeVcfs \
        I=$vcf I=${group}.fixinfo.vcf O=${group}.concatmito.vcf
    """
    
}


process annotate_vep {
	container = '/fs1/resources/containers/ensembl-vep_latest.sif'
	cpus 54
	tag "$group"
	memory '20 GB'
	time '5h'
	scratch true
	stageInMode 'copy'
	stageOutMode 'copy'

	input:
		set group, file(vcf) from hmtnote_vcfs

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

// run haplogrep 2 on resulting vcf
process run_haplogrep {
    time '10m'
    memory '16 GB'
    cpus '2'

    input:
        set group, id, file(ms_vcf) from ms_vcfs_2

    output:
       // set group, id, file(".hg2.vcf.dot") into hg2_dots  Merge med montage efter for loop!!

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
    montage -mode concatenate -tile 3x1 *.png !{group}.png
    '''

}

// use eKLIPse for detecting mitochondrial deletions
process run_eklipse {
    cpus 2
    memory '10GB'
    time '20m'

    input:
        set group, id, file(bam), file(bai) from eklipse_bam


    """
    source activate htslib10
    echo "${bam}\tsample" > infile.txt
    python /eKLIPse/eKLIPse.py \
    -in infile.txt \
    -ref /eKLIPse/data/NC_012920.1.gb
    mv eKLIPse_*/eKLIPse_deletions.csv ./${id}_deletions.csv
    mv eKLIPse_*/eKLIPse_genes.csv ./${id}_genes.csv
    mv eKLIPse_*/eKLIPse_sample.png ./${id}_plot.png
    hetplasmid_frequency_eKLIPse.pl --bam ${bam} --in ${id}_deletions.csv
    """

}
