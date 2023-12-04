### OUTDIR

Brief overview of the various output files that are produced by the WGS-pipeline.

Some outputs are only produced for certain run settings (i.e. profile and trio/single mode).

An overview of the different processes is found in the main page of this README.

* annotsv
    * annotsv: AnnotSV annotations (.tsv)
* bam
    * markdup: Deduplicated BAM files with index (.bam + .bam.bai)
    * fetch_MTseqs: Mitochondria extracted from bam (.bam + .bam.bai)
* bed
    * svvcf_to_bed: BED file with SVs (.bed)
* bqsr
    * bqsr: Base quality recalibration scores (.table, text files)
* cov
    * chanjo_sambamba: Coverage calculations (.cov, text files)
    * depth_onco: Onco coverage
    * gatkcov: Coverage calculations from GATK (.tsv) <- FIXME: What is different with this coverage?
* madeline
    * madeline: Pedigree illustration (.xml)
* ped
    * create_ped: Pedigree file (.ped)
    * peddy: Pedigree fiels (.ped, .csv) <- FIXME: Double check
* plot_data
    * generate_gens_data: (FIXME: Check)
* plots
    * SMNCopyNumberCaller, /SMNcnc: Illustrations of copy numbers for SMN exons (.pdf)
    * reviewer, /reviewer: Illustrations of repeat expansions
    * run_haplogrep, /mito: Illustration of subtypes of mitochondria (? FIXME)
    * run_eklipse, /mito: Illustration of coverage in different parts of mitochondria (.png + .txt)
    * upd_table: Table with uniparental disomy information (.xls)
    * overview_plot: High-level overview of chromosome differences - uniparental disomy and homozygotic runs (FIXME: Check)
    * cnvkit_panel: Visualization of CNVKit calls (FIXME: Check)
* qc
    * merge_qc_json: QC information for samples (.json)
* pod
    * plot_pod: Parental origin of duplication (.pdf)
* sv_vcf
    * postprocessgatk: GATK-called CNV segments (.vcf.gz)
    * filter_merge_gatk: Further processed GATK segments, merging adjacent and cleaning up reference segments (.vcf.gz)
    * manta: SV calls from manta (.vcf.gz)
    * manta_panel: SV calls from manta for panel (.vcf.gz)
    * delly_panel: SV calls from Delly for panel (.vcf.gz)
    * cnvkit_panel: SV calls from CNVKit for panel (.vcf.gz)
    * svdb_merge_panel, /merged: Combine panel VCFs into one single
    * tiddit: SV calls from TIDDIT (.vcf)
    * svdb_merge: Combined SV calls from manta, tiddit and GATK CNVCaller (.vcf)
* vcf
    * vcfbreakmulti_expansionhunter: Variant calls from expansionhunter (.vcf)
    * melt: Mobile elements (.vcf)
    * intersect_melt: Mobile elements within regions of interest (.vcf)
    * run_mutect2: Mitochondria variants from Mutect2 (.vcf)
    * split_normalize: Variant calls combined with hmtnote output, and with regions of interest extracted
    * vcf_completion: Scored variant calls (.vcf.gz + .vcf.gz.tbi)
    * fastgnomad: Filtered VCF removing high-frequency variants
    * score_sv: Scored SV (.vcf.gz + .vcf.gz.tbi)
    * compound_finder: SNV and SV updated scoring based on compounding (multiple variants in same genes at different locations) (.vcf.gz + .vcf.gz.tbi)
* versions
    * combine_versions: Combined versions from all involved processes (.yml)
* yaml
    * create_yaml: Summarized key information for Scout (.yml) (FIXME: What is the `alt_effect` output?)

### CRONDIR

Transferred to a separate computer for further processing in CMD (Lennart).

* gens
    * generate_gens_data: Prepared output data for usage in Gens (FIXME: Check format)
* loqus
    * LoqusDB variants (artefacts, .locus format)
* qc
    * qc_to_cdm: Samples quality information for CDM
* scout
    * create_yaml: File used to specify loading locations for Scout (.yml)