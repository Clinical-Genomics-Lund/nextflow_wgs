When running the pipeline, two separate results folders are produced. One is named "OUTDIR", and contains the bulk of the results. The second is named "CRONDIR" and contains specific files that are transferred away from the computational cluster in the Lund setup and used to load the files into Scout. It is called "CRONDIR" as the transfer from the computational cluster is triggered by a cron-job.

Here the output files, their origin and their meaning are very briefly outlined.

For a high-level overview of the processes, see the overview visualizatoin in the [running the pipeline](2_running_the_pipeline.md) section.

### OUTDIR

| Subfolder      | Function    | Process                | Name                                              | Description                
|----------------|-------------|------------------------|---------------------------------------------------|----------------------------
| annotsv        | Annotation  | `annotsv`              | `<group>_annotsv.tsv`                             | Annotation from AnnotSV    
| bam            | Mapping     | `markdup`              | `<group>_dedup.bam[.bai]`                         | Deduplicated BAM-files     
| bam            | Mapping     | `fetch_MTseqs`         | `<group>_mito.bam[.bai]`                          | Mapped mitochondrial reads  
| bed            | SV calling  | `svvcf_to_bed`         | `<group>.sv.bed`                                  | SV ranges                  
| bqsr           | SNV calling | `bqsr`                 | `<group>.bqsr.table`                              | Base-quality score recalibrations for Sentieon DNAScope
| cov            | Coverage    | `chanjo_sambamba`      | `<group>.bwa.chanjo.cov`                          | Coverage calculations
| cov            | Coverage    | `depth_onco`           | `<group>.lowcov.overlapping.bed`                  | Coverage calculations (onco only)
| cov            | Coverage    | `gatkcov`              | `<group>.denoisedCR.tsv`                          | Coverage calculations (wgs only)
| cov            | Coverage    | `gatkcov`              | `<group>.standardizedCR.tsv`                      | Coverage calculations (wgs only)
| ped            | Pedigree    | `create_ped`           | `<group>_[base|ma|fa].ped`                        | Pedigree information
| ped            | Pedigree    | `peddy`                | `<group>_[base|ma|fa].peddy.ped`                  | Pedigree information (FIXME)
| ped            | Pedigree    | `madeline`             | `<group>_[base|ma|fa].ped.madeline.xml`           | Pedigree visualization files
| ped            | Pedigree    | `FIXME`                | `<group>.sex_check.csv`                           | Control of assigned sex
| plot_data      | Visualize   | `generate_gens_data`   | Various files                                     | Files needed for visualization in Gens
| plots/SMNcnc   | Visualize   | `SMNCopyNumberCaller`  | `smn_<group>.pdf`                                 | Copy number reports for SMN1 and SMN2
| plots/reviewer | Visualize   | `reviewer`             | `<group>/<protein_ids>.svg`                       | Visualizations of repeat expansions
| plots/mito     | Visualize   | `run_haplogrep`        | `<group>.haplogrep.png`                           | Mitochondria haplogroup overview
| plots/mito     | Visualize   | `run_eklipse`          | `<group>_eklipse.png`                             | Mitochondrial deletions visualized
| plots          | Visualize   | `upd_table`            | `<group>.UPDTable.xls`                            | Uniparental disomy overview
| overview_plot  | Visualize   | `overview_plot`        | `<group>.genomic_overview.png`                    | High-level illustration of uniparental disomy and runs of homozygosity
| qc             | QC          | `merge_qc_json`        | `<group>.QC`                                      | Sample QC information
| pod            | Visualize   | `plot_pod`             | `<group>_POD_karyotype.pdf`                       | CNV origin illustration (for trio)
| pod            | Visualize   | `plot_pod`             | `<group>_POD_results.html`                        | CNV origin illustration (for trio)
| sv_vcf         | SV calling  | `postprocessgatk`      | `denoised-<group>-vs-cohort30.vcf.gz`             |
| sv_vcf         | SV calling  | `postprocessgatk`      | `genotyped-intervals-<group>-vs-cohort30.vcf.gz`  |
| sv_vcf         | SV calling  | `postprocessgatk`      | `genotyped-segments-<group>-vs-cohort30.vcf.gz`   |
| sv_vcf         | SV calling  | `filter_merge_gatk`    | `<group>.gatk.filtered.merged.vcf`                |
| sv_vcf         | SV calling  | `manta`                | `<group>.manta.vcf.gz`                            |
| sv_vcf         | SV calling  | `manta_panel`          | `<group>.manta.vcf.gz`                            |
| sv_vcf         | SV calling  | `delly_panel`          | `<group>.delly.vcf.gz`                            |
| sv_vcf         | SV calling  | `cnvkit_panel`         | `<group>.cnvkit_filtered.vcf`                     |
| sv_vcf         | SV calling  | `svdb_merge_panel`     | `merged/<group>.merged.filtered.melt.vcf`         |
| sv_vcf         | SV calling  | `tiddit`               | `<group>.tiddit.filtered.vcf`                     |
| sv_vcf         | SV calling  | `svdb_merge`           | `merged/<group>.merged.vcf`                       |
| sv_vcf         | SV calling  | `svdb_merge`           | `merged/<group>.merged.bndless.vcf`               |

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