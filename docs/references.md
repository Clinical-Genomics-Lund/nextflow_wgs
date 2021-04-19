## This is the reference docs

* Reference genome, hg38
* genome shards (in the repo)
* SNV-CADD (https://kircherlab.bihealth.org/download/CADD/v1.6/GRCh38/whole_genome_SNVs.tsv.gz)
* VEP-cache
* gnomAD-exomes (vcf)
* gnomAD-genomes (vcf)
* PhyloP (bigwig)
* PhastCons (bigwig)
* Mills_and_1000G_gold_standard.indels (vcf)
* ClinVar-vcf, keep this updated weekly
* Intersect-bed to reduce the amount of SNVs (basicly exome + interesting regions)
* GATK CNV bins (interval list) - ADD TO DOC!
* GATK CNV PON, males (hdf5) - ADD TO DOC!
* GATK CNV PON, females (hdf5) - ADD TO DOC!
* Fastgnomad SNP-list (dat) - ADD TO DOC!
* Gens SNP-list (txt) - ADD TO DOC!
* SVDB-databas (db) - might be going away, (https://github.com/J35P312/SVDB)
* loqusdb export (vcf) - Local observations, both artefacts and real variants. (https://github.com/moonso/loqusdb)
* Ensembl genes (bed)
* Gene panel-dump (json) - Needed for last step to create yaml for scout.
* Rankmodels
	* [SNV-single](https://github.com/Clinical-Genomics-Lund/nextflow_wgs/blob/documentation/rank_models/rank_model_v5.01_single.ini)
	* [SNV-family](https://github.com/Clinical-Genomics-Lund/nextflow_wgs/blob/documentation/rank_models/rank_model_v5.01.ini)
	* [SV-single](https://github.com/Clinical-Genomics-Lund/nextflow_wgs/blob/documentation/rank_models/svrank_single_v5.1.ini)
	* [SV-family](https://github.com/Clinical-Genomics-Lund/nextflow_wgs/blob/documentation/rank_models/svrank_model_v5.1.ini)

