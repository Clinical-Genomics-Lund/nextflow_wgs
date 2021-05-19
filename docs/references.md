## This is the reference docs

* Reference genome, hg38, also subdir with split per chromosome

	`wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz`

	`sed 's/^>chr/>/' GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz > GCA_000001405.15_GRCh38_no_alt_analysis_set_nochr.fna.gz`

	`bwa index GCA_000001405.15_GRCh38_no_alt_analysis_set_nochr.fna.gz`

	`java -jar picard.jar CreateSequenceDictionary REFERENCE=GCA_000001405.15_GRCh38_no_alt_analysis_set_nochr.fna OUTPUT=GCA_000001405.15_GRCh38_no_alt_analysis_set_nochr.dict`

	`samtools faidx GCA_000001405.15_GRCh38_no_alt_analysis_set_nochr.fna.gz`

* genome shards (in the repo)
* [SNV-CADD](https://kircherlab.bihealth.org/download/CADD/v1.6/GRCh38/whole_genome_SNVs.tsv.gz)
* [VEP-cache](https://m.ensembl.org/info/docs/tools/vep/script/vep_cache.html) Match your vep version.
* gnomAD-exomes (vcf)
* gnomAD-genomes (vcf)
* [PhyloP](http://hgdownload.cse.ucsc.edu/goldenpath/hg38/)
* [PhastCons](http://hgdownload.cse.ucsc.edu/goldenpath/hg38/)
* [Mills_and_1000G_gold_standard.indels](http://hgdownload.cse.ucsc.edu/goldenpath/hg38/)
* ClinVar-vcf, keep this updated weekly - ADD TO DOC!
* Intersect-bed to reduce the amount of SNVs (basicly exome + interesting regions)
* GATK CNV bins (interval list) - ADD TO DOC!
* GATK CNV PON, males (hdf5) - ADD TO DOC!
* GATK CNV PON, females (hdf5) - ADD TO DOC!
* Fastgnomad SNP-list (dat) - ADD TO DOC!
* Gens SNP-list (txt) - ADD TO DOC!
* SVDB-databas (db) - might be going away, (https://github.com/J35P312/SVDB) SVs
* [loqusdb export](https://github.com/moonso/loqusdb) (vcf) - Local observations, both artefacts and real variants. SNVs
* Ensembl genes (bed)
* Gene panel-dump (json) - Needed for last step to create yaml for scout.
* Rankmodels
	* [SNV-single](https://github.com/Clinical-Genomics-Lund/nextflow_wgs/blob/documentation/rank_models/rank_model_v5.01_single.ini)
	* [SNV-family](https://github.com/Clinical-Genomics-Lund/nextflow_wgs/blob/documentation/rank_models/rank_model_v5.01.ini)
	* [SV-single](https://github.com/Clinical-Genomics-Lund/nextflow_wgs/blob/documentation/rank_models/svrank_single_v5.1.ini)
	* [SV-family](https://github.com/Clinical-Genomics-Lund/nextflow_wgs/blob/documentation/rank_models/svrank_model_v5.1.ini)

