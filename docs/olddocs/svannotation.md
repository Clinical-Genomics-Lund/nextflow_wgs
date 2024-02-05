## Annotation and scout manipulation for structural variants

Several annotation steps are done after variant calling. Here all annotations of SNVs are documentated. If not specified, all SNV annotation is relevant for both mitochondrial and nuclear variants.


### AnnotSV
* input - merged SV VCF
  * mutation taster software
  * outputs tsv
    * score
    * gnomad frequencies
    * dbvar status (pathogenicity)
    * OMIM disease status (similar to later OMIM highest of these chosen)

### Variant effect predictor
* input - merged SV VCF
  * --merged
  * --everything
  * --synonyms $params.SYNONYMS
  * --vcf
  * --no_stats
  * --force_overwrite
  * --plugin LoFtool
  * --fasta $params.VEP_FASTA
  * --dir_cache $params.VEP_CACHE
  * --dir_plugins $params.VEP_CACHE/Plugins
  * --max_sv_size 50000000
  * --distance 200 -cache

### postprocess VEP
* input - VEP annotated SV VCF
  * cleanVCF.py, fixes csq-field of VEP
  * second svdb --merge 
    * 90% overlap to catch within sample overlaps (no --no_intra flag)
    * --notag (tags solved in first svdb merge)
  * add INFO header for VARID
  * add_omim.pl, adds OMIM disease annotation

  ### artefact
  * input - postprocessed VEP SV VCF, svdb database
      * add counts and frequency observations per variant


### prescore
* input - artefact SV VCF + annotsv tsv + pedigree file
  * prescore_sv.pl
    * adds annotsv annotations
    * adds genetic models (not used for scoring as of 3.0)
    * presents data for rankmodel

### score_sv
* input - prescored SV VCF +  rankmodel (single/family)
  * genmod score, single or trio
  * bcftools sort

### compound finder
* input - scored SV VCF and scored SNV VCF
  * finds SVs that are compounds to SNVs, adjusts scores of SNVs
  * only for trios
  * outputs adjusted SNV vcf