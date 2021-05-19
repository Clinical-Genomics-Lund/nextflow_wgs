## Create YAML for scout loading

Various files and data needs to be loaded per sample into scout. Information is gathered throughout pipeline execution. It's gathered by create_yml.pl


### create_yml.pl
* input - group.INFO file, genelist.json, pedigree-file and various meta-info
  * files created in pipeline gets stored into ${group}.INFO
    * BAM per sample
    * SV_VCF
    * SNV_VCF
    * STR_VCF
    * madeline2 xml
    * peddy outputs
  * scout json dump of gene-panels is read and added to YAML
  * What institute should sample belong to --assay
  * create_yml.pl collects all information, has hash-defined at start to separate sample into correct scout institute