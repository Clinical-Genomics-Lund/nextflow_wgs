#CHANGELOG

###X.X.X

####Features
- Allow for non-distributed BWA (new default). For urgent cases use --shardbwa

####Fixes
- Change back to adding intersected vcf for loqusdb instead of full genomic vcf


###0.1.2

####Features

####Fixes
- Add 1000G in a special field for positions missing gnomAD (typically non-exonic) and add it to rankmodel

###0.1.1

####Features

####Fixes
- Properly add selected gene panels to YAML
- Add STR-vcf to YAML
- Rename sample in STR-vcf to agree with sample name instead of bam filename
