# CHANGELOG

### 0.1.4

#### Features

#### Fixes
- Increase allocated  memory for Clinvar SnpSift process due to occasional crashes
- Further fixes to output folders
- Retry create_yaml and loqus process up to 5 times

### 0.1.3

#### Features
- Allow for non-distributed BWA (new default). For urgent cases use --shardbwa
- Add diagnosis field to CDM

#### Fixes
- Change back to adding intersected vcf for loqusdb instead of full genomic vcf
- Change name of 1000G INFO field to make it show up in Scout
- Removed some hardcoded paths to scripts and added them the the bin/ folder

### 0.1.2

#### Features

#### Fixes
- Add 1000G in a special field for positions missing gnomAD (typically non-exonic) and add it to rankmodel

### 0.1.1

#### Features

#### Fixes
- Properly add selected gene panels to YAML
- Add STR-vcf to YAML
- Rename sample in STR-vcf to agree with sample name instead of bam filename
