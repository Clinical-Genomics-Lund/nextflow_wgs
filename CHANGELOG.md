# CHANGELOG

### 3.1.0
#### Performance Improvements
- Removed all distribution of sentieon except first alignment step option
- Removed bam.toRealPath() from all processes. Bai files are now given along with bam files if alignment is to be skipped. More down below
#### New/Updated Features
- VCF start removed temporarily
- BAM start now work better. Add headers bam + bai to csv with corresponding files
- BATCH start now available for onco-samples. Thorough channel joining and removal of distributed sentieon made it possible (does not work for wgs profile!)

### 3.0.9
#### bug fixes
- fixed grep for multi-allelic FP loci

### 3.0.8
#### feature
- added hash element rstreshold (rankscore threshold), that if defined overwrites defualt -1 to createyml.pl

### 3.0.7
#### feature
- added panel depth as alternative to chanjo for panel-data

### 3.0.6
#### bug fixes
- correctly assigned theanoflag for gatk coverage and ploidy, would in rare cases cause crashes

### 3.0.5
#### new funcion
- GENS middleman command added to generate_gens_data. Needed for loading of data into GENS thorugh cron and middleman
#### bug fixees
- REViewer now loops through a perl shell script instead of bash. Low covered loci error no longer crash all other svg-image generation
- fixed a typo which named all svgs as 7156, a validation and verification sample

### 3.0.4
- recurring multi-allelic variant @MT:955 keeps vcfmultibreak in a never ending loop
- grep -v ^MT 955

### 3.0.3
- ignore errors of REViewer
- loqusdb faulty input caused wrong imports to loqusdb, now fixed

### 3.0.2
- GAV replaced with REViewer
- Stranger 0.8 with updated variant catalogue

### 3.0.1
- source activate gatk to all gatk processes that use 4.1.9 and set +eu. unbound errors
- increased memory allocation for several gatk and mito processes
- sharded merged bam and non-sharded bam now produces output for dedup too. Locuscollector no longer redirects bam. This saves upto 70% of temporary files!

### 3.0.0
#### Summary of changes
new functions
### main.nf
- added mito-calling
  - mutect2
  - hmtnote
  - haplogrep
  - eklipse
  - modifications to filter_indels (new VEP fields)
  - modifications to modify_vcf_scout.pl, ignore maxentscan for M
- added SMNCopyNumberCalling
- New VEP container and version (103)
- gatk cnv calling
  - adjustments to all affected scripts
### container
- new container specifically for madeline2
- main container now includes all software except madeline2 and VEP
  - new conda environments
  - updates to Expansionhunter
  - updates to Stranger
  - updated GATK version
  - updated Sentieon
  - added: haplogrep, hmtnote, eklipse, melt, graphalignmentviewer, SMNcopynumbercaller, CNVkit and imagemagick

### minor improvements
- group and sample IDs of outputs re-thought
- contig synonyms for VEP
- BAM start working better


### 2.1.12
#### params.panelsdef
- added path to symlinked latest weekly definition.

### 2.1.11
#### qc-values added
- added pf missmatch and error rates to qc-json

### 2.1.10
#### scout presentation
- rankmodels now separate VEP-consequence from AnnotSVrank and dbvar in Consequence and Clinical_significance respectively


### 2.1.9
#### bugfix
- rescore.nf had wrongly named variable in output for bamfiles 

### 2.1.8
- create_yml.pl now recieved gene_panel content from hopper-json. no longer require scout-vm connectivity

### 2.1.7
- clincalwes now has correct loqusdb not piggybacking of onco

### 2.1.6
- timelimit increases, scratch and stage in/out for processes

### 2.1.5

- create_yml.pl added ahus analysis for wgs_hg38 assay. Stinking mess initiated, please correct

### 2.1.4

- create_yml.pl added hemato analysis for clinicalwesv1-0 assay. Correct institute for myeloid normals

### 2.1.3

#### Features
- create_yml.pl now has a hash with all definable scout import-fields per assay, allowing easier additions and modifications to/of assays.

### 2.1.2

#### Fixes
- Fix a bug that generated corrupt Gens json files...

### 2.1.1

#### Features
- Generate a json with data for the Gens overview plot, to allow quicker loading in Gens.


### 2.1.0

#### Features
- added rescoring function through rescore.nf
#### misc
- fixed naming of expansionhunter vcfs
- new rankmodels for wgs profile. 5.1 (loqusdb cutoffs and VEP-csq scoring)


### 2.0.2

#### Features
- added specific delly filtering script
  - now correctly filters breakpoints outside panel
- filter_panel_cnv.pl now only annotates for scout
  - delly precise/imprecise annotation
- new artifact database for SV-calling for WGS and oncogenetics

### 2.0.1

#### Features
- added container and git-hash to logging

#### Fixes
- pathing through freebayes cause nextflow to not recieve completion status, now wgs is run through freebayes with touch-command only

### 2.0

#### Features
- Optional input files, fastq/bam/vcf
- hg38 alignment and annotations
- profiles, wgs/onco/exome
- several new sv-variant callers, melt, cnvkit, delly
- POD-tool for duplication events in trios
- Freebayes calling for difficult homopylomers in onco-samples
- Yaml-creation for scout import overhauled
- new container with needed software
- new rank-models for onco (both snv and sv) and wgs (sv-rank not live yet)
- Various small improvements of code

#### Fixes
- Optimization of cpu/memory/time for each process
- Numerous small improvements of several scripts


### 1.5.4

#### Features
- Last hg19 version

#### Fixes
- Minor file-paths issues resolved

### 1.5.3

#### Fixes
- Fix incorrect filtering of UPD calls in genomeplotter

### 1.5.2

#### Fixes
- Only plot UPDs in overview plot if > 100 informative sites
- Don't run UPD process on duos

### 1.5

#### Features

#### Fixes
- Use the correct scout server for create_yaml and loqusdb annotation
- Removed some hardcoded assumptions in create_yml.pl
- Always create gvcfs, and publish
- Publish chanjo coverage file
- Fix ID mixup in gatkcov process
- Added retry strategies to cdm, locuscollector and bqsr processes
- create_yml: Only add each panel once and sort panels alphabetically

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
