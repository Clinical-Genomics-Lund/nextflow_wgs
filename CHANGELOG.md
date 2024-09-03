# CHANGELOG

### TBD
* Extend the update_bed.pl script to handle multiple input files
* Rewrite to Python and add tests

### 3.9.8
* Reverted removed code in gene panel matches, caused missing gene panels for onco samples

### 3.9.7
* Solved trio eklipse image being wrongly added to yaml
* Removed outdated regex matches for genepanel, would remove important gene panels
* General clean-up of create_yml.pl

### 3.9.6
* Fix bug where wrong tuple value unpacked as group and sample id in `bqsr` when starting run from bam


### 3.9.5
* Fixed faulty if-condition for annotsv, would result in empty annotsv tsv everytime

### 3.9.4
* Use -K flag in bwa-mem for consistent results

### 3.9.3

* Re-optimized profiles wgs and onco. More memory allocations
* Added flag for reanalyze for bjorn to hook into

### 3.9.2

* Add updated and more communicative deploy script
* Remove or rename other deploy scripts

### 3.9.1

* Update MODY-cf configs to use the same as onco
* Clean up in MODY-cf config post merge

### 3.9.0
	
* Give the Sentieon container path by a parameter in the config file
* Update the Sentieon container to 202308 version
* Split out the `sentieon_qc` post-processing into its own process `sentieon_qc_postprocess`
* Update the Perl script used in `sentieon_qc_postprocess` to take input parameters as explicit arguments
* Update intersect file to latest used version of ClinVar (20231230)
* Update fastp to 0.23.4 and move to own container to fix reproducibility issue ([#143](https://github.com/Clinical-Genomics-Lund/nextflow_wgs/issues/143))
* Update CADD to v1.7
* Increase `inher_models` processing time
* Updated VEP from 103.0 to 111.0
* Updated VEP fasta from 98.0 to 111.0
* Updated VEP cache from 103.0 to 111.0
* Moved VEP parameters from processes to config
* Disabled vep `--everything` to disable VEP annotation w/ GNOMAD
* Removed deprecated `--af_esp` from `--everything`
* Tentative update of scout ranking.  
* cleanVCF.py now removes records missing CSQ-field.
* Add `SVTYPE` VEP 111 bug workaround in `vep_sv` process. (See  [Ensembl/ensembl-vep#1631](https://github.com/Ensembl/ensembl-vep/issues/1631#issuecomment-1985973568))
* Add VEP105 - 111 annotations to all rank models in use
* Fix onco model filename version (v5 rank model was misnamed as v4 in production)

### 3.8.2
* Re-enable D4 file generation (for Chanjo2)

### 3.8.1
* Disable Chanjo2

### 3.8.0
* Add mody-cf profile

### 3.7.14
* Run D4 coverage for full file

### 3.7.13
* Further simplifications of the checklist template

### 3.7.12
* Trim down size of checklist template, and add check for entering used test samples

### 3.7.11
* Add d4 file path directly to Scout YAML

### 3.7.10
* Tag Mitochondrial variants with GQ, loqusdb enabling

### 3.7.9
* Add CRON file to load Chanjo2

### 3.7.8
* added csv-file to onComplete function to accomodate CCCP

### 3.7.7
#### create_yml.pl
* small name change for myeloid constitutional to match clarity
* removed custom_images header for samples without images as pydantic would crash in scout load

### 3.7.6
* Add d4 coverage calculations to the workflow

### 3.6.6
* Fix genmod caller-penalty bug for GATK GQC vals ([#170](https://github.com/Clinical-Genomics-Lund/nextflow_wgs/issues/170))

### 3.6.5
* Remove bgzip and gunzip from versions
* Some cleanup in version documentation and code

### 3.6.4
* Use new docs as main entry point in repo
* Start removing old docs

### 3.6.3
* Update software responsible list in docs

### 3.6.2
* Added changelog reminder to github workflows

### 3.6.1
* Adding a new variant catalogue for expansionhunter/stranger/reviewer

### 3.5.11
* Add `documentation` to change type category in PR template.

### 3.6.0
* Changed melt configs, added flags: exome, removed flags: cov (was being used improperly)
* Added priors to `mei_list`, and changed `mei_list` to a new location in config
* Changes has been verified, report can be found internally

### 3.5.10
* Changed path to normal-pool-refs for gens. Uses masked hg38 references

### 3.5.9
* Add first iteration of updated documentation

### 3.5.8
* Move out resource files from `main.nf` to `nextflow.config`
* Move the selected fields for PHYLOP and PHASTCONS in vep to be specified in the process, similarly to the other plugins/custom fields

### 3.5.7
* Clean out unused files in repo root directory

### 3.5.6
* Add Github PR template/test documentation

### 3.5.5
* Update the cron log directory to use the `params.crondir` folder as base

### 3.5.4
* Add version outputs from all processes that use external software.
* Add stubs to processes to allow performing stub runs.

### 3.5.3
- Hotfix, increase melt sensivity by increasing amount of reads melt are alowed to use in RAM. 

### 3.5.2
- MELT is no longer filtered on location based upon regex names INTRONIC/null/PROMOTER, instead added a intersect towards bedfile. This will show splice site variants

### 3.5.1
* Add REVEL (Rare Exome Variant Ensemble Learner) Scores to VEP annotations (VEP `REVEL_rankscore` and `REVEL_score`)
	
### 3.5.0

#### Added

* Two processes for computing mitochondrial seq QC data from mt bam files and saving to JSON:
* Script `bin/merge_json_files.py` to merge 1 or more JSON files into one JSON.  Used to generate the final `{id}.QC` from the json output of the processes `sentieon_qc` and `sentieon_mitochondrial_qc`.
* Script `bin/mito_tsv_to_json.py` to extract and convert mtQC data from `sentieon_mitochondrial_qc` process output to json

#### Changed

* process `sentieon_qc` outputs to intermediate `{id}_qc.json` file instead of the final `{id}.QC`

### 3.4.3
#### patch genes
- added two more genes to expansionhunter variant catalogue.

### 3.4.2
#### hotfix
- dont print Mitochondrion, we handle the mitochondrion seperatly in the pipeline, caused loqusdb errors

### 3.4.1
#### minor additions/edits
- fixed filepaths for access-dir for myeloid profile in nextflow.config
- fixed assay name for create_yml.pl so yaml-file gets correct institute owner for myeloid samples

### 3.4.0
#### new features
- added a script to update wgs-bed file with current clinvar intron + intergenic regions. Also produces a log file of what's been added and removed
- added support to dry run vcf for testing scoring

### 3.3.3
#### minor improvements
- merged cnv2bed branch, small updates to color scheme for Alamut import files for CNVs
- increase time limit of create_pedigree

### 3.3.2
#### performance improvements
- added retries to vcfanno, file-system caching bug out?
- removed deep caching from freebayes(onco only) weird bug?

### 3.3.1
#### bug fixes
- alt affect type was lost for SNV<->SV compound, would get mixed up
  - added type and joined upon the value
- compounds for only SNVs for alt affect duos was wrongly renamed, added a sed-command

### 3.3.0
#### new features
- oncov2-0 and wgs profiles now both use loqusdb dumps for SV artefact annotations
- create pedigree has completely changed, now it is a separate perl-skript
  - creates one pedigree per affections status of parent, i.e in a trio, three ped-files with mother affect/father affected/no parent affected(default loaded into scout)
  - will calculate all states per genomod score, per vcf
  - optionally load these cases into scout, located in a subcategory in yaml-output-folder
- oncov2-0 now implemented, uses the old onco profile and oncov1-0 uses oncov1-0 profile (to be discontinued)
  - No longer use delly SV-caller, instead use GATK + CNVkit + manta
  - new version of MELT that catch much more important variation
  - indicator of onco-version in rankmodel-name of yaml-file, visable in scout case page

### 3.2.6
#### enhancement
- Added regex to support wgs-hg38-XXXX. suffix to run wgs-profile with different flags. ie --noupload true, no cdm/loqusdb upload for reruns

### 3.2.5
#### Bugfix
- Fixed a serious bug in prescore_sv.pl, would randomly chose proband-id for duos


### 3.2.4
#### Feature improvement
- Added SVs to loqusdb load. Using scored snv-vcf for correct MT->M notation

### 3.2.3
#### Bugfixes
- genmod patch not taking effect in singularity, switched to a smaller genmod container with patch
- updated processes in main.nf for above container

### 3.2.2
#### Bugfixes and improvement
- Changes to custom_images in yaml, case/str
- Added support for reviewer, activate when scout is updated

### 3.2.1
#### Bugfixes and improvement
- Image sizes for mitochondrial plots in yaml
- resource management in processes
- gatk-ref moved to cached directory

### 3.2.0
### Minor release
- Added support for loading images into scout, each process generating a plot can now be added as a path to scout-yaml
- Some support for Grace (new cluster)

### 3.1.2
#### Bugfixes
- CDM load-file only created for one individual of family, fixed (join function corrected)
- increased memory allocation for onco_depth
- removed shards from dedup, caused malformed output for dedup_metrics. Works as intended still

### 3.1.1
#### Bugfixes
- params.assay for onco has historical name, depth_onco process used wrong value
- using shard-specification for depup caused faulty dedupmetrics file

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
- GENS middleman command added to `generate_gens_data`. Needed for loading of data into GENS thorugh cron and middleman
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
  - modifications to `modify_vcf_scout.pl`, ignore maxentscan for M
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
- create\_yml.pl now recieved gene_panel content from hopper-json. no longer require scout-vm connectivity

### 2.1.7
- clincalwes now has correct loqusdb not piggybacking of onco

### 2.1.6
- timelimit increases, scratch and stage in/out for processes

### 2.1.5

- create\_yml.pl added ahus analysis for wgs_hg38 assay. Stinking mess initiated, please correct

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
- `filter_panel_cnv.pl` now only annotates for scout
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
