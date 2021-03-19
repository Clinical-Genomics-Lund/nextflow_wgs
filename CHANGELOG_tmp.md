mito_and_containerupdates

new functions
## main.nf ##
- added mito-calling
  - mutect2
  - hmtnote
  - haplogrep
  - eklipse
- added SMNCopyNumberCalling
- New VEP container and version (103)
- gatk cnv calling
  - adjustments to all affected scripts
## container ##
- new container specifically for madeline2
- main container now includes all software except madeline2 and VEP
  - new conda environments
  - updates to Expansionhunter
  - updates to Stranger
  - updated GATK version
  - updated Sentieon
  - added: haplogrep, hmtnote, eklipse, melt, graphalignmentviewer, SMNcopynumbercaller, CNVkit and imagemagick
