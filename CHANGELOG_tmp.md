main.nf updated to hg38
config updated to hg38
rank_models updated to hg38 spidex/gerp removed

cosmetics
 - removed most variables from main.nf, called through $params.variable from config file
 - renamed several input channels to more accurately describe where their destination is
 - small indendation fixes

new functions
 - MAJOR, added choices to input-files. Now input can be bam and vcf in addition to fastq. read1 and read2 in input-csv can now be bam + bqsr table or vcf + vcf.idx respectively.
 - flags --align --varcall --annotate can be combined with above inputs to further modify what results is desired.
 - added --varcall --gatkcov flags. Defaulted to false.


