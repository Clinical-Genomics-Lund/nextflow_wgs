#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use Getopt::Long;
use List::MoreUtils qw(first_index);
use File::Basename;
use Backticks;
use Text::Unidecode;

my $container = $ARGV[0];

my $exec = "singularity exec --bind /data ".$container;
my $cmd;
my $status;
my $results;
my $count_tot;
my $count_success;
### fastp ### fastp
$count_tot++;
print "testing fastp...\n";
$cmd = $exec." fastp";
$results = `$cmd`;
$status = $results->success;
if ($status) { print "OK\n"; $count_success++;}
### sentieon ###
$count_tot++;
print "testing sentieon...\n";
$cmd = $exec." sentieon driver";
$results = `$cmd`;
$status = $results->success;
if ($status) { print "OK\n";$count_success++;}
### sambamba ### sambamba depth
$count_tot++;
print "testing sambamba...\n";
$cmd = $exec." sambamba depth";
$results = `$cmd`;
$status = $results->success;
if ($status) { print "OK\n";$count_success++;}
### ExpansionHunter ### ExpansionHunter
$count_tot++;
print "testing ExpansionHunter...\n";
$cmd = $exec." ExpansionHunter";
$results = `$cmd`;
$status = $results->success;
if ($status) { print "OK\n";$count_success++;}
### stranger ### source activate py3-env stranger
$count_tot++;
print "testing stranger...\n";
$cmd = $exec." bash -c \"source activate py3-env; stranger --help\"";
$results = `$cmd`;
$status = $results->success;
if ($status) { print "OK\n";$count_success++;}
### break multi ### vcfbreakmulti
$count_tot++;
print "testing vcfbreakmulti...\n";
$cmd = $exec." vcfbreakmulti --help";
$results = `$cmd`;
$status = $results->success;
if ($status) { print "OK\n";$count_success++;}
### picard ###
$count_tot++;
print "testing picard...\n";
$cmd = $exec." java -jar /opt/conda/envs/CMD-WGS/share/picard-2.21.2-1/picard.jar CollectHsMetrics --version 2> /dev/null";
$results = system($cmd);
if ($results =~ /256/) { print "OK\n";$count_success++;}
### melt ### java -jar  /opt/MELT.jar source activate java8-env
$count_tot++;
print "testing melt...\n";
$cmd = $exec." bash -c \"source activate java8-env; java -jar /opt/MELT.jar  Single -help \" > /dev/null";
$results = system($cmd);
if ($results =~ /256/) { print "OK\n";$count_success++;}
### ped_parser ###
$count_tot++;
print "testing ped_parser...\n";
$cmd = $exec." ped_parser --help";
$results = `$cmd`;
$status = $results->success;
if ($status) { print "OK\n";$count_success++;}
### madeline2 ###
$count_tot++;
print "testing madeline2...\n";
$cmd = $exec." madeline2 --help";
$results = `$cmd`;
$status = $results->success;
if ($status) { print "OK\n";$count_success++;}
### freebayes ###
$count_tot++;
print "testing freebayes...\n";
$cmd = $exec." freebayes -help";
$results = `$cmd`;
$status = $results->success;
if ($status) { print "OK\n";$count_success++;}
### bcftools ### 
$count_tot++;
print "testing bcftools...\n";
$cmd = $exec." bash -c \"bcftools roh \" 2> /dev/null";
$results = system($cmd);
if ($results =~ /256/) { print "OK\n";$count_success++;}
### bedtools ### bedtools
$count_tot++;
print "testing bedtools...\n";
$cmd = $exec." bedtools ";
$results = `$cmd`;
$status = $results->success;
if ($status) { print "OK\n";$count_success++;}
### vcfanno_linux64 ###
$count_tot++;
print "testing vcfanno...\n";
$cmd = $exec." vcfanno_linux64 ";
$results = `$cmd`;
$status = $results->success;
if ($status) { print "OK\n";$count_success++;}
### genmod ### genmod
$count_tot++;
print "testing genmod...\n";
$cmd = $exec." genmod ";
$results = `$cmd`;
$status = $results->success;
if ($status) { print "OK\n";$count_success++;}
### gunzip ### gunzip/bgzip
$count_tot++;
print "testing bgzip...\n";
$cmd = $exec." bgzip ";
$results = `$cmd`;
$status = $results->success;
if ($status) { print "OK\n";$count_success++;}
### tabix ###
$count_tot++;
print "testing tabix...\n";
$cmd = $exec." bash -c \"tabix \" 2> /dev/null";
$results = system($cmd);
if ($results =~ /256/) { print "OK\n";$count_success++;}
### peddy ### source activate py3-env python -m peddy
$count_tot++;
print "testing peddy...\n";
$cmd = $exec." bash -c \"source activate py3-env; python -m peddy -h\"";
$results = `$cmd`;
$status = $results->success;
if ($status) { print "OK\n";$count_success++;}
### upd ### upd
$count_tot++;
print "testing upd...\n";
$cmd = $exec." upd";
$results = `$cmd`;
$status = $results->success;
if ($status) { print "OK\n";$count_success++;}
### gatk ### source activate gatk4-env gatk
$count_tot++;
print "testing gatk...\n";
$cmd = $exec." bash -c \"source activate gatk4-env; gatk\"";
$results = `$cmd`;
$status = $results->success;
if ($status) { print "OK\n";$count_success++;}
### manta ### configManta.py
$count_tot++;
print "testing manta...\n";
$cmd = $exec." bash -c \" configManta.py --version\"";
$results = `$cmd`;
$status = $results->success;
if ($status) { print "OK\n";$count_success++;}
### delly ### delly
$count_tot++;
print "testing delly...\n";
$cmd = $exec." delly";
$results = `$cmd`;
$status = $results->success;
if ($status) { print "OK\n";$count_success++;}
### cnvkit ### cnvkit.py
$count_tot++;
print "testing cnvkit...\n";
$cmd = $exec." cnvkit.py -h";
$results = `$cmd`;
$status = $results->success;
if ($status) { print "OK\n";$count_success++;}
### svdb ### source activate py3-env svdb 
$count_tot++;
print "testing svdb...\n";
$cmd = $exec." bash -c \"source activate py3-env svdb -h\" 2> /dev/null";
$results = system($cmd);
if ($results =~ /256/) { print "OK\n";$count_success++;}
### vcf-concat ### vcf-concat
$count_tot++;
print "testing vcf-concat...\n";
$cmd = $exec." bash -c \"vcf-concat -h\" 2> /dev/null";
$results = system($cmd);
if ($results =~ /65280/) { print "OK\n";$count_success++;}
### TIDDIT ### TIDDIT.py
$count_tot++;
print "testing TIDDIT...\n";
$cmd = $exec." TIDDIT.py -h";
$results = `$cmd`;
$status = $results->success;
if ($status) { print "OK\n";$count_success++;}
### mongodb ### create_yml.pl
$count_tot++;
print "testing mongodb...\n";
$cmd = $exec." ./perlmoduletester.pl";
$results = system($cmd);
$results = `$cmd`;
$status = $results->success;
if ($status) { print "OK\n";$count_success++;}

print $count_success." passed out of ".$count_tot."\n";