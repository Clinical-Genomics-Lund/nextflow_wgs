#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename;
use lib dirname (__FILE__);
use vcf2;
use strict;
use Data::Dumper;

my %opt = ();
GetOptions( \%opt, 'mergedvcf=s', 'callers=s'  );



my $vcf = CMD::vcf2->new('file'=>$opt{mergedvcf} );

my $used_callers = $opt{callers};

my @header = split/\n/,$vcf->{header_str};
my $count = 1;
foreach my $header (@header) {
    $header =~ s/\'\\\|\'/pipe-sign/g;
    
    if ($header =~ /^##INFO/ && $count == 1) {
        print '##INFO=<ID=SCOUT_CUSTOM,Number=.,Type=String,Description="Custom annotations for scout">'."\n";
        print '##INFO=<ID=MELT_RANK,Number=.,Type=String,Description="Evidence level 1-5, 5highest">'."\n";
        print '##INFO=<ID=MELT_QC,Number=.,Type=String,Description="Quality of call">'."\n";
        $count++;
    }
    else {
        print $header."\n";
    }

}
while ( my $a = $vcf->next_var() ) {

    my @str = split/\t/,$a->{vcf_str};
    my $start = $a -> {POS};

    #print Dumper($a);
    my $manta = 0;
    my $cnvkit = 0;
    my $gatk = 0;
    my @callers = split/-/,$a->{INFO}->{set};
    foreach my $caller (@callers) {

        if ($caller =~ /manta/) {
            $manta = 1;
        }
        elsif ($caller =~ /cnvkit/) {
            $cnvkit = 1;
        }
        elsif ($caller =~ /gatk/) {
            $gatk = 1;
        }
        elsif ($caller =~ /Intersection/) {
            if ($used_callers =~ /manta/) {
                $manta = 1;
            }
            if ($used_callers =~ /gatk/) {
                $gatk = 1;
            }
            if ($used_callers =~ /cnvkit/) {
                $cnvkit = 1;
            } 
        }
    }

    next if ($a->{INFO}->{SVTYPE} eq 'BND');

    ## Print first 7 columns of vcf
    print join("\t",@str[0..6])."\t";
    
    ## ADD END to info and change SVTYPE to INS
    my @INFO = split/;/,$str[7];
    #print $a->{vcf_str}."\n";

    my @foundin;
    if ($manta) { push @foundin,"manta"; }
    if ($cnvkit) { push @foundin,"cnvkit"; }
    if ($gatk) { push @foundin,"gatk"; }
    #print Dumper($a);
    push @INFO,"SCOUT_CUSTOM=Caller|".join('&',@foundin);
    print join(';',@INFO)."\t";
    #print join("\t",@str[8..$#str]);
    print $str[8]."\t";
    ## 0/0 genotypes cannot be loaded as artefacts, if 0/0 set as 0/1
    my $gtmod = $str[9];
    $gtmod =~ s/^0\/0/0\/1/;
    print $gtmod;
    print "\n";

}
