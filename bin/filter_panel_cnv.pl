#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename;
use lib dirname (__FILE__);
use vcf2;
use strict;
use Data::Dumper;

my %opt = ();
GetOptions( \%opt, 'manta=s', 'delly=s'  );



my $vcf = CMD::vcf2->new('file'=>$ARGV[0] );

my $bed = $ARGV[1];
open (BED, $bed) or die $!;
my %bed;
my $probe = 1;
while (<BED>) {
    my @tmp = split/\t/;

    $bed{$tmp[0]}->{$probe}->{START} = $tmp[1];
    $bed{$tmp[0]}->{$probe}->{END} = $tmp[2]; 
    $probe++;
}
my $ref = \%bed;

my @header = split/\n/,$vcf->{header_str};
my $count = 1;
foreach my $header (@header) {
    $header =~ s/\'\\\|\'/pipe-sign/g;
    
    if ($header =~ /^##INFO/ && $count == 1) {
        print '##INFO=<ID=SCOUT_CUSTOM,Number=.,Type=String,Description="Custom annotations for scout">'."\n";
        print '##INFO=<ID=MELT_RANK,Number=.,Type=Number,Description="Evidence level 1-5, 5highest">'."\n";
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
    my $delly = 0;
    my $manta = 0;
    my $cnvkit = 0;
    my @callers = split/-/,$a->{INFO}->{set};
    foreach my $caller (@callers) {
        if ($caller =~ /delly/) {
            $delly = 1;
        }
        elsif ($caller =~ /manta/) {
            $manta = 1;
        }
        elsif ($caller =~ /cnvkit/) {
            $cnvkit = 1;
        }
        elsif ($caller =~ /Intersection/) {
            $delly = 1; $manta = 1; $cnvkit = 1;
        }
    }
    next if ($a->{INFO}->{SVTYPE} eq 'BND');
    ## Filter delly-only variants
    my $check = 1;
    if ($delly == 1 && $manta == 0 && $cnvkit == 0) {
        $check = delly($a);
    }
    
    next if ($check == 0);
    #print join('_',@callers)."\n";
    
    ## Print first 7 columns of vcf
    print join("\t",@str[0..6])."\t";
    
    ## ADD END to info and change SVTYPE to INS
    my @INFO = split/;/,$str[7];
    #print $a->{vcf_str}."\n";

    my @foundin;
    if ($manta) { push @foundin,"manta"; }
    if ($delly) { push @foundin,"delly"; }
    if ($cnvkit) { push @foundin,"cnvkit"; }
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


sub delly {
    my $a = shift;
    my $start = $a->{POS};
    my $end = $a->{INFO}->{END};
    my $len = $end - $start + 1;
    my $check = 1;
    ## remove all IMPRECISE called variants smaller than 1000bp
    if ($len <= 1000 && $a->{INFO}->{IMPRECISE}) {
        $check = 0;
    }
    ## remove inversions called with IMPRECISE method
    elsif ($a->{INFO}->{SVTYPE} eq 'INV' && $a->{INFO}->{IMPRECISE}) {
        $check = 0;
    }
    ## remove deletions with copy-number 2
    #elsif ($a->{INFO}->{SVTYPE} eq 'DEL' && $a->{GT}->[0]->{CN} == 2) {
    #    $check = 0;
    #}
    ## remove duplications with copy-number 2
    # elsif ($a->{INFO}->{SVTYPE} eq 'DUP' && $a->{GT}->[0]->{CN} == 2) {
    #     $check = 0;
    # }
    ## High-quality variant junction reads
    elsif ($a->{INFO}->{PRECISE} && $a->{GT}->[0]->{RV} <= 5) {
        $check = 0;
    }
    ## remove low-quality IMPRECISE variants   
    elsif ($a->{FILTER} ne 'PASS' && $a->{INFO}->{IMPRECISE}) {
        $check = 0;
    }
    
    #print Dumper($a);
    ## Remove variants where no breakpoints of precise variants are within design of panel
    else {
        my $s_ok = 0;
        my $e_ok = 0;
        foreach my $probe (keys %{ $ref->{ $a->{CHROM} } }) {
            my $start_bed = $ref->{ $a->{CHROM} }->{$probe}->{START};
            my $end_bed = $ref->{ $a->{CHROM} }->{$probe}->{END};
            #print "$a->{CHROM}:  $start - $end - $start_bed - $end_bed\n";

            if ($start >= $start_bed-150 && $start <= $end_bed+150) {
                $s_ok = 1;
            }
            if ($end <= $start_bed-150 && $end >= $end_bed+150) {
                $e_ok = 1;
            }    
        }
        if ($s_ok + $e_ok < 1) {
            $check = 0;
        }
    }
    ## remove IMPRECISE variants on non-covered chromosomes
    if (!defined $ref->{ $a->{CHROM} }) {
        $check = 0;
    }
    #if ($check == 1) {   print $len."____"; }
    return $check;
}