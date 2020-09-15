#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename;
use lib dirname (__FILE__);
use vcf2;
use strict;
use Data::Dumper;

my %opt = ();
GetOptions( \%opt, 'vcf=s', 'bed=s', 'tumor-id=s' );

my $vcf = CMD::vcf2->new('file'=>$opt{vcf} );

my $bed = $opt{bed};
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
    #$header =~ s/\'\\\|\'/pipe-sign/g;
    
    if ($header =~ /^##INFO/ && $count == 1) {
        $count++;
    }
    else {
        print $header."\n";
    }

}
while ( my $a = $vcf->next_var() ) {

    my @str = split/\t/,$a->{vcf_str};
    my $start = $a -> {POS};

    ## ignore BNDs
    next if ($a->{INFO}->{SVTYPE} eq 'BND');
    ## Filter delly-only variants
    my $check = 1;
    my $paired;
    if (scalar(@{$a->{GT}}) > 1 ) {
        $paired = "true";
    }
    # Filter function
    $check = delly($a);
    
    next if ($check == 0);
    
    ## Print first 7 columns of vcf
    print join("\t",@str[0..6])."\t";
    
    ## ADD END to info and change SVTYPE to INS
    my @INFO = split/;/,$str[7];

    print join(';',@INFO)."\t";
    #print join("\t",@str[8..$#str]);
    print $str[8]."\t";
    ## if variant survives filters and still is 0/0 
    # those genotypes cannot be loaded as artefacts, as such 0/0 is set to 0/1
    my $gtmod = $str[9];
    if ($a->{INFO}->{IMPRECISE}) {
        $gtmod =~ s/^0\/0/0\/1/;
    }
    print $gtmod;
    if ($paired) {
        my $gtmod2 = $str[10];
        if ($a->{INFO}->{IMPRECISE}) {
            $gtmod2 =~ s/^0\/0/0\/1/;
        }
        print "\t".$gtmod2;
    }
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
    ## remove precise called variants smaller than 300bp
    elsif ($len <= 300) {
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
    
    ## High-quality variant junction reads, AND if precise require actual genotype call!
    elsif ($a->{INFO}->{PRECISE} ) {
        my $RVOK = 0;
        my $GENOK = 0;
        foreach my $ind ( @{$a->{GT}} ) {
            if ($ind->{RV} >= 5) {
                $RVOK = 1;
            }
            if ($ind->{GT} eq '0/1' || $ind->{GT} eq '1/1' ) {
                $GENOK = 1;
            }
        }
        if ($RVOK == 0 || $GENOK == 0) {
            $check = 0;
        }
    }
    ## remove low-quality IMPRECISE variants   
    elsif ($a->{FILTER} ne 'PASS' && $a->{INFO}->{IMPRECISE}) {
        $check = 0;
    }
    ## High-quality variant reads, AND if imprecise require actual genotype call!
    elsif ($a->{INFO}->{IMPRECISE} ) {
        my $DVOK = 0;
        my $GENOK = 0;
        foreach my $ind ( @{$a->{GT}} ) {
            if ($ind->{DV} >= 20) {
                $DVOK = 1;
            }
            if ($ind->{GT} eq '0/1' || $ind->{GT} eq '1/1' ) {
                $GENOK = 1;
            }
        }
        if ($DVOK == 0 || $GENOK == 0) {
            $check = 0;
        }
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