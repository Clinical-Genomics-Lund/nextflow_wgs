#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename;
use lib dirname (__FILE__);
use vcf2;
use strict;
use Data::Dumper;


my $line = 1;
my $index;
my $vcf = CMD::vcf2->new('file'=>$ARGV[0] );
my $proband = $ARGV[1];
print $vcf->{header_str};
while ( my $a = $vcf->next_var() ) {

    
    if ($line == 1) {
        my $count = 0;
        foreach my $ind (@{ $a->{GT} }) {
            if ($ind->{_sample_id} eq $proband) {
                $index = $count;
            }
            $count++;
        }
    }

    ## variant exist in proband
    if ($a->{GT}->[$index]->{GT} =~ /1/) {
        ## variant has at least 50 depth    
        if ($a->{GT}->[$index]->{DP} > 50) {
            print $a->{vcf_str}."\n";
        }
    }
    $line++;
}