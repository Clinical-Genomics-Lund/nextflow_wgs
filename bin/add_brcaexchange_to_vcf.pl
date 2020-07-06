#!/usr/bin/perl -w
use strict;
use File::Basename;
use lib dirname (__FILE__);
use tsv qw(read_tsv);
use Data::Dumper;


my @data = read_tsv($ARGV[1]);

my %enigma_clinsig;
foreach my $var ( @data ) {
    my( $chr, $pos, $ref, $alt ) = ( $var->{'Genomic_Coordinate_hg38'} =~ /chr(.+?):g\.(\d+?):(.+?)>(.+?)$/ );
    $enigma_clinsig{$chr."_".$pos."_".$ref."_".$alt} = $var->{'Clinical_significance_ENIGMA'};
}

open( VCF, $ARGV[0] );

while( <VCF> ) {

    unless( /^#/ ) {
	my @var = split /\t/;
	my $short = $var[0]."_".$var[1]."_".$var[3]."_".$var[4];
	if( $enigma_clinsig{$short} ) {
	    $var[7] .= ";ENIGMA_CLNSIG=".$enigma_clinsig{$short};
		$var[7] .= ";SCOUT_CUSTOM="."ENIGMA|".$enigma_clinsig{$short};
	}
	print join( "\t", @var );
    }
    else {
	print ;
    }
}
