#!/usr/bin/perl -w
use strict;
use CMD::vcf qw( parse_vcf );
use List::Util qw( min sum );
use Data::Dumper;

my( $vcf, $run_dir, $assay ) = @ARGV;

my( $meta, $vars, $samples ) = parse_vcf($vcf);

my %p;
foreach my $var ( @$vars ) {
    foreach my $id ( keys %{$var->{GT}} ) {
	
	my @AD = split /,/, $var->{GT}->{$id}->{AD};
	my $DP = $var->{GT}->{$id}->{DP};
	next if !$DP or  $DP eq "." or $DP < 30 or $var->{CHROM} =~ /X|Y/;
	
	my $perc = min(@AD) / $DP;

	if( $perc < 0.30 ) {
	    push( @{$p{$id}}, $perc );
	}
    }
}

foreach my $id ( keys %p ) {
    system("/data/bnf/scripts/register_sample.pl --run-folder $run_dir --sample-id $id --assay $assay --contamination ". mean(@{$p{$id}}) );
    print $id."\t".mean(@{$p{$id}})."\n";
}



sub mean {
    if( @_ ) {
	return sum(@_)/@_;
    }
    else {
	return 0;
    }
}
