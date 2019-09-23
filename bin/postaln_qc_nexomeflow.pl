#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use IPC::Cmd qw[can_run run];
use JSON;# qw( encode_json );

my $hsmetrics = $ARGV[0];
my $reads = $ARGV[1];
my $ins = $ARGV[2];
my $depth = $ARGV[3];
my $SID = $ARGV[4];
my %results;
## HSMETRICS ##
open( HS, $hsmetrics );
while( <HS> ) {
    if( /^\#\# METRICS CLASS/ ) {
	    <HS>;
	    my $vals = <HS>;
	    my @a = split /\t/, $vals;
        $results{'pct_on_target'} = $a[18];
	    $results{'fold_enrichment'} = $a[26];
    	$results{'median_coverage'} = $a[23];
    	$results{'fold_80'} = $a[33];
	}
}
close HS;


## READS ##
my @flagstat;
open( READS, $reads );
while( <READS> ) {
    push @flagstat, $_;
}
my( $num_reads ) = ( $flagstat[0] =~ /^(\d+)/ );
my( $dup_reads ) = ( $flagstat[3] =~ /^(\d+)/ );
my( $mapped_reads ) = ( $flagstat[4] =~ /^(\d+)/ );
close READS;

## INSSIZE ##
open( INS, $ins );
while( <INS> ) {
	if( /^\#\# METRICS CLASS/ ) {
	    <INS>;
	    my $vals = <INS>;
	    my @a = split /\t/, $vals;
	    $results{'ins_size'} = $a[0];
	    $results{'ins_size_dev'} = $a[1];
	}
}
close INS;

my @thresholds = qw( 1 10 30 100 250 500 1000);
## DEPTH ##

my( $pct_above, $mean_cov ) = parse_basecov_bed( $depth, \@thresholds );

$results{pct_above_x} = $pct_above;
$results{tot_reads} = $num_reads;
$results{mapped_reads} = $mapped_reads;
$results{dup_reads} = $dup_reads;
$results{dup_pct} = $dup_reads / $mapped_reads;
$results{sample_id} = $SID;
$results{mean_cov} = $mean_cov;


my $json = JSON->new->allow_nonref;
print $json->pretty->encode( \%results );



sub parse_basecov_bed {
    my( $fn, $thresholds ) = @_;
    open( my $cov_fh, $fn );

    chomp( my $head_str = <$cov_fh> );
    
    $head_str =~ s/^#\s+//;
    my @head = split /\t/, $head_str;
    my $cov_field;
    for my $i ( 0..$#head ) {
	$cov_field = $i if $head[$i] eq "COV";
    }
 

    my $tot_bases = 0;
    my %above_cnt;
    my( $tot, $cnt ) = (0,0);
    while( <$cov_fh> ) {
	chomp;
	my @a = split /\t/;

	next if $a[0] =~ /^chr(Un|\d+_)/ ;

	$tot += $a[2];
	$cnt++;
	
	$tot_bases ++;
	foreach my $min ( @$thresholds ) {
	    $above_cnt{ $min } ++ if $a[ $cov_field ] >= $min; 
	}
    }

    my %above_pct;
    #print "Total: $tot_bases\n";
    foreach( sort {$a<=>$b} keys %above_cnt ) {
	$above_pct{ $_ } = 100 * ($above_cnt{$_} / $tot_bases);
    }

    my $mean_cov = $tot / $cnt;
    
    return \%above_pct, $mean_cov;
}



sub parse_cov_bed {
    my( $fn, $thresholds ) = @_;
    open( my $cov_fh, $fn );

    chomp( my $head_str = <$cov_fh> );
    
    $head_str =~ s/^#\s+//;
    my @head = split /\t/, $head_str;
    my $cov_field;
    for my $i ( 0..$#head ) {
	$cov_field = $i if $head[$i] eq "meanCoverage";
    }
 

    my $tot_exons = 0;
    my %above_cnt;
    while( <$cov_fh> ) {
	chomp;
	my @a = split /\t/;
	$tot_exons ++;
	foreach my $min ( @$thresholds ) {
	    $above_cnt{ $min } ++ if $a[ $cov_field ] >= $min; 
	}
    }

    foreach( sort {$a<=>$b} keys %above_cnt ) {
	print "$_\t";
	printf "%.2f%%", 100 *($above_cnt{$_} / $tot_exons);
	print "\n";
    }

}

