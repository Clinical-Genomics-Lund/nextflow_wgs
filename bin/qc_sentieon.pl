#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use JSON;# qw( encode_json );

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# General QC-script for sentieon-data. Takes two arguments: SAMPLE-ID, TYPE(panel or wgs)
# PANEL: Requires the following algos from sentieon qc:
#        --algo MeanQualityByCycle mq_metrics.txt \\
#        --algo QualDistribution qd_metrics.txt \\
#        --algo GCBias --summary gc_summary.txt gc_metrics.txt \\
#        --algo AlignmentStat aln_metrics.txt \\
#        --algo InsertSizeMetricAlgo is_metrics.txt \\
#        --algo HsMetricAlgo --targets_list $target_intervals --baits_list $target_intervals hs_metrics.txt \\
#        --algo CoverageMetrics --cov_thresh 1 --cov_thresh 10 --cov_thresh 30 --cov_thresh 100 --cov_thresh 250 --cov_thresh 500 cov_metrics.txt
# WGS: Requires the following algos from sentieon qc
#       --algo MeanQualityByCycle mq_metrics.txt \\
#        --algo QualDistribution qd_metrics.txt \\
#        --algo GCBias --summary gc_summary.txt gc_metrics.txt \\
#        --algo AlignmentStat aln_metrics.txt \\
#        --algo InsertSizeMetricAlgo is_metrics.txt \\
#        --algo WgsMetricsAlgo wgs_metrics.txt \\
#
# BOTH options also need the dedup_metrics from the locus_collector + dedup steps
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


my $SID = $ARGV[0];
my $type = $ARGV[1];

my %pct_above_x;
my $median;
my $pct_above_panel;
if ($type eq "panel") {

    ($pct_above_panel, $median) = coverage_calc();
    %pct_above_x = %$pct_above_panel;
}


my $align_metrics_file = "aln_metrics.txt";
my $insert_file = "is_metrics.txt";
my $dedup_metrics_file = "dedup_metrics.txt";

my %results;
my $metrics_file;
if ($type eq "wgs") {
    $metrics_file = "wgs_metrics.txt";
    open( HS, $metrics_file );
    while( <HS> ) {
        if( /^\#SentieonCommandLine/ ) {
	    <HS>;
            my $vals = <HS>;
            my @a = split /\t/, $vals;
            $results{'median_cov'} = $a[3];
            $results{'sd_coverage'} = $a[2];
            $results{'mean_coverage'} = $a[1];
            $pct_above_x{'1'} = $a[12];
            $pct_above_x{'5'} = $a[13];
            $pct_above_x{'10'} = $a[14];
            $pct_above_x{'15'} = $a[15];
            $pct_above_x{'20'} = $a[16];
            $pct_above_x{'25'} = $a[17];
            $pct_above_x{'30'} = $a[18];
            $pct_above_x{'40'} = $a[19];
            $pct_above_x{'50'} = $a[20];
            $pct_above_x{'60'} = $a[21];
            $pct_above_x{'70'} = $a[22];
            $pct_above_x{'80'} = $a[23];
            $pct_above_x{'90'} = $a[24];
            $pct_above_x{'100'} = $a[25];

	    }
    }
    close HS;

}
elsif ($type eq "panel") {
    $metrics_file = "hs_metrics.txt";
    open( HS, $metrics_file );
    while( <HS> ) {
        if( /^\#SentieonCommandLine/ ) {
	    <HS>;
            my $vals = <HS>;
            my @a = split /\t/, $vals;
            $results{'pct_on_target'} = $a[18];
            #print "pct_on_target: $a[18]\n";
            $results{'fold_enrichment'} = $a[25];
            #print "fold_enrichment: $a[25]\n";
            $results{'mean_coverage'} = $a[22];
            #print "median_coverage: $a[22]\n";
            $results{'fold_80'} = $a[32];
            #print "fold_80: $a[32]\n";
	    }
    }
    close HS;
    $results{'median_cov'} = $median;
}
else { print STDERR "panel or wgs\n"; exit;}


## INSSIZE ##
open( INS, $insert_file );
while( <INS> ) {
	if( /^\#SentieonCommandLine/ ) {
	    <INS>;
	    my $vals = <INS>;
	    my @a = split /\t/, $vals;
        my $ins_size = sprintf "%.0f", $a[4];
	    $results{'ins_size'} = $ins_size;
        #print "ins_size: $a[4]\n";
        my $ins_size_dev = sprintf "%.0f", $a[5];
	    $results{'ins_size_dev'} = $ins_size_dev;
        #print "ins_size_dev: $a[5]\n";
	}
}
close INS;


## DEDUP ##

open( DEDUP, $dedup_metrics_file );
while( <DEDUP> ) {
    if( /^\#SentieonCommandLine/ ) {
	    <DEDUP>;
	    my $vals = <DEDUP>;
	    my @a = split /\t/, $vals;
	    $results{'dup_reads'} = $a[6];
        #print "dup_reads: $a[6]\n";
	    $results{'num_reads'} = $a[2];
        #print "num_reads: $a[2]\n";
        $results{'dup_pct'} = $a[8];
        #print "dup_pct: $a[8]\n";
        my $mapped = $a[2]-$a[4];
        $results{'mapped_reads'} = $mapped;
        #print "mapped_reads: $mapped\n";
	}
}
close DEDUP;

sub coverage_calc {

    my $coverage_file_summary = "cov_metrics.txt.sample_summary";
    my $coverage_file = "cov_metrics.txt";
    my @cov;
    open( COV, $coverage_file );
    while( <COV> ) {

        unless( /^Locus/ ) {
            my $vals = <COV>;
            my @a = split /\t/, $vals;
            push @cov,$a[1];
        }
        
    }
    my @sorted_cov = sort @cov;

    my $len_cov = scalar(@cov);
    my $median;
    if (($len_cov / 2) =~ /\d+\.\d+/) { 
        my $median_index = ($len_cov / 2) + 0.5;
        $median = $sorted_cov[$median_index-1];
    }
    else {
        my $median_index = ($len_cov / 2);
        $median = ($sorted_cov[$median_index-1] + $sorted_cov[$median_index-1]) / 2;
    }


    ## cov_% ##
    my %pct_above_x;
    open( COV_THRESH, $coverage_file_summary );
    while( <COV_THRESH> ) {
        if( /^sample_id/ ) {
            my $vals = <COV_THRESH>;
            my @a = split /\t/, $vals;
            chomp @a;
            $pct_above_x{'1'} = $a[6];
            $pct_above_x{'10'} = $a[7];
            $pct_above_x{'30'} = $a[8];
            $pct_above_x{'100'} = $a[9];
            $pct_above_x{'250'} = $a[10];
            $pct_above_x{'500'} = $a[11];
        }
    }
    close COV_THRESH;
    return \%pct_above_x, $median;

}

$results{'pct_above_x'} = \%pct_above_x;
$results{'sample_id'} = $SID;
my $json = JSON->new->allow_nonref;
print $json->pretty->encode( \%results );
