#!/usr/bin/perl -w
use strict;
use MongoDB;


# Script that queries loqusdb, must be run on cmdscout1
# Annotates vcf with loqusdb frequency for each variant 
# in info-field as loqusdb_freq=

open(VCF, $ARGV[0]);
my @SAMPLE;
# Save all variants in vcf into array
while ( <VCF>   ) {
    if (/^#/) {
    }
    else {
        my @coords = split /\t/;
        my $vcf_var = join "_", @coords[0..1,3..4];
       # print $vcf_var,"\n";
        push @SAMPLE,$vcf_var;
    }


}
close VCF;

# Connect to mongodb
# Set default path
my $host = "mongodb://cmdscout2.lund.skane.se/loqusdb";

# If on hopper, reverse tunnel needs to bu used. Ports are gotten from an env var, given by second argument
# Do not use this feature anywhere else than on Hopper!
if( $ARGV[1] ) {
    if( $ENV{$ARGV[1]} ) {
	my $port = $ENV{$ARGV[1]};
	$host = "mongodb://localhost:$port/loqusdb";
    }
    else {
	die "No port envvar set for $ARGV[1]";
    }
}
my $client = MongoDB->connect($host);

my $CASES = $client->ns("loqusdb.case");
my $VARS = $client->ns("loqusdb.variant");

my $num_cases = $CASES->count;
print STDERR "$num_cases Cases in loqusdb\n";
my $cases = $CASES->find();
my $variants = $VARS->find( {'_id'=> { "\$in" => [@SAMPLE] }} );

# Count observations in locus_db
my %obs_per_var;
while( my $var = $variants->next ) {
    $obs_per_var{$var->{'_id'}} = $var->{'observations'};
}

open(VCF, $ARGV[0]);
my $c = 0;
while ( <VCF>   ) {
    if (/^#/) {
        print;
        $c++;
        if ( $c == 10 ) {
            print "##INFO\=<ID=loqusdb_freq\,Number=1\,Type=Float,Description=\"loqusdb allele frequency\">\n";
        }
    }
    else {
        my @coords = split /\t/;
        my $vcf_var = join "_", @coords[0..1,3..4];
        print join "\t", @coords[0..6];
        my $freq;
        print "\t";
        if ($obs_per_var{$vcf_var}) {
            $freq = $obs_per_var{$vcf_var}/$num_cases;
            $freq = 'loqusdb_freq='.$freq;
            print $coords[7].";".$freq;
            print "\t";
        }
        else {
            print $coords[7];
            print "\t";
        }

        #print everything after info field
        print join "\t", @coords[8..$#coords];

    }



}
