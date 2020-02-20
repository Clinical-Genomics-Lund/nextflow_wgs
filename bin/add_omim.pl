#!/usr/bin/env perl
use strict;
use warnings;
use File::Basename qw( basename dirname );
use File::Spec;
use Data::Dumper;

my $omim_db = File::Spec->rel2abs(dirname($0))."/omim_genes.txt";
my %omim = get_omim($omim_db);


my $hgncid_idx;
while(<>) {
    if( /^##INFO=<ID=CSQ/ ) {
        my( $vep_fields_str ) = ( $_ =~ /Format: (.*?)">/ );
        my @vep_fields = split /\|/, $vep_fields_str;
        ($hgncid_idx) = grep { $vep_fields[$_] eq 'HGNC_ID' } (0..@vep_fields-1);
    }
    if( /^#CHROM/ ) {
        print "##INFO=<ID=OMIM_GENES,Number=.,Type=String,Description=\"Overlapping OMIM genes\">\n";
    }
    if( /^#/ ) {
        print $_;
    }
    else {
        my @v = split /\t/;
        my $info = $v[7];
        my @info = split /;/, $info;
        my @omim_matches;

        foreach my $i ( @info ) {
            if( $i =~ /^CSQ=(.*)$/ ) {
                my @csqs = split /,/, $1;
                foreach my $csq ( @csqs ) {
                    my @fields = split /\|/, $csq;
                    my $hgnc_id = $fields[$hgncid_idx];
                    if($hgnc_id and $omim{$hgnc_id} ) {
                        push( @omim_matches, "$omim{$hgnc_id}" );
                    }
                }
            }
        }
        if( @omim_matches ) {
            $v[7] .= ";OMIM_GENES=".join(",", @omim_matches);
        }
        
        print join("\t", @v);
    }
}




sub get_omim {
    my $fn = shift;
    open( OMIM, $fn );
    my %omim;
    while(<OMIM>) {
        chomp;
        my( $symbol, $hgnc ) = split /\t/;
        $omim{$hgnc} = $symbol;
    }
    close OMIM;
    return %omim;
}
