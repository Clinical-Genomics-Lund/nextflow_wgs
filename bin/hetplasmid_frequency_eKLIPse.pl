#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;

our %opt;
GetOptions( \%opt, 'bam=s', 'in=s', 'help' );

print_usage("ERROR: --bam,--in required!\n") unless $opt{bam} and $opt{in};

unless( -e $opt{bam} ) { print "ERROR: BAM file not found\n"; exit; }

unless( -e $opt{in} ) { print "ERROR: infile not found\n"; exit; }

if( defined($opt{bam}) and  defined($opt{in}) ){

    #my %eKLIPse;
    open OUT, ">hetplasmid_frequency.txt" or die "cannot write to file!";
    print OUT "Sample\t5_break\t3_break\t5_blast\t3_blast\tdel_cov\trest_cov\tfrequency\n";
    my $e_data = &process_eKLIPse( $opt{in} );

    for my $var ( sort keys %{$e_data} ){
        my ( %data, $res, $del_mc, $rest_mc );
        if( $$e_data{$var}{'5 breakpoint'} >  $$e_data{$var}{'3 breakpoint'} ){
            $res = process_sambamba( 'pre', $$e_data{$var}{'3 breakpoint'}, $$e_data{$var}{'5 breakpoint'} ); @data{ keys %$res } = values %$res;
            $res = process_sambamba( 'del', $$e_data{$var}{'5 breakpoint'}, $$e_data{$var}{'3 breakpoint'} ); @data{ keys %$res } = values %$res;
            $del_mc = sprintf( "%.0f", ($data{del}{meanCoverage} / 2) );
            $rest_mc = sprintf( "%.0f", $data{pre}{meanCoverage} );
        }else{
            $res = process_sambamba( 'pre', 1, $$e_data{$var}{'5 breakpoint'} ); @data{ keys %$res } = values %$res;
            $res = process_sambamba( 'del', $$e_data{$var}{'5 breakpoint'}, $$e_data{$var}{'3 breakpoint'} ); @data{ keys %$res } = values %$res;
            $res = process_sambamba( 'post', $$e_data{$var}{'3 breakpoint'}, '-1' ); @data{ keys %$res } = values %$res;
            $del_mc = sprintf( "%.0f", $data{del}{meanCoverage} );
            $rest_mc = sprintf("%.0f", (($data{pre}{meanCoverage} + $data{post}{meanCoverage})/2) );
        }
        print OUT "$data{pre}{sampleName}\t".( $$e_data{$var}{'5 breakpoint'} + 1 )."\t".( $$e_data{$var}{'3 breakpoint'} - 1 )."\t";
        print OUT "$$e_data{$var}{'5 Blast'}\t$$e_data{$var}{'3 Blast'}\t$del_mc\t$rest_mc\t".sprintf( "%.2f", ( 1 - ($del_mc/$rest_mc) ) * 100 )."\n";
        #print Dumper( \%data );
    }
    close OUT;
}

sub process_eKLIPse{
    my ( $in ) = @_;
    my ( %e_data, @labels );
    open IN, "$in" or die "cannot read infile!";
    while ( my $line = <IN> ){
        chomp( $line );
        $line =~s/\"//g;
        $line =~s/\x27//g; 
        my @cols = split /;/, $line;
        if( $line =~ /Title/ ){
            @labels = @cols;
        }else{
            for my $index(0..$#labels ){ $e_data{ "$cols[ 1 ]-$cols[ 2 ]" }{ $labels[ $index ] } = $cols[ $index ]; }
        }
    }
    return \%e_data;
} # process_eKLIPse

sub process_sambamba {
    my ( $order, $start, $end ) = @_;
    my ( %data, @labels );
    if( $order eq 'del' and $start > $end ){ # special case
        open BED, ">regions.bed" or die "cannot write to file!";
        print BED "M\t$start\t-1\n";
        print BED "M\t1\t$end\n";
        close BED;
        my @depth = `sambamba depth region -L regions.bed $opt{bam}`;
        for my $line( @depth ){
            chomp( $line );
            my @cols = split /\t/, $line;
            if( $line =~ /^#/ ){
                @labels = @cols;
            }else{
                $data{ $order }{ 'chromStart' } = $start;
                $data{ $order }{ 'chromEnd' } = $end;
                for my $index(0..$#labels ){
                    if( $labels[ $index ] eq '# chrom' ){ $data{ $order }{ $labels[ $index ] } = $cols[ $index ]; }
                    if( $labels[ $index ] eq 'sampleName' ){ $data{ $order }{ $labels[ $index ] } = $cols[ $index ]; }
                    if( $labels[ $index ] eq 'readCount' ){ $data{ $order }{ $labels[ $index ] } += $cols[ $index ]; }
                    if( $labels[ $index ] eq 'meanCoverage' ){ $data{ $order }{ $labels[ $index ] } += $cols[ $index ]; } 
                }
            }
        }
    }else{
        my @depth = `sambamba depth region -L M:$start-$end $opt{bam}`;
        for my $line( @depth ){
            chomp( $line );
            my @cols = split /\t/, $line;
            if( $line =~ /^#/ ){
                @labels = @cols;
            }else{
                for my $index(0..$#labels ){ $data{ $order }{ $labels[ $index ] } = $cols[ $index ]; }
            }
        }
    }
    return \%data;
} # process_sambamba


sub print_usage {
    print "$_[0]\n\n" if $_[0];
    print "USAGE: hetplasmid_frequency_eKLIPse.pl --bam <MITO BAM FILE> --in<eKLIPse output>\n\n";
    print "    --bam     FILE       Path to BAM file \(required\)\n\n";
    print "    --in      FILE       eKLIPse deletions outfile \(required\)\n\n";
    print "    --help               Will print this message\n\n";
    exit(0);
} # print_usage
