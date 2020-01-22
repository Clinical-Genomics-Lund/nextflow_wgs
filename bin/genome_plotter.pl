#!/usr/bin/env perl
use warnings;
use strict;
use Data::Dumper;
use GD;
use List::Util qw( max min sum );
use POSIX 'strftime';
use Getopt::Long;

my $MIN_ROH_QUAL = 85;

my $date = strftime '%Y%m%d', localtime;

my( $roh_in, $upd_in, $cov_in, $sample_id, $png_fn, $dict_in );
$png_fn = "out.png";
GetOptions( "roh=s" => \$roh_in,
	    "upd=s" => \$upd_in,
	    "cov=s" => \$cov_in,
	    "dict=s" => \$dict_in,
	    "sample=s" => \$sample_id,
	    "out=s" => \$png_fn);

die "--sample must be given" unless $sample_id;
die "--dict must be given" unless $dict_in;

my @chr_order = qw(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y);
my %chr_size = read_chromosome_sizes($dict_in);
my $max_chr_size = max( values %chr_size );


my( %roh_data, %upd_data, %cov_data );


if( $roh_in ) {
    print STDERR "Reading ROH...";
    %roh_data = read_roh($roh_in);
}

if( $upd_in ) {
    print STDERR "\nReading UPD...";
    %upd_data = read_upd($upd_in);
}

if( $cov_in ) {
    print STDERR "\nReading coverage data...";
    %cov_data = read_cov($cov_in);
}

print STDERR "\nPlotting...";
###################
# PLOT EVERYTHING #
###################


my( $x_ofs, $y_ofs ) = ( 30, 30 );
my( $chr_spacing, $chr_H ) = ( 95, 10 );
my( $W, $H ) = ( 4000, 2300 );
my $im = new GD::Image($W+$x_ofs*2,$H+$y_ofs*2,1);
$im->alphaBlending( 1 );
$im->saveAlpha( 1 );
my $white   = $im->colorAllocate(255,255,255);
my $black   = $im->colorAllocate(0,0,0);
my $gray    = $im->colorAllocateAlpha(220,220,220,20);
my $covmid  = $im->colorAllocateAlpha(20,20,20,40);
my $orange  = $im->colorAllocateAlpha(255,165,00,30);
my $red     = $im->colorAllocateAlpha(255,20,20,30);
my $blue    = $im->colorAllocateAlpha(20,20,255,30);
my $lblue    = $im->colorAllocateAlpha(120,120,255,30);

my $loh_border  = $im->colorAllocate(190,130,0);
my $loh_col     = $im->colorAllocate(250,140,15);

my $scaleX = $W / $max_chr_size;

$im->filledRectangle( 0, 0, $W+$x_ofs*2, $H+$y_ofs*2, $gray );
my $y_pos = 0;
my $roh = $roh_data{$sample_id};
foreach my $chr ( @chr_order ) {

    $im->string( gdLargeFont,
		 $x_ofs - 20, 
		 $y_ofs + $y_pos * $chr_spacing + 2, 
		 $chr, $black);

    $im->filledRectangle( $x_ofs -1 ,
		    $y_ofs + $y_pos * $chr_spacing - 1,
		    $x_ofs + $chr_size{$chr} * $scaleX + 1,
		    $y_ofs + $y_pos * $chr_spacing + $chr_spacing - 15,
		    $white);

    # Plot ROH regions
    for my $rohreg ( @{ $roh->{$chr} } ) {

	$im->filledRectangle( $x_ofs + $rohreg->{start} * $scaleX,
			$y_ofs + $y_pos * $chr_spacing + 3,
			$x_ofs + $rohreg->{end} * $scaleX,
			$y_ofs + $y_pos * $chr_spacing + 17,
			$orange );

    }

    # Plot UPD regions
    for my $updreg ( @{ $upd_data{$chr} } ) {
	my $info = $updreg->{info};
	next if $info->{INFORMATIVE_SITES} < 100;
	$im->filledRectangle( $x_ofs + $updreg->{start} * $scaleX,
			$y_ofs + $y_pos * $chr_spacing + 61,
			$x_ofs + $updreg->{end} * $scaleX + 1,
			$y_ofs + $y_pos * $chr_spacing + 75,
			$info->{ORIGIN} eq "PATERNAL" ? $blue : $red );

    }


    my $cov_ypos = 38;
    my $cov_height = 21;
    if( !$upd_in and !$roh_in ) {
	$cov_height = 30;
    }
    # Plot coverage
    for my $c ( @{ $cov_data{$chr} } ) {
	$im->setPixel( $x_ofs + $c->{pos} * $scaleX,
		       $y_ofs + $y_pos * $chr_spacing + $cov_ypos- $c->{cov}*$cov_height,
		       $lblue );

    }
    $im->line( $x_ofs, 
	       $y_ofs + $y_pos * $chr_spacing + $cov_ypos,
	       $x_ofs + $chr_size{$chr}*$scaleX,
	       $y_ofs + $y_pos * $chr_spacing + $cov_ypos,
	       $covmid);
    
    $y_pos++;
}
    


open( PNG, ">".$png_fn );
print PNG $im->png;



sub read_chromosome_sizes {
    my $fn = shift;
    open( DICT, $fn );
    my %size;
    while(<DICT>) {
	if( /^\@SQ.*?SN:(.*?)\tLN:(.*?)\t/ ) {
	    $size{$1} = $2;
	}
    }
    return %size;
}


#RG      1443-14 1       827277  931748  104472  293     67.9
sub read_roh {
    my $fn = shift;

    open(ROH, $fn );
    my %roh;
    while(<ROH>) {
	next unless /^RG/;
	my( $rg, $sample, $chr, $start, $end, $size, $snps, $score ) = split /\t/;
	if( $score > $MIN_ROH_QUAL ) {
	    push( @{$roh{$sample}->{$chr}},  {'start'=>$start, 'end'=>$end, 'snps'=>$snps, 'score'=>$score} );
	}
    }
    close ROH;

    return %roh;
}

sub read_upd {
    my $fn = shift;

    open(UPD, $fn);
    my %upd;
    while(<UPD>) {
	my( $chr, $start, $end, $data ) = split /\t/;
	my %info;
	for my $pair (split /;/, $data) {
	    my( $key, $val ) = split /=/, $pair;
	    $info{$key} = $val;
	}
	push @{$upd{$chr}},  {'start'=>$start, 'end'=>$end, 'info'=>\%info};
    }
    return %upd;
}

sub read_cov {
    my $fn = shift;

    open(COV, $fn);

    my %cov;
    my $cnt;
    my $cn_sum;
    my $step = 11;
    my @pos;
    while( <COV> ) {
	$cnt++;
	next if /^\@/ or /^CONTIG/;
	my( $chr, $start, $stop, $cn ) = split /\t/;
	push( @pos, $start );
	$cn_sum += $cn;

	if( $cnt % $step == 0 ) {
	    push @{$cov{$chr}}, { 'pos'=>$pos[int $step/2], 'cov'=> $cn_sum/$step};
	    @pos = ();
	    $cn_sum = 0;
	}
    }
    return %cov;
}
