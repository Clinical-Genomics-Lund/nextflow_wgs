use strict;
use CMD::vcf qw( parse_vcf );
use Data::Dumper;
use GD;
use List::Util qw( max min sum );
use POSIX 'strftime';


my @chr_order = qw(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y);
my %chr_size = read_chromosome_sizes("/data/bnf/ref/b37/human_g1k_v37.dict");
my $max_chr_size = max( values %chr_size );
my $date = strftime '%Y%m%d', localtime;

# Read and parse VCF
my( $meta, $vars, $samples ) = parse_vcf( $ARGV[0] );

my( $table_fn, $png_fn ) =  ( $ARGV[1].".tsv", $ARGV[1].".png" );

# Get IDs of proband, mother and father (only proband if single sample)
my( $proband_id, $mother_id, $father_id ) = ( $ARGV[2], $ARGV[3], $ARGV[4] );

my $sample_id_string = "Proband: $proband_id";
$sample_id_string   .= " Mor: $mother_id" if $mother_id;
$sample_id_string   .= " Far: $father_id" if $father_id;

my %data;
my %plot_data;
my $in_loh = 0;
my $loh_p = 1;
my $loh_tot_size = 0;
my $prev_chr = "";
my $prev_pos = 0;

foreach my $v ( @$vars ) {

    next unless $v->{INFO}->{CSQ}->[0]->{MAX_AF};

    if( $v->{INFO}->{CSQ}->[0]->{MAX_AF} > 0.05 ) {

	my( $chr, $pos ) = ( $v->{CHROM}, $v->{POS} );
	
	# TRIO, filter low quality positions
	my( $pb, $m, $f );
	if( $mother_id and $father_id ) {
	    ( $pb, $m, $f ) = ( $v->{GT}->{$proband_id}->{GT}, $v->{GT}->{$mother_id}->{GT}, $v->{GT}->{$father_id}->{GT} );
	    next if $v->{GT}->{$proband_id}->{GQ} eq "." or $v->{GT}->{$mother_id}->{GQ} eq "." or $v->{GT}->{$father_id}->{GQ} eq ".";	
	    next if $v->{GT}->{$proband_id}->{GQ} < 30 or $v->{GT}->{$mother_id}->{GQ} < 30 or $v->{GT}->{$father_id}->{GQ} < 30;
	    next if $pb =~ /\./ or $m =~ /\./ or $f =~ /\./;
	}

	# SINGLE filter low quality positions
	else {
	    $pb = $v->{GT}->{$proband_id}->{GT};
	    next if $v->{GT}->{$proband_id}->{GQ} eq "." or $v->{GT}->{$proband_id}->{GQ} < 30;
	    next if $pb =~ /\./;
	}

	# TRIO, call UPD-informative positions
	my $call = ".";
	if( $m and $f ) {
	    if ($pb eq "0/0" ) {
		$call = "UF" if $m !~ /0/;
		$call = "UM" if $f !~ /0/;
	    }
	    elsif ($pb eq "1/1" ) {
		$call = "UF" if $m !~ /1/;
		$call = "UM" if $f !~ /1/;
	    }
	    elsif ($pb eq "0/1" ) {
		$call = "NU" if ( $m eq "1/1" and $f eq "0/0" ) or ( $f eq "1/1" and $m eq "0/0" );
	    }
	}

	# Count call types, for table
	$data{$chr}->{$call}++;

	my $alt_AF = max( 0.00001, ( $v->{INFO}->{CSQ}->[0]->{MAX_AF} or 0 ) );

	# Save plot data
	$plot_data{$chr}->{$pos}->{upd} = $call;
	my $ref_AF = max( 0.00001, 1 - $alt_AF );

	$loh_p = 1 if $prev_chr ne $chr;
	$loh_p = min(1, $loh_p*10000) if $pb eq "0/1";
	$loh_p *= $ref_AF if $pb eq "0/0";
	$loh_p *= $alt_AF if $pb eq "1/1";
	$plot_data{$chr}->{$pos}->{lohp} = $loh_p;

	if( $loh_p < 0.000000001 ) {
	    unless( $in_loh ) {
		$in_loh = $pos;
	    }
	}
	else {
	    if ( $in_loh ) {
		unless( $prev_chr eq "X" ) {
		    $loh_tot_size += ( $prev_pos - $in_loh );
		}
		$in_loh = 0;
	    }
	}
	
	$prev_chr = $chr;
	$prev_pos = $pos;
    }
}

open( TABLE, '>'.$table_fn );
print TABLE "Proband: $proband_id\n";
print TABLE "Analysdatum: $date\n";

if( $mother_id and $father_id ) {
    # Print UPD table
    print TABLE "#Chromosome\tTotal\tUninformative\tUPD Mother ($mother_id)\tUPD Father ($father_id)\tanti-UPD\n";
    foreach my $chr ( @chr_order ) {
	my $tot = (($data{$chr}->{'.'} or 0)+($data{$chr}->{'UM'} or 0)+($data{$chr}->{'UF'} or 0)+($data{$chr}->{'NU'} or 0) or 1);

	my $uninformative_percent = sprintf("%.2f", 100*($data{$chr}->{'.'} or 0)/$tot);
	my $UM_percent =  sprintf("%.2f",100*($data{$chr}->{'UM'} or 0)/$tot);
	my $UF_percent = sprintf("%.2f",100*($data{$chr}->{'UF'} or 0)/$tot);
	my $NU_percent = sprintf("%.2f",100*($data{$chr}->{'NU'} or 0)/$tot);
	
	print TABLE $chr."\t$tot\t$uninformative_percent\% (n=".($data{$chr}->{'.'} or 0).")\t$UM_percent\% (n=".($data{$chr}->{'UM'} or 0).")\t$UF_percent\% (n=".($data{$chr}->{'UF'} or 0).")\t$NU_percent\% (n=".($data{$chr}->{'NU'} or 0).")\n";
    }
}


#print TABLE "\nLOH: $loh_tot_size (". 100*($loh_tot_size / ( sum(values %chr_size) - $chr_size{'X'} ) )."%)\n";
printf TABLE "\nLOH: %d (%.2f%%)", $loh_tot_size, 100*($loh_tot_size / ( sum(values %chr_size) - $chr_size{'X'} ) );



###################
# PLOT EVERYTHING #
###################


my( $x_ofs, $y_ofs ) = ( 30, 30 );
my( $chr_spacing, $chr_H ) = ( 30, 8 );
my( $W, $H ) = ( 1300, 700 );
my $im = new GD::Image($W+$x_ofs*2,$H+$y_ofs*2,1);
$im->alphaBlending( 1 );
$im->saveAlpha( 1 );
my $white   = $im->colorAllocate(255,255,255);
my $black   = $im->colorAllocate(0,0,0);
my $gray    = $im->colorAllocateAlpha(200,200,200,60);
my $red     = $im->colorAllocateAlpha(255,20,20,60);
my $blue    = $im->colorAllocateAlpha(20,20,255,60);

my $loh_border  = $im->colorAllocate(190,130,0);
my $loh_col     = $im->colorAllocate(250,140,15);


$im->filledRectangle( 0, 0, $W+$x_ofs*2, $H+$y_ofs*2, $white );
$im->string( gdMediumBoldFont, $x_ofs, 2, $sample_id_string . "    Analysdatum: ". $date, $black );
my %color = ( 'UM'=>$red, 'UF'=>$blue, '.'=>$gray, 'NU'=>$black );
my $scaleX = $W / $max_chr_size;
my $y_pos = 0;
foreach my $chr ( @chr_order ) {

    $im->string( gdSmallFont,
		 $x_ofs - 15, 
		 $y_ofs + $y_pos * $chr_spacing - 2, 
		 $chr, $black);

    my( $prev_pos, $prev_LOH_P );
    foreach my $pos ( sort {$a<=>$b} keys %{$plot_data{$chr}} ) {

	my $UPD = $plot_data{$chr}->{$pos}->{upd};
	my $LOH_P = $plot_data{$chr}->{$pos}->{lohp};
	$LOH_P = 1 unless $LOH_P;
	
	if( $prev_pos ) {
	    $im->line( $x_ofs + $prev_pos*$scaleX,
		       $y_ofs + $y_pos * $chr_spacing + $chr_H + 3 + (1 - min(1,log($prev_LOH_P)/-40) ) * 11,
		       $x_ofs + $pos*$scaleX,
		       $y_ofs + $y_pos * $chr_spacing + $chr_H + 3 + (1 - min(1,log($LOH_P)/-40) ) * 11,
		       $loh_border
		);
	    if( $LOH_P < 0.000000001 ) {
		$im->filledRectangle( $x_ofs + $prev_pos*$scaleX,
				      $y_ofs + $y_pos * $chr_spacing + $chr_H + 3 + (1 - min(1,log($prev_LOH_P)/-30) ) * 11+1,
				      $x_ofs + $pos*$scaleX,
				      $y_ofs + $y_pos * $chr_spacing + $chr_H + 14, $loh_col );
	    }
	}
	
	if( $UPD eq "." ) {
	    $im->filledRectangle( $x_ofs + $pos*$scaleX,
				  $y_ofs + $y_pos * $chr_spacing,
				  $x_ofs + $pos*$scaleX,
				  $y_ofs + $chr_H + $y_pos * $chr_spacing,
				  $color{$UPD} )
	}
	($prev_pos, $prev_LOH_P) = ( $pos, $LOH_P );
    }

    
    foreach my $pos ( sort keys %{$plot_data{$chr}} ) {
	my $UPD = $plot_data{$chr}->{$pos}->{upd};
	if( $UPD ne "." ) {
	    $im->filledRectangle( $x_ofs + $pos*$scaleX,
				  $y_ofs + $y_pos * $chr_spacing,
				  $x_ofs + $pos*$scaleX+1,
				  $y_ofs + $chr_H + $y_pos * $chr_spacing,
				  $color{$UPD} );    
	}
    }

    $im->rectangle( $x_ofs -1 ,
		    $y_ofs + $y_pos * $chr_spacing - 1,
		    $x_ofs + $chr_size{$chr} * $scaleX + 1,
		    $y_ofs + $chr_H + $y_pos * $chr_spacing + 1,
		    $black);
    
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

