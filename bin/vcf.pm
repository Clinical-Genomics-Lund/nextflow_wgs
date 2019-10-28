package vcf;
use strict;

use Exporter qw(import);
our @EXPORT_OK = qw( parse_vcf parse_metainfo parse_variant parse_genotype parse_info ); 



# Parse VCF file and return
sub parse_vcf {
    
    my $fn = shift;
    my %vcf_meta;
    my @vcf_data;
    my @head ;
    my $full_header_str;

    my $vcf_fh;
    if( is_gzipped($fn) ) {
#    if( $fn =~ /\.gz/ ) {
	open( $vcf_fh, "zcat $fn |" );
    }
    else {
	open( $vcf_fh, $fn );
    }
	
    while( <$vcf_fh> ) {
	chomp;

	# Skip empty lines
	next if /^\s*$/;

	$full_header_str .= $_."\n" if /^#/;
	
	# Header row with meta data.
	if( /^##/ ) {
	    my( $type, $meta ) = parse_metainfo( $_ );
	    $vcf_meta{$type}->{$meta->{ID}} = $meta if defined $type;
	}

	# Header with column description.
	elsif( /^#/ ) {
	    #print STDERR $_;
	    $_ =~ s/^#//;
	    @head = split /\t/;
	}

	# Actual variant data
	else {
	    die "Malformed VCF: No column description header." unless @head;
	    my $variant = parse_variant( $_, \@head, \%vcf_meta );
	    if( $variant->{CHROM} ) {
		push @vcf_data, $variant if defined $variant;
	    }
	}
    }

    my @samples = splice @head, 9;
    $vcf_meta{header_str} = $full_header_str;
    
    return \%vcf_meta, \@vcf_data, \@samples;
}




# Parse VCF meta info line (only FORMAT and INFO)
sub parse_metainfo {
    my $comment = shift;
    
    $comment =~ s/^##//;
    my( $type, $data ) = ( $comment =~ /^(.*?)=(.*)$/ );


    if( $type eq "FORMAT" or $type eq "INFO" or $type eq "SAMPLE" or $type eq "FILTER" ) {
	$data = remove_surrounding( $data, '<', '>' );
	my $pairs = keyval( $data, '=', ',' );
	return $type, $pairs;
    }

    return undef, undef;
}


# Parse VCF variant line
sub parse_variant {
    my( $var_str, $head, $meta ) = @_;
    my @var_data = split /\t/, $var_str;
    my %var;

    $var{ vcf_str } = $var_str;
    
    # First seven fields
    for ( 0..6 ) {
	$var{ $head->[$_] } = $var_data[$_];
    }

    # Eigth field, INFO
    $var{ INFO } = parse_info( $var_data[7] );

    # Parse VEP annotation field, if any
    if( $var{ INFO }->{ CSQ } ) {
	$var{ INFO }->{ CSQ } = parse_VEP_CSQ( $var{INFO}->{CSQ}, $meta->{INFO}->{CSQ} );
    }
    
    # Genotypes for each sample
    for ( 9 .. (@var_data-1) ) {
	$var{ GT } -> { $head->[$_] } = parse_genotype( $var_data[8], $var_data[$_] );
    }

    return \%var;
}


# Parse genotype field of VCF
sub parse_genotype {
    my( $format, $data ) = @_;

    my @format = split ':', $format;
    my @data   = split ':', $data;

    my %gt;
    @gt{@format} = @data;

    return \%gt;
}


# Parse info column of VCF file
sub parse_info {
    my $str = shift;
    my $info = keyval( $str, "=", ";" );

    return $info;
}


sub parse_VEP_CSQ {
    my( $CSQ_var, $CSQ_meta ) = @_;

    $CSQ_meta->{Description} =~ /Consequence annotations from Ensembl VEP\. Format: (.*?)$/;

    my @field_names = split /\|/, $1;
    
    my @transcripts = split /,/, $CSQ_var;

    my @data_transcripts;
    foreach my $transcript_CSQ ( @transcripts ) {
	my @values = split /\|/, $transcript_CSQ;

	my %data;
	for( 0 .. $#field_names ) {
	    if( $field_names[$_] eq "Consequence" ) {
		my @conseq_array = split '&', $values[$_];
		$data{ $field_names[$_] } = \@conseq_array;
	    }
	    else {
		$data{ $field_names[$_] } = ( $values[$_] or "" );
	    }
	}

	push( @data_transcripts, \%data )
    }
    return \@data_transcripts;
}

# Removes character(s) defined in arg2 if first in string, and arg3 if last in string.
sub remove_surrounding {
    my( $str, $before, $after ) = @_;
    $str =~ s/^$before//;
    $str =~ s/$after$//;
    return $str;
}


# Parse string with key value pairs. Return hash. 
#  * Keys and values separated by 2nd argument. 
#  * Pairs separated by 3rd argument
#  * Handles commas in values if surrounded by double quotes
sub keyval {
    my( $str, $keyval_sep, $pair_sep ) = @_;

    my @pair_str = split /$pair_sep/, $str;
    my %pairs;
    foreach( @pair_str ) {

	# If key-value separator exists, save the value for the key
	if( /$keyval_sep/ ) {
	    my( $key, $val ) = split /$keyval_sep/;
	    $pairs{$key} = $val;
	}

	# Otherwise treat the whole string as a flag and set it to one (true).
	else {
	    $pairs{$_} = 1;
	}
    }
    return \%pairs;
}

sub excel_float {
    my $val = shift;

    return 0 if $val eq ".";
    
    $val =~ s/\./,/;
    return $val;
}

sub is_gzipped {
    my $fn = shift;

    my $file_str = `file -L $fn`;
    return 1 if $file_str =~ /gzip compressed/;
    return 0;
}


1;
