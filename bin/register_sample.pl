#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use MongoDB;
use MongoDB::BSON;
use MongoDB::OID;
use DateTime;
use Data::Dumper;
use JSON::Parse 'json_file_to_perl';


# Connect to mongodb
my $client = MongoDB->connect();
 
my $DB_FILES = $client->ns("CMD.files");



my %opt;
GetOptions( \%opt, 
	    'sample-id=s', 'run-folder=s', 'lims-id=s', 'fq-se=s', 'fq-fwd=s', 'fq-rev=s', 'bam-desc=s', 
	    'bam-file=s', 'vcf-desc=s', 'vcf-file', 'overwrite', 'assay=s', 'desc=s', 'subassay=s',
            'family=s', 'qc=s', 'sex=s', 'contamination=s' );

print_usage() unless $opt{'sample-id'} or $opt{'run-folder'} or $opt{'assay'};


# Ugly hack to fix assay name of myeloid samples 
$opt{'assay'} = "myeloid" if $opt{'assay'} eq "myeloid-nextera" or $opt{'assay'} eq "myeloid-panel";


my $OVERWRITE = 1 if $opt{'overwrite'};

# Add fq-data
if( $opt{'fq-se'} and ( $opt{'fq-fwd'} or $opt{'fq-rev'} ) ) {
    print_usage("ERROR: Single-end (--fq-se) and paired-end (--fq-fwd/--fq-rev) cannot be combined!");
}
elsif( $opt{'fq-fwd'} xor $opt{'fq-rev'} ) {
    print_usage( "ERROR: Both --fq-fwd and --fq-rev must be specified. For single-end data, use --fq-se!");
}
elsif( $opt{'fq-se'} ) {
    add_se_fq( $opt{'sample-id'}, $opt{'run-folder'}, $opt{'assay'}, $opt{'fq-se'} );
}
elsif( $opt{'fq-fwd'} and $opt{'fq-rev'} ) {
    add_pe_fq( $opt{'sample-id'}, $opt{'run-folder'}, $opt{'assay'}, $opt{'fq-fwd'}, $opt{'fq-rev'} )
}

# Add clarity ID
if( $opt{'lims-id'} ) {
    add_lims_id( $opt{'sample-id'}, $opt{'run-folder'}, $opt{'assay'}, $opt{'lims-id'} );
}

# Add description
if( $opt{'desc'} ) {
    add_description( $opt{'sample-id'}, $opt{'run-folder'}, $opt{'assay'}, $opt{'desc'} );
}

# Add family
if( $opt{'family'} ) {
    add_family( $opt{'sample-id'}, $opt{'run-folder'}, $opt{'assay'}, $opt{'family'} );
}

# Add subassay
if( $opt{'subassay'} ) {
    add_subassay( $opt{'sample-id'}, $opt{'run-folder'}, $opt{'assay'}, $opt{'subassay'} );
}

# Add sex
if( $opt{'sex'} ) {
    add_sex( $opt{'sample-id'}, $opt{'run-folder'}, $opt{'assay'}, $opt{'sex'} );
}

# Add contaminations value
if( $opt{'contamination'} ) {
    add_contamination_value( $opt{'sample-id'}, $opt{'run-folder'}, $opt{'assay'}, $opt{'contamination'} );
}

# Add QC-data
if( $opt{'qc'} ) {
    add_qc( $opt{'sample-id'}, $opt{'run-folder'}, $opt{'assay'}, $opt{'qc'} );
}

sub print_usage {
    print $_[0]."\n\n" if $_[0];
    print "USAGE: coming soon...\n";
    exit(0);
}



sub add_se_fq {
    my( $id, $run, $assay, $fq ) = @_;

    my @fq = split /,/, $fq;

    # Find the sample from this run, if it already exists
    my $results = $DB_FILES->find_one( {'sample-id'=>$id, 'run'=>$run } );

    # Check if fq-se data already exists for this sample in this run (ie. duplicate add)
    my %exists;
    foreach( @{ $results->{'fq-se'} } ) {
	$exists{$_} = 1;
    }

    # Add files to data
    my @fq_to_add;
    foreach( @fq ) {
	if( $exists{$_} and !$OVERWRITE ) {
	    print STDERR "WARNING: SE fastq-file $_ is already added to this sample/run. Not adding again!\n";
	}
	else {
	    push @fq_to_add, $_;
	}
    }

    if( @fq_to_add ) {
	if( !$OVERWRITE ) {
	    $DB_FILES->update_one( {'sample-id'=>$id, 'run'=>$run, 'assay'=>$assay }, { '$set'=>{'date_modified' => DateTime->now}, '$push'=> {'fq_se'=>\@fq_to_add}}, {'upsert'=>1} );
	}
	else {
	    $DB_FILES->update_one( {'sample-id'=>$id, 'run'=>$run, 'assay'=>$assay }, { '$set'=>{'date_modified' => DateTime->now, 'fq_se'=>\@fq_to_add} }, {'upsert'=>1} );	    
	}
    }
    
}


sub add_pe_fq {
    my( $id, $run, $assay, $fq_fwd, $fq_rev ) = @_;

    my @fq_fwd = split /,/, $fq_fwd;
    my @fq_rev = split /,/, $fq_rev;

    if( scalar @fq_fwd != scalar @fq_rev ) {
	print STDERR "ERROR: Different number of files for fwd and rev\n";
	exit();
    }

    # Find the sample from this run, if it already exists
    my $results = $DB_FILES->find_one( {'sample-id'=>$id, 'run'=>$run } );

    # Check if fq-se data already exists for this sample in this run (ie. duplicate add)
    my %exists;
    foreach( @{ $results->{'fq_fwd'} } ) {
	$exists{$_} = 1;
    }
    foreach( @{ $results->{'fq_rev'} } ) {
	$exists{$_} = 1;
    }

    # Add files to data
    my( @fq_to_add_fwd, @fq_to_add_rev );
    foreach( @fq_fwd ) {
	if( $exists{$_} and !$OVERWRITE ) {
	    print STDERR "ERROR: PE fastq-file $_ is already added to this sample/run!\n";
	    exit;
	}
	else {
	    push @fq_to_add_fwd, $_;
	}
    }
    foreach( @fq_rev ) {
	if( $exists{$_} and !$OVERWRITE ) {
	    print STDERR "ERROR: PE fastq-file $_ is already added to this sample/run!\n";
	    exit;
	}
	else {
	    push @fq_to_add_rev, $_;
	}
    }


    if( @fq_to_add_rev and @fq_to_add_fwd ) {
	if( !$OVERWRITE ) {
	    $DB_FILES->update_one( {'sample-id'=>$id, 'run'=>$run, 'assay'=>$assay }, { '$set'=>{'date_modified' => DateTime->now}, '$push'=>{'fq_fwd'=>{'$each'=>\@fq_to_add_fwd}, 'fq_rev'=>{'$each'=>\@fq_to_add_rev}}}, {'upsert'=>1} );
	}
	else {
	    $DB_FILES->update_one( {'sample-id'=>$id, 'run'=>$run, 'assay'=>$assay }, { '$set'=>{'fq_fwd'=>\@fq_to_add_fwd, 'fq_rev'=>\@fq_to_add_rev, 'date_modified' => DateTime->now} }, {'upsert'=>1} );
	}
    }
    
}


sub add_lims_id {
    my( $id, $run, $assay, $lims_id ) = @_;

    # Find the sample from this run, if it already exists
    my $results = $DB_FILES->find_one( {'sample-id'=>$id, 'run'=>$run } );

    if( $results->{'lims-id'} and !$OVERWRITE ) {
	print STDERR "ERROR: LIMS-ID already exists. Use --overwrite\n";
	exit;	    
    }
    else {
	$DB_FILES->update_one( {'sample-id'=>$id, 'run'=>$run, 'assay'=>$assay }, { '$set'=>{'lims-id'=>$lims_id}}, {'upsert'=>1} );
    }
}


sub add_description {
    my( $id, $run, $assay, $desc) = @_;

    # Find the sample from this run, if it already exists
    my $results = $DB_FILES->find_one( {'sample-id'=>$id, 'run'=>$run } );

    if( $results->{'desc'} and !$OVERWRITE ) {
	print STDERR "ERROR: Description already exists. Use --overwrite\n";
	exit;	    
    }
    else {
	$DB_FILES->update_one( {'sample-id'=>$id, 'run'=>$run, 'assay'=>$assay }, { '$set'=>{'desc'=>$desc}}, {'upsert'=>1} );
    }
}

sub add_family {
    my( $id, $run, $assay, $family) = @_;

    # Find the sample from this run, if it already exists
    my $results = $DB_FILES->find_one( {'sample-id'=>$id, 'run'=>$run } );

    if( $results->{'family'} and !$OVERWRITE ) {
	print STDERR "ERROR: Family ID already exists. Use --overwrite\n";
	exit;	    
    }
    else {
	$DB_FILES->update_one( {'sample-id'=>$id, 'run'=>$run, 'assay'=>$assay }, { '$set'=>{'family'=>$family}}, {'upsert'=>1} );
    }
}

sub add_subassay {
    my( $id, $run, $assay, $subassay) = @_;

    # Find the sample from this run, if it already exists
    my $results = $DB_FILES->find_one( {'sample-id'=>$id, 'run'=>$run } );

    if( $results->{'subassay'} and !$OVERWRITE ) {
	print STDERR "ERROR: Subassay already exists. Use --overwrite\n";
	exit;	    
    }
    else {
	$DB_FILES->update_one( {'sample-id'=>$id, 'run'=>$run, 'assay'=>$assay }, { '$set'=>{'subassay'=>$subassay}}, {'upsert'=>1} );
    }
}


sub add_sex {
    my( $id, $run, $assay, $sex ) = @_;

    # Find the sample from this run, if it already exists
    my $results = $DB_FILES->find_one( {'sample-id'=>$id, 'run'=>$run } );

    if( $results->{'sex'} and !$OVERWRITE ) {
	print STDERR "ERROR: Sex already exists. Use --overwrite\n";
	exit;	    
    }
    else {
	$DB_FILES->update_one( {'sample-id'=>$id, 'run'=>$run, 'assay'=>$assay }, { '$set'=>{'sex'=>$sex}}, {'upsert'=>1} );
    }
}

sub add_contamination_value {
    my( $id, $run, $assay, $contamination ) = @_;

    # Find the sample from this run, if it already exists
    my $results = $DB_FILES->find_one( {'sample-id'=>$id, 'run'=>$run } );

    if( $results->{'contamination'} and !$OVERWRITE ) {
	print STDERR "ERROR: Contamination value already exists. Use --overwrite\n";
	exit;	    
    }
    else {
	$DB_FILES->update_one( {'sample-id'=>$id, 'run'=>$run, 'assay'=>$assay }, { '$set'=>{'contamination'=>$contamination}}, {'upsert'=>1} );
    }
}



sub add_qc {
    my( $id, $run, $assay, $json_fn) = @_;

    my %qc = %{json_file_to_perl ($json_fn)};
    
    my $results = $DB_FILES->find_one( {'sample-id'=>$id, 'run'=>$run } );
    if( $results->{'json'} and !$OVERWRITE ) {
	print STDERR "ERROR: QC data already exists. Use --overwrite\n";
	exit;	
    }
    else {
	$DB_FILES->update_one( {'sample-id'=>$id, 'run'=>$run, 'assay'=>$assay }, { '$set'=>{'qc'=>\%qc }}, {'upsert'=>1} );
    }
}





