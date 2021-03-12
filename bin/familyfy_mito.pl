#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename;
use lib dirname (__FILE__);
use vcf2;
use strict;
use Data::Dumper;


my %opt = ();
GetOptions( \%opt, 'vcfs=s', 'proband=s', 'mother=s', 'father=s', 'out=s' );
my $proband = $opt{proband};
my $mother = $opt{mother};
my $father = $opt{father};
my $out = $opt{out};

my @vcfs = check_options( \%opt );
my @extra_individuals = ($proband);
if ($father ne "null") {
	push @extra_individuals,$father;
}
if ($mother ne "null") {
	push @extra_individuals,$mother;
}
my $test = aggregate(@vcfs);
#print Dumper($test);
foreach my $variant (keys %$test) {
	my $print_vcf = vcfstr($test->{$variant},[]);
	print $print_vcf;
}


# ## Add the relatives to VCF-header
# foreach my $line (@header) {
# 	if ($line =~ /^#CHROM/) {
# 		$line = $line."\t".join("\t",@extra_individuals);
# 	}
# 	print OUT $line,"\n";
# }


sub aggregate {
    my @vcfs = @_;
    my %agg;
    foreach my $fn (@vcfs) {
        my $vcf = CMD::vcf2->new( 'file'=>$fn );
        my @header = split/\n/,$vcf->{header_str};
        while ( my $a = $vcf->next_var() ) {
			my $ind = $a->{GT}->[0]->{_sample_id};
            my $simple_id = $a->{CHROM}."_".$a->{POS}."_".$a->{REF}."_".$a->{ALT};
            ## found in other ind
            if ($agg{$simple_id}) {
                $agg{$simple_id}->{INFO}->{more} .= "|$fn";
				my $exist = 1;
				my $index = 0;
				## Add gt to agg, if individual exist as remove dummy and add real
				foreach my $fm (@{ $agg{$simple_id}{GT}}) {	
					if ($fm->{_sample_id} eq $ind) {
						if ($fm->{GT} eq '0') {
							splice(@{ $agg{$simple_id}{GT}},$index, 1);
							push @{ $agg{$simple_id}{GT}},$a->{GT}->[0];
						}
					}
					$index++;
				}		
            }
            ## not found before
            else {
				add_info( $a, "more", $fn );
				add_gt($a,$ind);
				$agg{$simple_id} = $a;                    
            }
			
        }
    }
    return \%agg;
}

sub add_gt {
	my ($var,$ind) = @_;
	
	my @gts;
	foreach my $fm (@extra_individuals) {
		my %gt_add;
		if ($fm ne $ind) {
			foreach my $key ( @{ $var->{FORMAT} }) {
				$gt_add{$key} = 0;
			}
			$gt_add{_sample_id} = $fm;
			push @{ $var->{GT}},\%gt_add;
		}
	}

}

sub add_info {
    my( $var, $key, $val ) = @_;
    push( @{$var->{INFO_order}}, $key );
    $var->{INFO}->{$key} = $val;
}

sub vcfstr {
	my( $v, $sample_order ) = @_;
	
	my @all_info;
	my $tot_str = $v->{CHROM}."\t".$v->{POS}."\t".$v->{ID}."\t".$v->{REF}."\t".$v->{ALT}."\t".$v->{QUAL}."\t".$v->{FILTER}."\t";

	# Generate and print INFO field
	for my $info_key (@{$v->{INFO_order}}) {
		if($info_key eq "CSQ") {
			push @all_info, $info_key."=".$v->{_CSQstr};
		}
		else {
			push @all_info, $info_key."=".$v->{INFO}->{$info_key};
		}
	}
	$tot_str = $tot_str.join(";", @all_info)."\t";

	# Print FORMAT field
	$tot_str = $tot_str.join(":", @{$v->{FORMAT}})."\t";


	my %order;
	my $i=0;
	if( @$sample_order > 0 ) {
		$order{$_} = $i++ foreach @{$sample_order};
	}
	else {
		$order{$_->{_sample_id}} = $i++ foreach @{$v->{GT}};
	}
	# Print GT fields for all samples
	for my $gt ( sort {$order{$a->{_sample_id}} <=> $order{$b->{_sample_id}}} @{$v->{GT}}) {
		my @all_gt;
		for my $key ( @{$v->{FORMAT}} ) {
			push @all_gt, ( defined $gt->{$key} ? $gt->{$key} : "");
		}
		$tot_str = $tot_str.join(":", @all_gt)."\t";
	}
	$tot_str = $tot_str."\n";
	return $tot_str;
}

sub check_options {
    my %opt = %{ $_[0] };
    
    my @files;
    foreach my $opt_key ( sort keys %opt ) {

	if( $opt_key eq "vcfs" ) {
	    my @vcfs = split /,/, $opt{$opt_key};
	    foreach( @vcfs ){
		die "VCF does not exist: $_" unless -s $_;
		push @files, $_;
	    }
	}
    }
    return @files;
}