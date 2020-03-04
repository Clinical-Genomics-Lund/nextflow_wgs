#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename;
use lib dirname (__FILE__);
use vcf2;
use strict;
use Data::Dumper;



my $vcf = CMD::vcf2->new('file'=>$ARGV[0] );

print $vcf->{header_str};
while ( my $a = $vcf->next_var() ) {
    my %called = ();
    my @callset = split/-/,$a -> {INFO} -> {set};
    
    
    foreach (@callset){
        if ($_ =~ /manta/) {
            $called{manta} = 1;
        }
        if ($_ =~ /tiddit/) {
            $called{tiddit} = 1;
        }
        if ($_ =~ /cnvnator/) {
            $called{cnvnator} = 1;
        }
        if ($_ =~ /Intersection/) {
            $called{Intersection} = 1;
        }
    }
    my @set;
    foreach my $key (keys %called) {
        if (defined $called{$key}) {
            if ($key eq "Intersection") {
                push @set,"manta","tiddit","cnvnator";
            }
            else {
                push @set,$key;
            }
            
        }
    }

    $a->{INFO}->{set} = join('-',@set);
    my $print = vcfstr($a,[]);
	## extra tab at end crashes bedtools
	$print =~ s/\t$//;
    print $print;
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
