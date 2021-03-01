#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename;
use lib dirname (__FILE__);
use vcf2;
use strict;
use Data::Dumper;


my $vcf = CMD::vcf2->new('file'=>$ARGV[0] );
my $fasta = $ARGV[1];
my @nondel;
my @dels;

my @header = split/\n/,$vcf->{header_str};

## Add the relatives to VCF-header
foreach my $line (@header) {
	if ($line =~ /^##contig/) {
		print "##contig=<ID=M,length=16569>\n";
	}
    else {
        print $line,"\n";
    }

}


while ( my $a = $vcf->next_var() ) {
    my $REF = $a->{REF};
    my $ALT = $a->{ALT};
    ## find deletions, save positions to bed-file
    if ($ALT =~ /\*/) {
        my $pos2 = $a->{POS}-2;
        my $pos1 = $a->{POS}-1;
        open (BED, '>', "del.bed");
        print BED "M\t".$pos2."\t".$pos1."\t".$a->{POS}.":".$REF.":*\n";
        close BED;
        my ($pos,$ref,$prebase) = bedtools();
        if ($pos != $a->{POS} or $ref ne $a->{REF}) { print STDERR "missmatching from bedtools..exiting"; exit;}  ## sanity check, bedtools and perl desync?
        $a->{CHROM} = "M";      ## instead of chrM, consistent with fasta
        $a->{POS} = $a->{POS}-1; ## move pos one base left
        $a->{REF} = $prebase.$a->{REF}; ## concat ref with prebase, ref is deleted
        $a->{ALT} = $prebase; ## alt now becomes prebase only, as ref is deleted
        my $print = vcfstr($a,[]);  ## print changed vcf
        print $print;
    }
    else {
        $a->{CHROM} = "M";      ## instead of chrM, consistent with fasta
        my $print = vcfstr($a,[]); 
        print $print;
    }

}


sub bedtools {
    my $cmd = "bedtools getfasta -fi $fasta -name -bed del.bed -tab -fo prebase_del";
    #print "Running..".$cmd."\n";
    my $results = system($cmd);
    if ($results == 0) { ## if bedtools ran successfully
        open (FASTABED, "prebase_del") or die $!;
        my $prebase;
        my $ref;
        my $pos;
        while (<FASTABED>) {
            chomp;
            my @row = split/\t/;
            $prebase = $row[1];
            @row = split/:/,$row[0]; 
            $pos = $row[0]; ## position to match
            $ref = $row[1]; ## sanity check
        }
        #$prebase =~ s/\r\n//;
        return $pos,$ref,$prebase;
    }

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