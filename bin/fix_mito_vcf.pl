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
my @variants;
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
        push @variants,$print;     
    }
    else {
        $a->{CHROM} = "M";      ## instead of chrM, consistent with fasta
        my $print = vcfstr($a,[]); 
        push @variants,$print;
    }

}
#print join("\n",@variants);
print_and_mergedel(@variants);

sub print_and_mergedel {
    my @variants = @_;
    my $lastpos;
    my $withindel = 1;
    my $stretch = 0;
    my $newref = ""; my $newalt = "";
    for (my $i = 0; $i <= scalar(@variants)-1; $i++) {
        my $variant = $variants[$i];
        for (my $j = 1; $j <= scalar(@variants)-1; $j++) {
            my @var = split("\t",$variants[$i]); my @var2 = split("\t",$variants[$j]);
            my $pos = $var[1]; my $pos2 = $var2[1];
            my $ref = $var[3]; my $alt = $var[4]; my $ref2 = $var2[3]; my $alt2 = $var2[4];
            ## if next position is diff == 1, next time 2 and so on
            if ( abs($pos - $pos2) == $withindel) {
                ## if it is a deletion
                if  (length($ref) - length($alt) == 1 and length($ref2) - length($alt2) == 1) {
                    $withindel++; #print "$withindel \t $i \t $j $ref $alt $ref2 $alt2\n"; 
                    $newref = $newref.substr($ref2,1,1  ); $newalt = $alt;
                    
                    $variant = "M\t".$pos."\t".".\t".$ref.$newref."\t".$newalt."\t".join("\t",@var[5..$#var]);
                    next;
                }
                ## if no deletion end aggregation and jump forward $i to the variant after stetch of del
                else {
                    $i = $i + $withindel -1;
                    $withindel = 1; 
                }
            }
            ## if no deletion end aggregation and jump forward $i to the variant after stetch of del
            else {
                $i = $i + $withindel -1;
                $newref = ""; $newalt = "";
                $withindel = 1; 
            }
        }
        print $variant;
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