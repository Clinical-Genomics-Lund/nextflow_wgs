#!/usr/bin/perl -w

use strict;

use Getopt::Long;
my %opt = ();
GetOptions( \%opt, 'i=s', 'o=s' );
my $infile = $opt{i};
my $outfile = $opt{o};


open (YAML, $infile) or die "Cannot open $infile\n";
open (OUT, '>' ,$outfile) or die "Cannot open $outfile\n";
my $snv;
my @bam;
my $ped;
my $id;
my $group;
print OUT "id,group,assay,bam,ped,snv\n";
while (<YAML>) {

	if ($_=~ /^vcf_snv:\s?(.+)/) {
		$snv = $1;
	}
	if ($_=~ /^\s+bam_path:\s?(.+)/) {
		push @bam,$1;
	}
	if ($_=~ /^peddy_ped:\s?(.+)/) {
		$ped = $1;
		$ped =~ s/peddy\.//;
		my @ped = split'/',$ped;
		$group = $ped[-1];
		$group =~ s/\.ped//;
	}
}

foreach my $bam (@bam) {
	my @id = split'/',$bam;
	$id = $id[-1];
	$id =~ s/_merged_dedup\.bam//;
	print OUT "$id,$group-sv,wgs_cnv38,$bam,$ped,$snv\n";
}
