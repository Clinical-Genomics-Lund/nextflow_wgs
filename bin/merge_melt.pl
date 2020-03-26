#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename;
use lib dirname (__FILE__);
use vcf2;
use strict;
use Data::Dumper;

my $vcfheader = $ARGV[0];
my $id = $ARGV[1];


my $alu = "ALU.final_comp.vcf";
my $line1 = "LINE1.final_comp.vcf";
my $sva = "SVA.final_comp.vcf";

my %qcval = (
    'lc' => 'Low Complex Region',
    's25' => 'Greater than 25.0% of samples do not have data',
    'hDP' => '"More than the expected number of discordant pairs at this site are also split',
    'PASS' => 'PASS',
    'rSD' => "Ratio of LP to RP is greater than 2.0 standard deviations"
);


my @files;
push @files,$alu,$line1,$sva;

my @non_empty;
foreach my $file (@files) {
    my $size = -s $file;
    if ($size < 1) { next; }

    my $vcf = CMD::vcf2->new('file'=>$file);
    my $out = $file."mod";
    open (OUT, '>' ,$out);
    my @header = split/\n/,$vcf->{header_str};
    my $count = 1;
    foreach my $header (@header) {
        $header =~ s/\'\\\|\'/pipe-sign/g;
        
        if ($header =~ /^##INFO/ && $count == 1) {
            print OUT '##INFO=<ID=END,Number=.,Type=Integer,Description="END position set to start position for insertions">'."\n";
            $count++;
        }
        if ($header =~ /^#CHROM/) {
            print OUT "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$id\n";
        }
        else {
            print OUT $header."\n";
        }

    }

    while ( my $a = $vcf->next_var() ) {

        my @str = split/\t/,$a->{vcf_str};
        my $start = $a -> {POS};

        ## Filter variants in introns and invalid positions
        if ($a->{INFO}->{INTERNAL} =~ /INTRONIC|null/) { next; }
        ## Filter variants due to low quality for sample
        if ($a->{FILTER} =~ /ac0/) { next; }

        ## Print first 7 columns of vcf
        print OUT join("\t",@str[0..6])."\t";
        
        ## ADD END to info and change SVTYPE to INS
        my @INFO = split/;/,$str[7];
        splice @INFO, 5, 0, "END=$start";
        my $ref = \@INFO;
        ${$ref}[3] = "SVTYPE=INS";

        my @alt = split/:/,$a->{ALT};
        my $alt = $alt[2];
        $alt =~ s/>//g;
        my $event = $a->{INFO}->{INTERNAL};
        $event =~ s/\,/_/g;


        push @$ref,"SCOUT_CUSTOM=Type|".$alt.","."Filter|".$qcval{$a->{FILTER}}.","."Consequence|".$event;
        print OUT join(';',@$ref)."\t";
        print OUT join("\t",@str[8..$#str]);
        print OUT "\n";

    }

    close OUT;
    push @non_empty,$out;

}

## If no variants in either of three VCFs, print empty VCF with only header
if (scalar(@non_empty) == 0) { 
    my $command = `cat $vcfheader > $id\.melt.merged.vcf`;
    my $sed = `sed -i 's/FORMAT/FORMAT\t$id/' $id\.melt.merged.vcf`;
}
## Concat all with values
else {
    my $vcfs = join(' ',@non_empty);
    my $command = `vcf-concat $vcfs | vcf-sort -c > $id\.melt.merged.vcf`;
}
