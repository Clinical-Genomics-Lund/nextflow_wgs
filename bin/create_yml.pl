#! /usr/bin/perl -w
use MongoDB;
use strict;
use Data::Dumper;
use Getopt::Long;
my %opt = ();
GetOptions( \%opt, 'b=s', 'g=s', 'd=s', 'p=s', 'snv=s', 'sv=s', 'str=s', 'out=s', 'genome=s', 'antype=s', 'ped=s', 'dir=s' );
my $vcf = $opt{snv};
my $vcf_str = $opt{str};
my $vcf_sv = $opt{sv};
my $out = $opt{out};
open (OUT,'>',$out);



my $antype = "wgs";
my $genome = "38";
my $kit = "Intersected WGS";

my $BAMS = $opt{b};
my @bams = split/,/, $BAMS;
my $group = $opt{g};

my $basedir = $opt{dir};
my $diagnosis = $opt{d};

my $PED = $opt{ped};
my $xml = $group.".ped.madeline.xml";
my $peddy_ped = $group.".peddy.ped";
my $ped_check = $group.".ped_check.csv";
my $sexcheck = $group.".sex_check.csv";

open (PED, $PED) or die "Cannot open $PED\n";
my @ped;
while ( <PED> ) {

    push @ped, $_;

}
close PED;
print OUT "---\n";
my $institute;
if ($diagnosis eq "validering") {
    print OUT "owner: wgsvalidering\n";
    $institute = "klingen";
}
else {
    print OUT "owner: klingen\n";
    $institute = "klingen";
}
print OUT "family: '$group'\n";
print OUT "samples: \n";

foreach my $sample (@ped) {
    my @pedline = split/\t/,$sample;
    print OUT "  - analysis_type: $antype\n";
    print OUT "    sample_id: '$pedline[1]'\n";
    print OUT "    sample_name: '$pedline[1]'\n";
    print OUT "    mother: '$pedline[3]'\n";
    print OUT "    father: '$pedline[2]'\n";
    print OUT "    capture_kit: $kit\n";
    if ($pedline[5] == 1) {
        print OUT "    phenotype: unaffected\n";
    }
    elsif ($pedline[5] == 2) {
        print OUT "    phenotype: affected\n";
    }
    else { print STDERR "not a valid phenotype!\n" }
    if ($pedline[4] == 1) {
        print OUT "    sex: male\n";
    }
    elsif ($pedline[4] == 2) {
        print OUT "    sex: female\n";
    }
    else { print STDERR "not a valid sex!\n" }
    my @match_bam = grep(/^$pedline[1]/, @bams);
    unless (scalar(@match_bam) == 1) { print STDERR "no matching bam"; exit; }
    print OUT "    bam_path: $basedir/bam/@match_bam\n";

}
print OUT "vcf_snv: $basedir/vcf/$vcf\n"; 
print OUT "vcf_str: $basedir/vcf/$vcf_str\n";
print OUT "vcf_sv: $basedir/vcf/$vcf_sv\n";

if (scalar(@bams) > 1 ) {
    print OUT "madeline: $basedir/ped/$xml\n";
}
print OUT "peddy_ped: $basedir/ped/$peddy_ped\n";
print OUT "peddy_check: $basedir/ped/$ped_check\n";
print OUT "peddy_sex: $basedir/ped/$sexcheck\n";
my $gene_panels = get_genelist($institute);
print OUT "gene_panels: [";
print OUT join ",", @$gene_panels;
print OUT "]\n";
if ($diagnosis eq "pediatrics") {
    print OUT "default_gene_panels: []\n";
}
else {
    my @panels = split /\+/, $diagnosis;
    my $panels_str = '"'. join('","', @panels). '"';
    print OUT "default_gene_panels: [$panels_str]\n";
}
print OUT "rank_model_version: 4.1\n";
print OUT "rank_score_threshold: -1\n";
print OUT "human_genome_build: $genome\n";


sub get_genelist {
    my $institute = shift;
    my $host = 'mongodb://cmdscout2.lund.skane.se/scout';
    if( $opt{p} ) {
        if( $ENV{$opt{p}} ) {
	        my $port = $ENV{$opt{p}};
	        $host = "mongodb://localhost:$port/loqusdb";
        }
        else {
	        die "No port envvar set for $opt{p}";
        }
    }
    my $client = MongoDB->connect($host);
    my $PANELS = $client->ns("scout.gene_panel");
    my $panels = $PANELS->find( {'institute'=>$institute} );

    my %ok_panels;
    my $counter = 0;
    while( my $panel = $panels->next ) {
        unless ( $panel->{'display_name'} =~ /ERSATT|TEST|test|Test/ ){
            my $tmp = $panel->{'panel_name'};
            $tmp = '"'.$tmp.'"';
	    $ok_panels{$tmp} = 1;
            $counter++;
        }
        #print $panel->{'display_name'},"\n";
    }


    #print $counter,"\n";
    my @ok_panels = sort {lc $a cmp lc $b } keys %ok_panels;
    return \@ok_panels;
}
