#! /usr/bin/perl -w
use MongoDB;
use strict;
use Data::Dumper;
my $antype = "wes";
my $genome = "37";
my $kit = "Agilent_SureSelectCRE.V2";

my $BAMS = $ARGV[0];
my @bams = split/,/, $BAMS;
my $PED = $ARGV[1];
my $group = $ARGV[2];
my $vcf = $ARGV[3];
my $xml = $ARGV[4];
my $peddy_ped = $ARGV[5];
my $ped_check = $ARGV[6];
my $sexcheck = $ARGV[7];
my $basedir = $ARGV[8];
my $diagnosis = $ARGV[9];
open (PED, $PED) or die "Cannot open $PED\n";
my @ped;
while ( <PED> ) {

    push @ped, $_;

}
close PED;
print "---\n";
my $institute;
if (scalar(@bams) > 1 ) {
    print "owner: klingen\n";
    $institute = "klingen";
}
else {
    print "owner: klingen_genlista\n";
    $institute = "klingen";
}
print "family: '$group'\n";
print "samples: \n";
my $count = 0;
foreach my $sample (@ped) {
    my @pedline = split/\t/,$sample;
    print "  - analysis_type: $antype\n";
    print "    sample_id: '$pedline[1]'\n";
    print "    sample_name: '$pedline[1]'\n";
    print "    mother: '$pedline[3]'\n";
    print "    father: '$pedline[2]'\n";
    print "    capture_kit: $kit\n";
    if ($pedline[5] == 1) {
        print "  phenotype: unaffected\n";
    }
    elsif ($pedline[5] == 2) {
        print "  phenotype: affected\n";
    }
    else { print STDERR "not a valid phenotype!\n" }
    if ($pedline[4] == 1) {
        print "  sex: male\n";
    }
    elsif ($pedline[4] == 2) {
        print "  sex: female\n";
    }
    else { print STDERR "not a valid sex!\n" }
    my @match_bam = grep(/^$pedline[1]/, @bams);
    unless (scalar(@match_bam) == 1) { print STDERR "no matching bam"; exit; }
    print "  bam_path: $basedir/bam/wgs/@match_bam\n";
    $count++;
}
print "vcf_snv: $basedir/vcf/wgs/$vcf\n"; ##obs skapa en version för exome specifikt
if (scalar(@bams) > 1 ) {
    print "madeline: $basedir/ped/wgs/$xml\n";
}
print "peddy_ped: $basedir/ped/wgs/$peddy_ped\n";
print "peddy_check: $basedir/ped/wgs/$ped_check\n";
print "peddy_sex: $basedir/ped/wgs/$sexcheck\n";
my $gene_panels = get_genelist($institute);
print "gene_panels: [";
print join ",", @$gene_panels;
print "]\n";
if ($diagnosis eq "pedriatrics") {
    print "default_gene_panels: []\n";
}
else {
    print "default_gene_panels: [$diagnosis]\n";
}
print "rank_model_version: 3.0\n";
print "rank_score_threshold: 0\n";
print "human_genome_build: $genome\n";


sub get_genelist {
    my $institute = shift;
    my $host = 'mongodb://cmdscout2.lund.skane.se/scout';
    if( $ARGV[10] ) {
        if( $ENV{$ARGV[10]} ) {
	        my $port = $ENV{$ARGV[10]};
	        $host = "mongodb://localhost:$port/loqusdb";
        }
        else {
	        die "No port envvar set for $ARGV[1]";
        }
    }
    my $client = MongoDB->connect($host);
    my $PANELS = $client->ns("scout.gene_panel");
    my $panels = $PANELS->find( {'institute'=>$institute} );

    my @ok_panels;
    my $counter = 0;
    while( my $panel = $panels->next ) {
        unless ( $panel->{'display_name'} =~ /ERSATT|TEST|test|Test/ ){
            my $tmp = $panel->{'panel_name'};
            $tmp = '"'.$tmp.'"';
            push @ok_panels, $tmp;
            $counter++;
        }
        #print $panel->{'display_name'},"\n";
    }


    #print $counter,"\n";
    return \@ok_panels;
}