#! /usr/bin/perl -w
use MongoDB;
use strict;
use Data::Dumper;
use Getopt::Long;
my %opt = ();
GetOptions( \%opt, 'g=s', 'd=s', 'p=s', 'out=s', 'genome=s', 'antype=s', 'ped=s', 'assay=s', 'files=s' );

### Define scout institute per assay and/or analysis. ###
my %assays = (
    'oncov1-0' => {
        'rankm' => 'SNV-RM-v5.0',
        'svrankm' => 'SV-Panel-RM-v1.0',
        'screening' => {
            'institute' => 'oncogen',
            'institute_owner' => 'onkogenetik',
        },
        'predictive' => {
            'institute' => 'oncogen',
            'institute_owner' => 'onkogenetik'
        },
        'PARP_inhib_normals' => {
            'institute' => 'oncogen',
            'institute_owner' => 'PARP_inhib_normals'
        }
    },
    'wgs-hg38' => {
        'rankm' => '5.1_SNV',
        'svrankm' => '5.1_SV',
        'ph' => {
            'institute' => 'klingen_38',
            'institute_owner' => 'klingen_38'
        }
    },
    'exome' => 'dummy'
);


### Analysis type ###
my $antype = "wgs";
if ($opt{antype}) { $antype = $opt{antype}; }

### Genome build ###
my $genome = "38";
if ($opt{genome}) { $genome = $opt{genome}; }

### Assay ###
my $assay = "wgs";
my $analysis = "";
if ($opt{assay}) { 
    my @a_a = split/,/,$opt{assay};
    $assay = $a_a[0];
    if ($a_a[1]) {
        $analysis = $a_a[1];
    }
    else { $analysis = 'ph';}
    
}

### Group ###
if (!defined $opt{g}) { print STDERR "need group name"; exit;}
my @g_c = split/,/,$opt{g};
my $group = $g_c[0];
my $clarity_id = $g_c[1];

### Open out, default $group.yaml ###
my $out = "$group.yaml";
if ($opt{out}) { $out = $opt{out}; }
open (OUT,'>',$out);

### Read ped, save individuals ####################
my $files = $opt{files};
open (INFO, $files) or die "Cannot open $files\n";
my %INFO;
my @bams;
while ( <INFO> ) {

    my @tmp = split/\s+/,$_;
    print STDERR $_,"\n";
    if ($tmp[0] eq "BAM") {
        $INFO{BAM}->{$tmp[1]} = $tmp[2];
        
    }
    elsif ($tmp[0] eq "TISSUE") {
        $INFO{TISSUE}->{$tmp[1]} = $tmp[2];
    }
    else {
        $INFO{$tmp[0]} = $tmp[1];
    }

}
close INFO;
####################################################

my $kit = "Intersected WGS"; ## placeholder, does not change for panels
my $diagnosis = $opt{d};

### Read ped, save individuals ####################
my $PED = $opt{ped};
open (PED, $PED) or die "Cannot open $PED\n";
my @ped;
while ( <PED> ) {
    push @ped, $_;
}
close PED;
####################################################

### PRINT YAML ####
print OUT "---\n";
my $institute = "klingen";
my $institute_owner = "klingen";
if ($opt{assay}) { 
    $institute = $assays{$assay}{$analysis}{institute};
    $institute_owner = $assays{$assay}{$analysis}{institute_owner};
}
### ASSAY DECIDE OWNER? ####
print OUT "owner: $institute_owner\n";
print OUT "family: '$group'\n";
print OUT "lims_id: '$clarity_id'\n";
print OUT "samples: \n";

### MATCH ped inidividuals with bams ###
foreach my $sample (@ped) {
    my @pedline = split/\t/,$sample;
    print OUT "  - analysis_type: $antype\n";
    print OUT "    sample_id: '$pedline[1]'\n";
    print OUT "    sample_name: '$pedline[1]'\n";
    print OUT "    mother: '$pedline[3]'\n";
    print OUT "    father: '$pedline[2]'\n";
    print OUT "    capture_kit: $kit\n";
    # unless ($INFO{TISSUE}{$pedline[1]} eq 'false') {
    #     print OUT "    tissue_type: '$INFO{TISSUE}{$pedline[1]}'\n";
    # }

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
    if ($INFO{BAM}{$pedline[1]}) {
        print OUT "    bam_path: $INFO{BAM}{$pedline[1]}\n";
    }
}
##########################################

### Print optional variables ###

## If trio and both SV and SNV calling has been done, SVc should be in INFO-file
## This contains both SV and SNV info
if ($INFO{SVc}) {
    my @tmp = split/,/,$INFO{SVc};
    print OUT "vcf_snv: $tmp[1]\n"; 
    print OUT "vcf_sv: $tmp[0]\n";
}
## If SNV single, check for SNV, if exist, check for SV
elsif ($INFO{SNV}) {
    print OUT "vcf_snv: $INFO{SNV}\n";
    if ($INFO{SV}) {
        print OUT "vcf_sv: $INFO{SV}\n";
    } 
}
## If only SV calling
elsif ($INFO{SV}) {
    print OUT "vcf_sv: $INFO{SV}\n";
}
else {
    print STDERR "need at least one VCF, SV/SNV"; exit;
}
if ($INFO{STR}) {
    print OUT "vcf_str: $INFO{STR}\n";
}


if ($INFO{MADDE}) {
    print OUT "madeline: $INFO{MADDE}\n";
}

if ($INFO{PEDDY}) {
    my @tmp = split/,/,$INFO{PEDDY};
    print OUT "peddy_ped: $tmp[1]\n";
    print OUT "peddy_check: $tmp[0]\n";
    print OUT "peddy_sex: $tmp[2]\n";
}



### Print gene panels ###
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
print OUT "rank_model_version: ".$assays{$assay}{rankm}."\n";
print OUT "sv_rank_model_version: ".$assays{$assay}{svrankm}."\n";
print OUT "rank_score_threshold: -1\n";
print OUT "human_genome_build: $genome\n";

close OUT;









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
