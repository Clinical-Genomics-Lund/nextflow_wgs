#! /usr/bin/perl -w
#use MongoDB;
use strict;
use Data::Dumper;
use Getopt::Long;
use JSON;
use List::MoreUtils qw(uniq);

my %opt = ();
GetOptions( \%opt, 'g=s', 'd=s', 'p=s', 'out=s', 'genome=s', 'antype=s', 'ped=s', 'assay=s', 'files=s', 'panelsdef=s' );

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
            'institute_owner' => 'onkogenetik'
        }
    },
    'wgs-hg38' => {
        'rankm' => '5.2.1_SNV',
        'svrankm' => '5.2_SV',
        'ph' => {
            'institute' => 'klingen_38',
            'institute_owner' => 'klingen_38'
        },
        'ahus' => {
            'institute' => 'ahus',
            'institute_owner' => 'ahus',
            'rstreshold' => -10000
        },
        'kit' => {
            'institute' => 'KIT',
            'institute_owner' => 'KIT'
        },
        'wgsvalid' => {
            'institute' => 'klingen_38',
            'institute_owner' => 'klingen-genlista'
        }
    },
    'wgs_dev' => {
        'rankm' => '5.2.1_SNV',
        'svrankm' => '5.2_SV',
        'ph' => {
            'institute' => 'klingen_38',
            'institute_owner' => 'klingen_38'
        },
        'ahus' => {
            'institute' => 'ahus',
            'institute_owner' => 'ahus'
        },
        'wgsvalid' => {
            'institute' => 'klingen_38',
            'institute_owner' => 'klingen-genlista'
        }
    },
    'clinicalwesv1-0' => {
        'rankm' => 'SNV-RM-v5.0',
        'svrankm' => 'SV-Panel-RM-v1.0',
        'hemato' => {
            'institute' => 'myeloid',
            'institute_owner' => 'myeloid',
        }
    }
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
    if ($a_a[1] ne 'false' && $a_a[1]) {
        $analysis = $a_a[1];
    }
    elsif ($opt{d} eq 'ahus') { ## beginning of stinking mess, deadline for fix 2021-03-01
        $analysis = 'ahus';
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
    elsif  ($tmp[0] eq "mtBAM") {
        $INFO{mtBAM}->{$tmp[1]} = $tmp[2];
    }
    elsif  ($tmp[0] eq "IMG") {
        $INFO{IMG}->{$tmp[1]} = $tmp[2];
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
    if ($INFO{mtBAM}{$pedline[1]}) {
        print OUT "    mt_bam: $INFO{mtBAM}{$pedline[1]}\n";
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
    if ($analysis ne 'kit') {
        print OUT "vcf_str: $INFO{STR}\n";
    }
}
if ($INFO{SMN}) {
    print OUT "smn_tsv: $INFO{SMN}\n";
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

## If IMGage is available
if ($INFO{IMG}) {

    print OUT "custom_images:\n";
    foreach my $img_type (keys %{ $INFO{IMG} }) {
        print OUT "  $img_type:\n";
        print OUT "    - title: $INFO{IMG}{$img_type}\n";
        print OUT "      description: Genome overview plot, UPD and ROH\n";
        # print OUT "      width: 4000\n";
        # print OUT "      height: 2300\n";
        print OUT "      path: $INFO{IMG}{$img_type}\n";
    }

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
if ($assays{$assay}{$analysis}{rstreshold}) {
    print OUT "rank_score_threshold: ". $assays{$assay}{$analysis}{rstreshold}."\n";
}
else {
    print OUT "rank_score_threshold: -1\n";
}
print OUT "human_genome_build: $genome\n";

close OUT;





sub get_genelist {
    my $institute = shift;
    
    my $file = $opt{panelsdef};
    my $data;
    my @ok_panels;
    open (JSON, $file);
    while (<JSON>) {
        $data = decode_json($_);
    }
   # print Dumper($data);
    foreach my $key (@{$data}) {
        if (ref $key->{institute} eq 'ARRAY') {
            foreach my $inst (@{ $key->{institute} }) {
                next if $key->{'display_name'} =~ /ERSATT|TEST|test|Test/;
                push @ok_panels,$key->{panel_name} if $inst eq $institute;
            }
        }
        elsif ($key->{institute} eq $institute) {
            next if $key->{'display_name'} =~ /ERSATT|TEST|test|Test/;
            push @ok_panels,$key->{panel_name};
        }
    }
    @ok_panels = uniq(@ok_panels);
    return \@ok_panels;
}
