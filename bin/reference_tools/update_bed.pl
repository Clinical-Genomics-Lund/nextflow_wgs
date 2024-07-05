#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use File::Basename;
use lib dirname (__FILE__);
use vcf2 qw( parse_vcf );
use Getopt::Long;

my %opt = (
    skip_download = 0
);
my @incl_bed_files;
GetOptions( \%opt, 'old=s', 'new=s', 'build=s', 'clinvardate=s', 'incl_bed=s@' => \@incl_bed_files, 'skip_download'  );

my @required_params = qw(old new build clinvardate);
my @missing_params;

foreach my $param (@required_params) {
    push @missing_params, $param unless defined $opt{$param};
}

if (@missing_params) {
    die "Error: Missing required parameters: " . join(", ", @missing_params) . "\nUsage: update_bed.pl --old <filepath> --new <filepath> --build <108> --clinvardate <20240701> [--incl_bed <filepath>] [--skip_download]\n";
}

if (@incl_bed_files) {
    print "Additional included BED files: " . join(", ", @incl_bed_files) . "\n";
}

### Logic ###
# have a baseline bed that defines exonic regions, padded with 20(?) bp
# read clinvar vcf (likely pathogenic + pathogenic), convert it to bed, intersect. Regions not in baseline padded (5 bp?) and merge

#my $clinvar_vcf = "/fs1/resources/ref/hg38/annotation_dbs/bedupdatertest/new_chrom8.vcf"; 
#my $clinvar_vcf_old = "/fs1/resources/ref/hg38/annotation_dbs/bedupdatertest/old_chrom8.vcf";
#my $clinvar_vcf = "/fs1/resources/ref/hg38/annotation_dbs/clinvar38_latest.vcf.gz";
#my $clinvar_vcf_old = "/data/bnf/ref/annotations_dbs/clinvar38_20200106.vcf.gz";
my $clinvar_vcf = $opt{'new'};
my $clinvar_vcf_old = $opt{'old'};
# my $agilent = "agilient_hg38_nochr_noalt_1-3.bed";

### BASE BED ########
my $release = $opt{'build'}; # as argument
my $gtf = "Homo_sapiens.GRCh38.".$release.".gtf.gz"; 
my $gtf_request = "https://ftp.ensembl.org/pub/release-".$release."/gtf/homo_sapiens/$gtf";
my $base_bed = "exons_hg38_".$release.".bed";

## fetch $release from ensembl and get all exons for coding genes 
get_base($gtf_request, $gtf, $base_bed, $release, $opt{'skip_download'});

### BED TO ADD DATA TO ###
my $clinvar = $opt{'clinvardate'};
# final name of bed, will concat to this
my $final_bed_fp = "exons_".$release."padded20bp_clinvar-".$clinvar."padded5bp."."bed";

# if recreating same clinvarversion and ensembl build, delete old so it won't be concatenated
if(-e $final_bed_fp) {
    unlink($final_bed_fp);
}

## ADD DATA ##
add_to_bed($final_bed_fp, $base_bed, "EXONS-$release");

foreach my $incl_bed (@incl_bed_files) {
    my $suffix = $incl_bed;
    add_to_bed($final_bed_fp, $incl_bed, $suffix);
}

# add_to_bed($final_bed_fp, $agilent,"AGILENT-EXOME");

## clinvar variants ##
my ($clinvar_new, $new_benign) = read_clinvar($clinvar_vcf);
my ($clinvar_old, $old_benign) = read_clinvar($clinvar_vcf_old);

my $num_clinvar_new = scalar keys %$clinvar_new;
my $num_clinvar_old = scalar keys %$clinvar_old;
my $num_benign_new = scalar keys %$new_benign;
my $num_benign_old = scalar keys %$old_benign;

print("$num_clinvar_new $num_clinvar_old $num_benign_new $num_benign_old\n");

my %old_benign = %$old_benign;
my %new_benign = %$new_benign;
my ($newtoadd, $oldtoremove) = compare_clinvar($clinvar_new, $clinvar_old, $final_bed_fp);
my $clinvarlog = "clinvar_".$clinvar.".log";
my $clinvar_final_bed_fp = "clinvar_".$clinvar.".bed";
if (-e $clinvarlog) { 
    unlink( $clinvarlog ); 
}
if (-e $clinvar_final_bed_fp) { 
    unlink( $clinvar_final_bed_fp ); 
}
clinvar_bed_and_info($newtoadd, "new", $clinvarlog, $new_benign);
clinvar_bed_and_info($oldtoremove, "old", $clinvarlog, $new_benign);
clinvar_bed_and_info($newtoadd, "bed", $clinvar_final_bed_fp, $new_benign);
add_to_bed($final_bed_fp, $clinvar_final_bed_fp, ".");

sort_merge_output($final_bed_fp);

## takes clinvar-vcf saves pathogenicity status, returns one hash with important variants
## and one with benign
sub read_clinvar {
    my $vcf = shift;
    my $vcf_hash = CMD::vcf2->new('file'=>$vcf);
    my %clinvar_variants;
    my %benign_clinvar_variants;
    my $clinvar_bed = "clinvar.bed";

    ## save all variants that is marked as Pathogenic in some way
    while ( my $var = $vcf_hash->next_var() ) {

        my $sig = ""; my $confidence = "";
        my $haplo = "";
        if ($var->{INFO}->{CLNSIG}) {
            $sig = $var->{INFO}->{CLNSIG};
        }
        if ($var->{INFO}->{CLNSIGINCL}) {
            $haplo = $var->{INFO}->{CLNSIGINCL};
        }
        if ($var->{INFO}->{CLNSIGCONF}) {
            $confidence = $var->{INFO}->{CLNSIGCONF};
        }        
        next if ($sig eq "" && $haplo);
        my $keep = 0;
        my $reason = "";
        if ( $sig =~ /Pathogenic/ || $sig =~ /Likely_pathogenic/ ) {
            $keep = 1;
            $reason = $sig;
        }
        elsif ( $sig =~ /Conflicting_interpretations_of_pathogenicity/ ) {
            if ( $confidence =~ /Pathogenic/ || $confidence =~ /Likely_pathogenic/ ) {
                $keep = 1;
                $reason = $confidence;
            }
        }
        elsif ( $haplo =~ /athogenic/ ) {
            $keep = 1;
            $reason = $haplo;
        }


        if ($keep) {
            next if ($var->{CHROM} eq "MT");
            $clinvar_variants{$var->{CHROM}.":".$var->{POS}."_".$var->{REF}."_".$var->{ALT}}{INFO} = $var->{INFO};
            $clinvar_variants{$var->{CHROM}.":".$var->{POS}."_".$var->{REF}."_".$var->{ALT}}{CHROM} = $var->{CHROM};
            $clinvar_variants{$var->{CHROM}.":".$var->{POS}."_".$var->{REF}."_".$var->{ALT}}{POS} = $var->{POS};
            $clinvar_variants{$var->{CHROM}.":".$var->{POS}."_".$var->{REF}."_".$var->{ALT}}{REASON} = $reason;

        }
        else {
            $benign_clinvar_variants{$var->{CHROM}.":".$var->{POS}."_".$var->{REF}."_".$var->{ALT}}{INFO} = $var->{INFO};
            $benign_clinvar_variants{$var->{CHROM}.":".$var->{POS}."_".$var->{REF}."_".$var->{ALT}}{REASON} = $reason;
        }
    }    
    return \%clinvar_variants, \%benign_clinvar_variants;
}

## fetches $release ensembl gtf and returns bed with coding genes exons
sub get_base {
    my $gtf_req = shift;
    my $gtf = shift;
    my $base = shift;
    my $release = shift;
    my $skip_download = shift;
    unless ($skip_download) {
        system("wget $gtf_req -O tmp.gtf.gz");
        system("gunzip tmp.gtf.gz");
    }
    open(BASE,'>',$base);
    open(GTF, "tmp.gtf");
    my $keep = 0;
    while(<GTF>) {
        chomp;
        next if /^#/;
        my @g =split /\t/;
        next if length($g[0]) >2;
        if( $g[2] eq "transcript" ) {
            $keep = 0;
            $keep = 1 if $g[8] =~ /transcript_biotype "protein_coding"/;
        }
        if( $g[2] eq "exon" and $keep ) {
            print BASE $g[0]."\t".($g[3]-20)."\t".($g[4]+20)."\n";
        }
    }
    ## All of mitochondria ##
    close GTF;
    close BASE;    
    
}
## concats bed files and adds 4th column describing why it was added
sub add_to_bed {
    my $bed = shift;
    my $bed2add = shift;
    my $forthcolumn = shift;
    open(BED2ADD, $bed2add);
    open(BED, '>>', $bed);
    while (<BED2ADD>) {
        chomp;
        my @line = split('\t');
    next if $line[0] eq "M";
        if ($line[3]) {
            print BED $line[0]."\t".$line[1]."\t".$line[2]."\t".$line[3]."\n";    
        }
        else {
            print BED $line[0]."\t".$line[1]."\t".$line[2]."\t".$forthcolumn."\n";
        }
        
    }    
    close BED;
    close BED2ADD;
}
## merges bed-file targets and collapses 4th column
sub sort_merge_output {
    my $bed = shift;
    system( "bedtools sort -i $bed > tmp.sort.bed" );
    system( "bedtools merge -i tmp.sort.bed -c 4 -o collapse > $bed" );
    unlink( "tmp.sort.bed");
}
## finds calculates which variants to add to new bed, and which to remove
## depending on variant status. Take clinvar-hashes (old and new) returns
## list of variants to create bed
sub compare_clinvar {
    my $new_clinvar = shift;
    my $old_clinvar = shift;
    my $final_bed_fp = shift;
    my %new_clinvar = %$new_clinvar;
    my %old_clinvar = %$old_clinvar;

    my $new_clinvar_count = scalar keys %new_clinvar;
    my $old_clinvar_count = scalar keys %old_clinvar;

    my $new_bed_fp = "clinvar_new.bed";
    my $old_bed_fp = "clinvar_old.bed";
    open (NEW, '>', $new_bed_fp);
    open (OLD, '>', $old_bed_fp);
    my @clinvar_in_common = ();
    foreach my $key (keys %new_clinvar ) {
        push(@clinvar_in_common, $key) if exists $old_clinvar{$key};
        print NEW $new_clinvar{$key}->{CHROM}."\t";
        print NEW ($new_clinvar{$key}->{POS} - 5)."\t";
        print NEW ($new_clinvar{$key}->{POS} + 5)."\t";
        my $info = clinvar_info($new_clinvar{$key}, $key);
        print NEW "$info\n";
    }
    my $clinvar_in_common_nbr = scalar @clinvar_in_common;

    my @new_added = ();
    foreach my $key (keys %new_clinvar ) {
        push(@new_added, $key) unless exists $old_clinvar{$key};
    }
    my @old_removed = ();
    foreach my $key (keys %old_clinvar ) {
        unless (exists $new_clinvar{$key}) {
            push(@old_removed, $key);
            print OLD $old_clinvar{$key}->{CHROM}."\t";
            print OLD ($old_clinvar{$key}->{POS} - 5)."\t";
            print OLD ($old_clinvar{$key}->{POS} + 5)."\t";
            my $info = clinvar_info($old_clinvar{$key}, $key);
            print OLD "$info\n";
        }

    }
    close NEW;
    close OLD;

    # FIXME: Does this make sense? Looks like it always will add back the new
    # to the already existing ones
    # Commenting until further clarification ...
    # push(@clinvar_in_common, @new_added);

    my $new_to_add = intersect($new_bed_fp, $final_bed_fp);
    unlink($new_bed_fp);
    my $old_to_remove = intersect($old_bed_fp, $final_bed_fp);
    unlink($old_bed_fp);
    print "ClinVar in common between versions :".scalar(@clinvar_in_common)."\n";
    print "Added new(unique targets)          :".scalar(@new_added)."(".scalar(@{$new_to_add}).")"."\n";
    print "Removed old(unique targets)        :".scalar(@old_removed)."(".scalar(@{$old_to_remove}).")"."\n";
    return $new_to_add, $old_to_remove;
}
## finds unique targets, important to keep track of what variants gets added
## and what variants get removed
sub intersect {
    my $clinvar = shift;
    my $bed = shift;
    my @notinbed = `bedtools intersect -a $clinvar -b $bed -v`;
    return \@notinbed;
}
## concats INFO-field into 4th column for clinvarbed
sub clinvar_info {
    my $info = shift;
    my $var = shift;
    my %info = %$info;
    my @fourth;

    push (@fourth, $var);
    push (@fourth, $info->{REASON});
    push (@fourth, $info->{INFO}->{CLNACC});

    my $clndn = "Undefined";
    if (defined $fourth[3]) {
        $clndn = $info->{INFO}->{CLNDN};
    }

    push (@fourth, $clndn);

    my $fourth = join('~', @fourth);
    return $fourth;
}
## creates a log for added and removed clinvar variants as well
## as creates the bed that gets added to the final big bed-file
sub clinvar_bed_and_info {
    my $clinvar = shift;
    my $new_or_old = shift;
    my $clinvarlog = shift;
    my $new_benign = shift;

    open (LOG, '>>', $clinvarlog);
    foreach my $clinvar_line (@{ $clinvar }) {
        chomp($clinvar_line);


        my @clinvar_fields = split('\t', $clinvar_line);
        my @clinvarreason = split("~", $clinvar_fields[3]);

        while (@clinvarreason < 3) {
            push @clinvarreason, "[missing value]" 
        }

        my @pos = split('_' ,$clinvarreason[0]);
        if ($new_or_old eq "old") {
            my $innew = "MISSING";
            if ($new_benign{$clinvarreason[0]}->{INFO}->{CLNSIG}) {
                $innew = $new_benign{$clinvarreason[0]}->{INFO}->{CLNSIG};
            }
            print LOG "REMOVED:".$pos[0].":".$clinvarreason[2].":".$clinvarreason[1]."=>".$innew."\n";
        }
        elsif ($new_or_old eq "new") {
            print LOG "ADDED:".$pos[0].":".$clinvarreason[2].":".$clinvarreason[1]."\n";
        }
        else {
            print LOG $clinvar_fields[0]."\t".$clinvar_fields[1]."\t".$clinvar_fields[2]."\t"."CLINVAR-".$clinvarreason[1]."\n";
        }
        
    }
    close LOG;
}

