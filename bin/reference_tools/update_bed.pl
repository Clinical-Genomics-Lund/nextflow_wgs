#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use File::Basename;
use lib dirname (__FILE__);
use vcf2 qw( parse_vcf );
use Getopt::Long;

my %opt = ();
my @incl_bed_files;
GetOptions( \%opt, 'old_clinvar=s', 'new_clinvar=s', 'build=s', 'clinvardate=s', 'incl_bed=s@' => \@incl_bed_files );

my @required_params = qw(old_clinvar new_clinvar build clinvardate);
my @missing_params;

foreach my $param (@required_params) {
    push @missing_params, $param unless defined $opt{$param};
}

if (@missing_params) {
    die "Error: Missing required parameters: " . join(", ", @missing_params) . "\nUsage: update_bed.pl --old_clinvar <filepath> --new_clinvar <filepath> --build <108> --clinvardate <20240701> [--incl_bed <filepath>]\n";
}

if (@incl_bed_files) {
    print "Additional included BED files: \n" . join("\n", @incl_bed_files) . "\n";
}

### Logic ###
# have a baseline bed that defines exonic regions, padded with 20(?) bp
# read clinvar vcf (likely pathogenic + pathogenic), convert it to bed, intersect. Regions not in baseline padded (5 bp?) and merge

#my $clinvar_vcf = "/fs1/resources/ref/hg38/annotation_dbs/bedupdatertest/new_chrom8.vcf"; 
#my $clinvar_vcf_old = "/fs1/resources/ref/hg38/annotation_dbs/bedupdatertest/old_chrom8.vcf";
#my $clinvar_vcf = "/fs1/resources/ref/hg38/annotation_dbs/clinvar38_latest.vcf.gz";
#my $clinvar_vcf_old = "/data/bnf/ref/annotations_dbs/clinvar38_20200106.vcf.gz";
my $clinvar_vcf = $opt{'new_clinvar'};
my $clinvar_vcf_old = $opt{'old_clinvar'};
# my $agilent = "agilient_hg38_nochr_noalt_1-3.bed";

### BASE BED ########
my $release = $opt{'build'}; # as argument
my $gtf = "Homo_sapiens.GRCh38.".$release.".gtf.gz"; 
my $gtf_request = "https://ftp.ensembl.org/pub/release-".$release."/gtf/homo_sapiens/$gtf";
my $base_bed = "exons_hg38_".$release.".bed";
## fetch $release from ensembl and get all exons for coding genes 

# FIXME: Uncomment
# get_base($gtf_request,$gtf,$base_bed,$release);

### BED TO ADD DATA TO ###
my $clinvar = $opt{'clinvardate'}; # as argument
my $final_bed = "exons_".$release."padded20bp_clinvar-".$clinvar."padded5bp."."bed";  # final name of bed, will concat to this

# if recreating same clinvarversion and ensembl build, delete old so it wont be concatenated
if(-e $final_bed) {
    unlink($final_bed);
}

## ADD DATA ##
# add_to_bed($final_bed, $base_bed, "EXONS-$release");

foreach my $incl_bed (@incl_bed_files) {
    print("Iterating with " . $incl_bed . "\n");
    my $suffix = $incl_bed;
    add_to_bed($final_bed, $incl_bed, $suffix);
}

# add_to_bed($final_bed,$agilent,"AGILENT-EXOME");

## clinvar variants ##
print("Clinvar 1\n");
my ($clinvar_new, $new_benign) = read_clinvar($clinvar_vcf);
print("Clinvar 2\n");
my ($clinvar_old, $old_benign) = read_clinvar($clinvar_vcf_old);
print("Clinvar 3\n");
my %old_benign = %$old_benign;
my %new_benign = %$new_benign;
my ($newtoadd, $oldtoremove) = compare_clinvar($clinvar_new, $clinvar_old, $final_bed);
my $clinvarlog = "clinvar_".$clinvar.".log";
my $clinvar_final_bed = "clinvar_".$clinvar.".bed";
if (-e $clinvarlog) { 
    unlink( $clinvarlog ); 
}
if (-e $clinvar_final_bed) { 
    unlink( $clinvar_final_bed ); 
}
clinvar_bed_and_info($newtoadd, "new", $clinvarlog);
clinvar_bed_and_info($oldtoremove, "old", $clinvarlog);
clinvar_bed_and_info($newtoadd, "bed", $clinvar_final_bed);
add_to_bed($final_bed, $clinvar_final_bed, ".");

print($final_bed, "\n");

sort_merge_output($final_bed);

print("Done\n");

## takes clinvar-vcf saves pathogenicity status, returns one hash with important variants
## and one with benign
sub read_clinvar {
    my $vcf = shift;
    my $vcf_hash = CMD::vcf2->new('file'=>$vcf);
    my %clinvar_variants;
    my %benign_clinvar_variants;
    my $clinvar_bed = "clinvar.bed";

    print("read_clinvar 1 " . $vcf . "\n");

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
    return \%clinvar_variants,\%benign_clinvar_variants;
}

## fetches $release ensembl gtf and returns bed with coding genes exons
sub get_base {
    my $gtf_req = shift;
    my $gtf = shift;
    my $base = shift;
    my $release = shift;
    system("wget $gtf_req -O tmp.gtf.gz");
    system("gunzip tmp.gtf.gz");
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
    my $new = shift;
    my $old = shift;
    my $final_bed = shift;
    my %new = %$new;
    my %old = %$old;
    my $new_bed = "clinvar_new.bed";
    my $old_bed = "clinvar_old.bed";
    open (NEW,'>',$new_bed);
    open (OLD,'>',$old_bed);
    my @incommon = ();
    foreach my $key (keys %new ) {
        push(@incommon, $key) if exists $old{$key};
        print NEW $new{$key}->{CHROM}."\t";
        print NEW ($new{$key}->{POS} - 5)."\t";
        print NEW ($new{$key}->{POS} + 5)."\t";
        my $info = clinvar_info($new{$key}, $key);
        # print("$info: " . $info . "\n");
        print NEW "$info\n";

    }
    my @new_added = ();
    foreach my $key (keys %new ) {
        push(@new_added, $key) unless exists $old{$key};
    }
    my @old_removed = ();
    foreach my $key (keys %old ) {
        unless (exists $new{$key}) {
            push(@old_removed, $key);
            print OLD $old{$key}->{CHROM}."\t";
            print OLD ($old{$key}->{POS} - 5)."\t";
            print OLD ($old{$key}->{POS} + 5)."\t";
            my $info = clinvar_info($old{$key}, $key);
            print OLD "$info\n";
        }

    }
    close NEW;
    close OLD;

    my @add2bed = push(@incommon, @new_added);
    my $newtoadd = intersect($new_bed, $final_bed);
    unlink($new_bed);
    my $oldtoremove = intersect($old_bed, $final_bed);
    unlink($old_bed);
    print "clinvar in common between versions :".scalar(@incommon)."\n";
    print "added new(unique targets)          :".scalar(@new_added)."(".scalar(@{$newtoadd}).")"."\n";
    print "removed old(unique targets)        :".scalar(@old_removed)."(".scalar(@{$oldtoremove}).")"."\n";
    return $newtoadd,$oldtoremove;
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
    if ( defined $fourth[3] ) {
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

    print("Inside clinvar_bed_and_info with $clinvar " . "\n");

    open (LOG, '>>', $clinvarlog);
    foreach my $clinvar_line (@{ $clinvar }) {
        chomp($clinvar_line);

        # print("With target: " . $target . "\n");

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

