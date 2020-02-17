#!/usr/bin/perl -w
use strict;


my $fn = $ARGV[0];

if( is_gzipped($fn) ) {
    open( VCF, "zcat $fn|");
}
else {
    open( VCF, $fn );
}

print 'track db="hg38" name="SVcalls" description="SVcalls"'."\n";
while(<VCF>) {
    next if /^#/;
    
    my @a = split /\t/;
    my( $chr, $start ) = ( $a[0], $a[1]-1 );
    my ($foo, $type, $bar) = ($a[7] =~ /(^|;)SVTYPE=(.*?)($|;)/);
    my ($foo2, $end, $bar2) = ($a[7] =~ /(^|;)END=(.*?)($|;)/);
    my( $id, $score, $foo3) = ($a[7] =~ /RankScore=(.*?):(.*?)($|;)/);

    next if $type eq "BND";
    next if $score < 0;
    $end = $start+1 unless $end;

    next if $end-$start < 100;
    print add_chr($chr)."\t".$start."\t".$end."\n";#."\t".$type."\n";

			       
}


sub add_chr {
    my $str = shift;
    unless($str =~ /^chr/) {
	return "chr$str";
    }
    return $str;
}

sub is_gzipped {
    my $fn = shift;

    my $file_str = `file -L $fn`;
    return 1 if $file_str =~ /gzip compressed/;
    return 0;
}
