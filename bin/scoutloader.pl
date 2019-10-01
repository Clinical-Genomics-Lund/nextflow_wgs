#!/usr/bin/perl -w

use strict;
my $directory = '/fs1/results/cron/scout';

opendir (DIR, $directory) or die $!;

while (my $file = readdir(DIR)) {

    if ( $file =~ /\.yaml/) {
        my $fullpath = $directory."/".$file;
        scoutcommand($fullpath);
    }

}

closedir(DIR);


sub scoutcommand {
    my $yaml_file = shift;
    my $command = "ssh viktor\@cmdscout1.lund.skane.se 'scout load case $yaml_file'";
    my $log = '/fs1/results/cron/scout/scout_upload.log';
    my $datestring = localtime();
    open(LOG, '>>' , $log) or die $!;
    print LOG "$datestring :: $yaml_file was loaded using: $command\n\n";
    my $go = `$command`;
    close(LOG);
    unlink $yaml_file;
    
}