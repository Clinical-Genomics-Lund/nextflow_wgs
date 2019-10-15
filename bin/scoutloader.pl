#!/usr/bin/perl -w
use Backticks;
use strict;
my $directory = '/fs1/results/cron/scout';

opendir (DIR, $directory) or die $!;

while (my $file = readdir(DIR)) {

    if ( $file =~ /\.yaml$/) {
        scoutcommand($directory,$file);
    }

}

closedir(DIR);


sub scoutcommand {
    my ($directory, $file) = @_;

    my $yaml_file = $directory."/".$file;

    my $out_file = "~/scout_upload.stdout.".int rand 100000000;
    rename($yaml_file, $yaml_file.".inprogress");
    my $command = "ssh viktor\@cmdscout1.lund.skane.se 'scout load case $yaml_file.inprogress' &>> $out_file";
    my $log = '/fs1/results/cron/scout/scout_upload.log';
    my $errlog = '/fs1/results/cron/scout/scout_upload.errlog';
    my $datestring = localtime();
    open(LOG, '>>' , $log) or die $!;
    
    my $results = `$command`;
    my $status = $results->success;

    if ( $status ) {
        print LOG "$datestring :: $yaml_file was loaded using: $command\n";
        unlink $yaml_file.".inprogress";
    }
    else {
        open(ERRLOG, '>>' , $errlog) or die $!;
	print ERRLOG "$datestring :: $yaml_file: could not be loaded.\n";
        close(ERRLOG);
	rename($yaml_file.".inprogress", $yaml_file.".failed");
    }
    close(LOG);
    
    
}
