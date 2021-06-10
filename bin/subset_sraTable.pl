#!/usr/local/bin/perl

use strict;
use warnings;
use Getopt::Long;

my ($acc_list_file, $sra_table_file);

&GetOptions(
    "l=s" =>\$acc_list_file,
    "s=s" =>\$sra_table_file
    );

($acc_list_file & $sra_table_file) ||
    die "usage: $0 OPTIONS
where options are:\n -l <input list file containing the disiable accession number>\n -s <sra metadata table file download from the sra selector>\n";

my %names = ();
open(LIST, $acc_list_file) or die $!;
while(<LIST>){
	if(/^(\S+)/){
	    $names{$1} = 1;
	    print STDERR "$1\n";
	}


    }
    close(LIST);


open(SRA, "<$sra_table_file") or die "cannot open $sra_table_file for reading: $!";
while(<SRA>){

    if(my ($acc) =  /^(\S+?),/){

	if(exists $names{$acc}){
	    print $_;
	}

    }

}

close(SRA);
