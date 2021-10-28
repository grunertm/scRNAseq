#!/usr/bin/perl -w

use strict; 
use warnings;
use lib '.';

my @temp = ();
my $in;
my $out;
my %IDs;
my %IDs_redundant;

#gene names
open ($in, '<', $ARGV[0]) or die "Input file $ARGV[0] not found!\n";
while(<$in>) {
	chomp;
	@temp = split("\t", $_);
		
	if(exists($IDs_redundant{$temp[1]})) {
		$IDs_redundant{$temp[1]} += 1;
	} else {
		$IDs_redundant{$temp[1]} = 1;
	}
	
	#Key:ENSEMBL ID and value:gene name
	$IDs{$temp[0]} = $temp[1];	
}
close $in;

#gene names
open ($in, '<', $ARGV[0]) or die "Input file $ARGV[0] not found!\n";
while(<$in>) {
	chomp;
	@temp = split("\t", $_);
		
	if($IDs_redundant{$temp[1]}>1) {
		$IDs{$temp[0]} = $temp[1]."_".$temp[0];
	} else {
		#Key:ENSEMBL ID and value:gene name
		$IDs{$temp[0]} = $temp[1];	
	}
}
close $in;


#values
open ($in, '<', $ARGV[1]) or die "Input file $ARGV[1] not found!\n";
while(<$in>) {
	chomp;
	@temp = split("\t", $_);
	
	#header
	if($temp[0] =~ m/ENSEMBL/) {
		$temp[0] = "Gene";
		next;
	}
	
	if(exists($IDs{$temp[0]})) {
		$temp[0] = $IDs{$temp[0]};
	}

	print $temp[0]."\t";	
	for(my $i=1; $i<@temp; $i++) {
		print $temp[$i]."\t";
	}
	print "\n";
	
}
close $in;


exit(0);









