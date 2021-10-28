#!/usr/bin/perl -w

use strict; 
use warnings;
use lib '.';

my @temp = ();

my $in;
my $out;
my %IDs;

#gene names
open ($in, '<', $ARGV[0]) or die "Input file $ARGV[0] not found!\n";
while(<$in>) {
	chomp;
	@temp = split("\t", $_);
	$IDs{$temp[0]} = $temp[1];	
}
close $in;

#values
open ($in, '<', $ARGV[1]) or die "Input file $ARGV[1] not found!\n";
while(<$in>) {
	chomp;
	@temp = split("\t", $_);
	
	#header
	if($temp[0] =~ m/pos/) {
		print $_;
		print "\n";
		next;
	}
	
	print $temp[0];
	$temp[0] =~ s/mt-//g;
	print "_".$IDs{$temp[0]};
	for(my $i=1; $i<@temp; $i++) {
		print "\t".$temp[$i];
	}
	print "\n";
	
}
close $in;



exit(0);









