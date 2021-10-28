#!/usr/bin/perl -w

use strict; 
use warnings;
use lib '.';

my @temp = ();
my @temp2 = ();
my $in;
my $in2;

my %matrix;
my %genes;
my %cells;

#input files
open ($in, '<', $ARGV[0]) or die "Input file $ARGV[0] not found!\n";
while(<$in>) {
	chomp;
	@temp = split(" ", $_);
	
	my @sample = split(/\./, $temp[0]);
	
	print STDERR $sample[0]."\n";
	
	open ($in2, '<', $temp[0]) or die "Input file $temp[0] not found!\n";
	while(<$in2>) {
		chomp;
		@temp2 = split("\t", $_);
		
		if($temp[0] =~ m/gene\_id/) {
			next;
		}
		
		#store TPM
		$matrix{$temp2[0]}{$sample[0]} = $temp2[4];
		
		$genes{$temp2[0]} = 1;
		$cells{$sample[0]} = 1;

	}
	close $in2;
	
}
close $in;


#print header
print "ENSEMBL_Gene_ID";
for my $cell ( sort keys %cells ) {
	print "\t".$cell;
}
print "\n";


for my $gene ( sort keys %genes ) {
	print $gene;
	for my $cell ( sort keys %cells ) {
		
		if(exists($matrix{$gene}{$cell})) {
			print "\t".$matrix{$gene}{$cell};
		} else {
			print "\t0";
		}	
	}
	print "\n";
}

exit(0);









