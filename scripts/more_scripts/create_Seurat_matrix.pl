#!/usr/bin/perl -w

use strict; 
use warnings;
use lib '.';

my @temp = ();
my @header = ();

my $in;
my $out;

my %matrix;
my %gene_in_matrix;
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
	if($temp[0] =~ m/GENE/) {
		@header = split("\t", $_);	
		next;
	}
	
	for(my $i=1; $i<scalar(@temp); $i++) {
		$gene_in_matrix{$temp[0]} = $IDs{$temp[0]};
		$matrix{$temp[0]}{$header[$i]} = $temp[$i];
	}
}
close $in;



#genes.tsv
open($out, '>', 'genes.tsv');
my $gene_count = 0;
my %genes;
for my $gene ( sort keys %gene_in_matrix ) {
	$gene_count++;
	$genes{$gene} = $gene_count;
	print $out $gene." ";
	print $out $gene_in_matrix{$gene};
	print $out "\n";
}
close $out;


# matrix.mtx
open($out, '>', 'matrix.mtx');
print $out "%%Matrix with TPM expression values\n";
print $out "%\n";
my $cell_count = 0;
for my $gene ( sort keys %matrix ) {
	
	$cell_count = 0;
	for my $cell ( sort keys %{$matrix{$gene}} ) {
		$cell_count++;
		print $out $genes{$gene};
		print $out " ";
		print $out $cell_count;
		print $out " ";
		print $out $matrix{$gene}{$cell};
		print $out "\n";
	}
	
	
}
close $out;

exit(0);









