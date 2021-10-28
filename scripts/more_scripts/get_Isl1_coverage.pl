#!/usr/bin/perl -w

use strict; 
use warnings;
use lib '.';

my @temp = ();
my @temp2 = ();
my $in;
my $in2;
my %coverage;

my $sum = 0;

#coverage file
open ($in, '<', $ARGV[0]) or die "Input file $ARGV[0] not found!\n";
while(<$in>) {
	chomp;
	@temp = split("\t| ", $_);
	my @sample = split(/\./, $temp[0]);
	open ($in2, '<', $temp[0]) or die "Input file $temp[0] not found!\n";
	while(<$in2>) {
		chomp;
		@temp2 = split("\t", $_);	
		$coverage{$temp2[1]}{$sample[0]} = $temp2[3];
	}
	close $in2;
}
close $in;

# E8.1 (E8_1_Set_A_*) 
# E8.2 (E8_2_Set_C_*)
# E9_anterior_1 (E9_anterior_1_Set_D_*) 
# E9_posterior_4 (Lib_E9_posterior_4_Set_B_*)

print "#Pos\tE8.1\tE8.2\tE0_ant\tE9_pos\n";

for my $pos ( sort keys %coverage ) {

	my $E8_1_sum = 0;
	my $E8_1_n = 0;
	
	my $E8_2_sum = 0;
	my $E8_2_n = 0;
	
	my $E9_ant_sum = 0;
	my $E9_ant_n = 0;
	
	my $E9_pos_sum = 0;
	my $E9_pos_n = 0;
	
	for my $sample ( sort keys %{$coverage{$pos}}) {
	
		if($sample =~ m/E8_1/) {
			$E8_1_n += 1;
			$E8_1_sum += $coverage{$pos}{$sample}
		}
		
		if($sample =~ m/E8_2/) {
			$E8_2_n += 1;
			$E8_2_sum += $coverage{$pos}{$sample}
		}
		
		if($sample =~ m/anterior/) {
			$E9_ant_n += 1;
			$E9_ant_sum += $coverage{$pos}{$sample}
		}
		
		if($sample =~ m/posterior/) {
			$E9_pos_n += 1;
			$E9_pos_sum += $coverage{$pos}{$sample}
		}
	}
	
	print $pos;
	print "\t";
	print ($E8_1_sum/$E8_1_n);
	print "\t";
	print ($E8_2_sum/$E8_2_n);
	print "\t";
	print ($E9_ant_sum/$E9_ant_n);
	print "\t";
	print ($E9_pos_sum/$E9_pos_n);
	print "\n";
	
}



exit(0);









