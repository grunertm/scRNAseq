#!/usr/bin/perl -w

use strict; 
use warnings;
use lib '.';

my @temp = ();

#Isl1 target region: chr13:116,295,000-116,315,000
for(my $i=116295000; $i<=116315001;$i++){
	print "chr13"."\t".$i."\t".($i+1)."\n";
}

exit(0);










