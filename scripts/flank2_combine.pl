#!/usr/bin/perl

use strict;
use warnings;
use 5.010;
use List::Util qw[min max];

my $col1;
my $col2;
my $name;
my $chr;
my $start;
my $end;
my $strand;
my $count = 0;
my $spliceSite;
my $firstSeq;
my $len;
my $genomeStrand;

while (<>) {
	chomp();
	$count++;
	($col1, $col2) = split ("\t",$_);
	($name, $chr, $start, $end, $strand) = split ("~",$col1);
	$len = $end - $start;
	if ( $count == 1 ) {
		$spliceSite = "";
		$firstSeq = $col2;
	} elsif ( $count == 2 ) {
		$count = 0;
		$spliceSite = "$firstSeq$col2";
		if ($spliceSite eq "AGGT"){
			$genomeStrand = "+";
		} elsif ($spliceSite eq "ACCT"){
			$genomeStrand = "-";
		} else {
			$genomeStrand = "Unknown";
		}
		printf "$chr\t$start\t$end\t$name\t$len\t$genomeStrand\t$spliceSite\n";
	}
}

