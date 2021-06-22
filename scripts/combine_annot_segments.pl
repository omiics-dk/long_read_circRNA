#!/usr/bin/perl

use strict;
use warnings;
use 5.010;
#use List::Util qw[min max];
#use Math::Round;

my $col1;
my $col2;
my $col3;
my $col4;
my $col5;
my $col6;
my $col7;
my $col8;
my $col9;
my $col10;
my $col11;
my $name;
my $first_name = "Morten";
my $first_col2;
my $sum_col5 = 0;
my $sum_col6 = 0;
my $sum_col7 = 0;
my $sum_col8 = 0;
my $sum_col9 = 0;
my $sum_col10 = 0;
my $sum_col11 = 0;

my $count = 0;
my $line;
my $first_line;
my $internal_count = 0;
my $round5;
my $round6;
my $round7;
my $round8;
my $round9;
my $round10;
my $round11;
my $prev_col1;
my $prev_col3;

while (<>) {
	chomp();
	$count++;
	$line = $_;
	($col1, $col2, $col3, $name, $col5, $col6, $col7, $col8, $col9, $col10, $col11) = split ("\t",$line);
	$internal_count++;
#	if ( $count == 1 ){
#		print "$line\n";
#		$internal_count = 0;
#	}
	if ( $internal_count == 1 ){
		$first_line = $line;
		$first_col2 = $col2;
		$first_name = $name;
	}
	if ( $name eq $first_name ){
		$sum_col5 = $sum_col5 + $col5;
		$sum_col6 = $sum_col6 + $col6;
		$sum_col7 = $sum_col7 + $col7;
		$sum_col8 = $sum_col8 + $col8;
		$sum_col9 = $sum_col9 + $col9;
#		$sum_col10 = $sum_col10 + $col10;
#		$sum_col11 = $sum_col11 + $col11;
		$prev_col1 = $col1;
		$prev_col3 = $col3;
	} elsif ( $count > 1 ) {
		$internal_count = 1;
		$round5 = int( $sum_col5 + 0.5 );
		$round6 = int( $sum_col6 + 0.5 );
		$round7 = int( $sum_col7 + 0.5 );
		$round8 = int( $sum_col8 + 0.5 );
		$round9 = int( $sum_col9 + 0.5 );
#		$round10 = int( $sum_col10 + 0.5 );
#		$round11 = int( $sum_col11 + 0.5 );
		printf "$prev_col1\t$first_col2\t$prev_col3\t$first_name\t$round5\t$round6\t$round7\t$round8\t$round9\n";
		$sum_col5 = $col5;
		$sum_col6 = $col6;
		$sum_col7 = $col7;
		$sum_col8 = $col8;
		$sum_col9 = $col9;
#		$sum_col10 = $col10;
#		$sum_col11 = $col11;
		$first_line = $line;
		$first_col2 = $col2;
		$first_name = $name;
	}
}
