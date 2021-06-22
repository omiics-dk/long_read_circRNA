#!/usr/bin/perl

use strict;
use warnings;
use 5.010;
#use List::Util qw[min max];

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
my $col12;
my $col13;
my $col14;
my $col15;
my $col16;
my $col17;
my $col18;
my $col19;
my $col20;
my $col21;
my $name;
my @blockStarts;
my @blockSizes;
my $count = 0;
my $line;
my $rgb;

while (<>) {
	chomp();
	$count++;
	$line = $_;
	# First 5 lines in psl format are header. They are just skipped.
	# From line 6 we parse data.
	if ( $count > 5 ) {
		($col1,$col2,$col3,$col4,$col5,$col6,$col7,$col8,$col9,$col10,$col11,$col12,$col13,$col14,$col15,$col16,$col17,$col18,$col19,$col20,$col21 ) = split ("\t",$line);
		@blockSizes = split (",",$col19);
		@blockStarts = split (",",$col20);

		$rgb = "0,0,0";
		$name = "$col10~$col14:$col16-$col17";
		print "$col10\t$col12\t$col13\t$name\t$col1\t$col9\t$col12\t$col13\t$rgb\t$col18\t";

		# Because BLAT psl format handles query start and end in a strange way, we need to do adjustments for negative strand reads
                if ( $col9 eq "-" ) {
			my @blockSizesRev = reverse @blockSizes;
                        foreach (@blockSizesRev) {
                                print "$_,";
                        }
			print "\t";


                        my @blockEnds = map { $blockStarts[$_] + $blockSizes[$_] } 0..$#blockStarts;
                        my @blockEndsRev = reverse @blockEnds;
                        foreach (@blockEndsRev) {
                                my $temp = $col11 - $_ - $col12;
                                print "$temp,";
                        }
		} else {
			print "$col19\t";
			foreach (@blockStarts) {
				my $temp = $_ - $col12;
				print "$temp,";
			}
		}
		print "\n";

	}
}
