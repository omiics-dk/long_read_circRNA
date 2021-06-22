#!/usr/bin/perl

use strict;
use warnings;
use 5.010;
use List::Util qw[min max];

my $col1;
my $col2;
my $col3;
my $col4;
my $col5;
my $col6;
my $col7;
my $col8;
my $strand;
my $name = "Veno";
my $col11;
my $read_start;
my $read_end;
my $chr;
my $col15;
my $start;
my $end;
my @rest;
my $first_strand;
my $first_name = "Morten";
my $first_read_start;
my $first_read_end;
my $first_chr;
my $first_start;
my $first_end;
my @first_rest;

my $count = 0;
my $line;
my $first_line;
my $min;
my $max;
my $abs;
my $internal_count = 0;
my $end_to_start;
my $start_to_end;
my $type;
my $prev_type;

while (<>) {
	chomp();
	$count++;
	$line = $_;
	# First 5 lines in psl format are header. They are just printed.
	# From line 6 we parse data.
	if ( $count > 5 ) {
		($col1, $col2, $col3, $col4, $col5, $col6, $col7, $col8, $strand, $name, $col11, $read_start, $read_end, $chr, $col15, $start, $end, @rest) = split ("\t",$line);
		$internal_count++;
		if ( $name eq $first_name ){
			$min = min $first_start,$start;
			$max = max $first_end,$end;
			$abs = $max - $min;
			# Has to be on the same chromosome and within a megabase
			if ( $first_chr eq $chr && $abs < 1000000 ) {
				if ( $internal_count > 1 ){
					if ( $first_strand ne $strand ) {
						$type = "Not_same_strand";
					} elsif ( $strand eq "+" ){
							# Check if the read segments are approximately contiguous, can be either circRNA or linear
							if ( $first_read_end > $read_start && ($read_end - $first_read_start) < 50 ) {
								if ( $first_start < $start ) {
									$type = "circRNA";
                                                                } elsif ( $first_start > $start) {
                                                                        $type = "linear";
								}
							} elsif ($read_end > $first_read_start && ($first_read_end - $read_start) < 50 ) {
								if ( $first_start > $start ) {
									$type = "circRNA";
								} elsif ( $first_start < $start) {
									$type = "linear";
								}
							} else {
								$type = "ambiguous";
							}
					} elsif ( $strand eq "-" ) {
                                                $end_to_start = abs($first_read_end - $read_start);
                                                $start_to_end = abs($first_read_start - $read_end);
                                                        # Check if the read segments are approximately contiguous, can be either circRNA of linear
                                                        if ( $first_read_end > $read_start && ($read_end - $first_read_start) < 50 ) {
                                                                if ( $first_start > $start ) {
                                                                        $type = "circRNA";
                                                                } elsif ( $first_start < $start) {
                                                                        $type = "linear";
                                                                }
                                                        } elsif ($read_end > $first_read_start && ($first_read_end - $read_start) < 50 ) {
                                                                if ( $first_start < $start ) {
                                                                        $type = "circRNA";
                                                                } elsif ( $first_start > $start) {
                                                                        $type = "linear";
                                                                }
                                                        } else {
                                                                $type = "ambiguous";
                                                        }
					}
					if ( $internal_count == 2 ){
						print "$first_line\t$type\n";
						print "$line\t$type\n";
						$prev_type = $type;
					} elsif ( $internal_count > 2 && !( $prev_type eq "circRNA") ) {
						print "$line\t$type\n";
					}
				}
			}
		} else {
			if ( $count > 6 ) {
				#print "$first_line\t$type\n";
				if ( $internal_count == 2 ){
	                                #print "yes, it gets here";
        	                        $type = "linear_1_fragment";
                	                print "$first_line\t$type\n";
                       		 }

			}
			# Read the first fragment of a read
			$first_line = $line;
			($col1, $col2, $col3, $col4, $col5, $col6, $col7, $col8, $first_strand, $first_name, $col11, $first_read_start, $first_read_end, $first_chr, $col15, $first_start, $first_end, @first_rest) = split ("\t",$first_line);
			$internal_count = 1;
			$type = "Potential_multi-round_circRNA";
		}
	} else {
		# This is where the first 5 lines get printed
		print "$line\n";
	}
}
