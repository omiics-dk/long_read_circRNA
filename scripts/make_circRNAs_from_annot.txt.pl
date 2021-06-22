#!/usr/bin/perl

use strict;
use warnings;
use 5.010;
use List::Util qw[min max];

my $col1;
my $col2;
my $col3;
my $col4;
my $length;
my $gene;
my $exon;
my $est;
my $intron;
my $name = "Morten";
my $prev_length;
my $prev_gene;
my $prev_exon;
my $prev_est;
my $prev_intron;
my $prev_name = "Veno";

my $region;
my $chrStart;
my $chr;
my $start;
my $end;
my $prev_chr;
my $prev_start;
my $prev_end;

my $start_min;
my $end_max;
my $length_sum = 0;
my $gene_sum = 0;
my $exon_sum = 0;
my $est_sum = 0;
my $intron_sum = 0;

my $count = 0;
my $line;
my $first_line;
my $internal_count = 0;
my $do_print = 0;
my $for_printer;

while (<>) {
	chomp();
	$count++;
	$line = $_;
	if ( $count > 1 ) {
		($col1, $col2, $col3, $col4, $length, $gene, $exon, $est, $intron) = split ("\t",$line);
		($name, $region) = split ("~",$col4);
		($chrStart, $end) = split ("-",$region);
		($chr, $start) = split (":",$chrStart);
		if ( $prev_name ne $name ) {
			$internal_count = 0;
			if ( $do_print == 1 ){
				# Print
				print "$for_printer\n";
			}
		}
	}
	$internal_count++;
	if ( $count > 2 && $internal_count == 1 ){
		$prev_name = $name;
		$prev_chr = $chr;
		$prev_start = $start;
		$prev_end = $end;
		$prev_length = $length;
		$prev_gene = $gene;
		$prev_exon = $exon;
		$prev_est = $est;
		$prev_intron = $intron;
		$do_print = 0;
	} elsif ( $internal_count == 2 ){
		if ($prev_name = $name ) {
			$start_min = min $prev_start,$start;
			$end_max = max $prev_end,$end;
			$length_sum = eval ($prev_length + $length);
			$gene_sum = eval ($prev_gene + $gene);
			$exon_sum = eval ($prev_exon + $exon);
			$est_sum = eval ($prev_est + $est);
			$intron_sum = eval ($prev_intron + $intron);
			# Combine for printing:
			$for_printer = "$chr\t$start_min\t$end_max\t$name\t$length_sum\t$gene_sum\t$exon_sum\t$est_sum\t$intron_sum";
			$do_print = 1;
#			print "$for_printer\n";
		} else {
			$do_print = 0;
		}
#		print "$do_print\t$internal_count\t$name\t$prev_name\n";
	} elsif ( $internal_count > 2 ) {
		$do_print = 0;
	}
#	print "$do_print\t$internal_count\t$name\t$prev_name\n";
	if ( $prev_name ne $name ) {
#		print "$do_print\t$internal_count\t$name\t$prev_name\n";
		$internal_count = 1;
                $prev_name = $name;
                $prev_chr = $chr;
                $prev_start = $start;
                $prev_end = $end;
                $prev_length = $length;
                $prev_gene = $gene;
                $prev_exon = $exon;
		$prev_est = $est;
                $prev_intron = $intron;
		$do_print = 0;
	}
}
if ( $internal_count == 2 && $do_print == 1) {
	# Print edge
	print "$for_printer\n";
}
