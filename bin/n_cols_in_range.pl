#!/usr/bin/perl -w
use strict;
use warnings;

my $colstr = shift;
my $range  = shift;

my @cols = split( ',', $colstr );
my ( $min, $max ) = split( ',', $range );

my $total_count = 0;
while (<>) {
    my $count = 0;
    next if (/^\s*#/);
    my @dcols = split( ' ', $_ );
    for my $acol (@cols) {
        my $value = $dcols[ $acol - 1 ];
	$count++ if ( $value >= $min  and  $value <= $max);
    }
	$total_count += $count;
    print "$count $total_count\n"; #  $_";
}
