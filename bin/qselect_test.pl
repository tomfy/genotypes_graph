#!/usr/bin/perl -w
use strict;
use warnings;
use List::Util qw ( min max sum shuffle );
use Time::HiRes qw( gettimeofday );

{
  my $N = shift // 100;
  my $k = shift // 10;
  my $alg = shift // 'old';
  my $seed = shift // 123457;

  srand($seed);
  
  my @list = ();
  for (1..1) {
    for my $i (0..$N) {
      my $x = $i; rand(100);
      push @list, ($x, $x) ;
    }
  }
  @list = shuffle(@list);

  my $t1 = gettimeofday();
  my $s = '';
  for (1..20) {
    #  my $p = qselect(\@list, $k);
    #  print join ' ', map( qselect(\@list, $_), (1 .. $k) ),  "\n";
    #  print join(' ', smallest_k(\@list, $k)), "\n";
    if ($alg eq 'old') {
      my $smallest_k = smallest_k(\@list, $k);
      print "A: ", join(' ', sort {$a <=> $b } @$smallest_k), "\n";
    } else {
      my @k_best = qselect_x(\@list, $k);
      print "A: ", join(' ', sort { $a <=> $b } @k_best), "\n";
    }
  }
  my $t2 = gettimeofday();
  print "time: ", $t2 - $t1, "\n";
}

sub qselect_x{
  my ($list, $k) = @_;
  #  print "list: ", join(' ', @$list), "  k: $k \n";
  my $pivot_index = int rand @{ $list } -1; #int (0.5* (rand @{ $list } + rand @{$list})) - 1;
  my $pivot = $list->[$pivot_index];

  my @lefts = (); my @rights = (); my @equals = ();
  if (1) {
    @lefts  = grep { $_ < $pivot } @$list;
    @rights = grep { $_ > $pivot } @$list;
    # my @equals = grep { $_ == $pivot } @$list;
    # my @equals = ();
    if (@lefts + @rights + 1 == scalar @$list) {
      push @equals, $pivot;
    } else {
      @equals = grep { $_ == $pivot } @$list;
    }
  } else {
 
    for (@$list) {
      if ($_ < $pivot) {
	push @lefts, $_;
      } elsif ($_ > $pivot) {
	push @rights, $_;
      } else {
	push @equals, $_;
      }
    }
  }
  if ($k < @lefts) {  # kth will be in @lefts, but lefts has too many.
    return qselect_x(\@lefts, $k);
  } elsif ($k > @lefts + @equals) { # kth will be in @rights
    return ((@lefts, @equals), qselect_x(\@rights, $k - @lefts - @equals));
  } elsif ($k == @lefts) {	# done 
    return @lefts;
  } elsif ($k <= @lefts + @equals) { # just @lefts plus 1 or more from @equals
    push @lefts, @equals[0..$k-@lefts-1];
    return @lefts;
  } else {
    die "???\n";
  }
}

sub qselect{ # doesn't work right when elements present more than once in array
  my ($list, $k) = @_;
  my $pivot_index = int rand @{ $list } - 1;
  my $pivot = $list->[$pivot_index];
  my @lefts  = grep { $_ < $pivot } @$list;
  my @rights = grep { $_ > $pivot } @$list;
  #  my $n_equals = @lists - (@lefts + @rights); # number of elements == to $pivot.
  if ($k <= @lefts) {
    return qselect(\@lefts, $k);
  } elsif ($k > @lefts+1) {
    return qselect(\@rights, $k - @lefts - 1);
  } else {			# @lefts == $k-1, i.e. kth is == pivot
    return $pivot;		# , $pivot_index);
  }
}

sub smallest_k{
  my ($list, $k) = @_;
  print "k: $k    list:  ", join(' ', @$list), "\n";
  my $kth_smallest = qselect($list, $k);
  my @smaller_elements = ();
  my @equal_elements = ();
  for (@$list) {
    if ($_ < $kth_smallest) {
      push @smaller_elements, $_;
    } elsif ($_ == $kth_smallest) {
      push @equal_elements, $_;
    }
  }
  my $n_smaller = scalar @smaller_elements;
  my @smallest_k = (@smaller_elements, @equal_elements[0..$k-$n_smaller-1]);
  return \@smallest_k;
}
