package Genotype;
use Moose;
use namespace::autoclean;
use Carp;
use List::Util qw ( min max sum );

use constant MISSING_DATA => 'X';

# a class for genotypes
# each has an identifier (integer),
# and a genotype represented as a string
# containing characters 0, 1, 2, MISSING_DATA
# 0 = AA, 1 = Aa, 2 = aa.

has id => (
           isa => 'Int',
           is => 'ro',
           required => 1,
          );

has generation => (
                   isa => 'Int',
                   is => 'ro',
                   default => sub { undef },
                  );

has pedigree => (
                 isa => 'Str',
                 is => 'ro',
                 default => sub { undef },
                );

has sequence => (
                isa => 'Str',
                is => 'ro',
                required => 1,
                );

around BUILDARGS => sub {
   my $orig = shift;
   my $class = shift;

   if (@_ == 1) {
      if (!ref $_[0]) { # argument should be string with >id generation pedigree \n sequence.
         my $arg = $_[0];
         my ($idline, $sequence_string) = split(/\n/, $arg);
         my ($id, $generation, $pedigree) = (undef, undef, undef);
         if ($idline =~ /^>(\S+)\s+(\d+)\s+(\S+)\s*$/x) {
            ($id, $generation, $pedigree) = ($1, $2, $3);
         } elsif ($idline =~ /^ > (\S+) /x) {
            $id = $1;
         } else {
            die "in Genotype BUILDARGS argument $arg does not have expected >id format.\n";
         }
         $sequence_string =~ s/\s+//gx;
      #   my @sequence_array = split(//, $sequence_string);
         return {id => $id, sequence => $sequence_string, generation => $generation, pedigree => $pedigree};
      } else {
         return $_[0];          # should be hash ref of arguments
      }
   }
   # otherwise, don't modify arguments, should be array ref with id, sequence, and optionally generation and pedigree
};


sub distance{ # calculate distance between this genotype obj. and another
   my $self = shift;
   my $other_genotype = shift;
   my $this_gt = $self->sequence();            # array ref of 0,1,2,3
   my $other_gt = $other_genotype->sequence(); # array ref of 0,1,2,3

   my $distance = 0;
   my $count_both = 0; # count of snps with data present in both sequences
   my $count_missing = 0; # count of snps with data absent in one or both sequences
   if (scalar @$this_gt == scalar @$other_gt) {
      while (my($i, $c1) = each @{$this_gt}) {
         my $c2 = $other_gt->[$i];
         if ( ($c1 eq MISSING_DATA)  or ($c2 eq MISSING_DATA) ) {
            $count_missing++;
         } else {
            $count_both++;
            if ($c2 != $c1) {
               if ($c1 == 0) {
                  $distance += $c2 - $c1;
               } elsif ($c1 == 1) {
                  $distance += 1;
               } elsif ($c1 == 2) {
                  $distance += $c1 - $c2;
               } else {
                  die "c1 has unhandled value: $c1 \n";
               }
            }
         }
      }
   }
   return ($distance, $count_both, $count_missing);
}

sub mean{
   my $self = shift;
   my $other_genotype = shift;
   my $this_gt = $self->sequence();            # array ref of 0,1,2,3
   my $other_gt = $other_genotype->sequence(); # array ref of 0,1,2,3
   my @mean_gt = ();
   if (scalar @$this_gt == scalar @$other_gt) { # check that lengths are equal
      while (my($i, $c1) = each @{$this_gt}) {
         my $c2 = $other_gt->[$i];
         if ($c1 eq MISSING_DATA) {
            push @mean_gt, $c2;
         } elsif ($c2 eq MISSING_DATA) {
            push @mean_gt, $c1;
         } else {
            push @mean_gt, 0.5*($c1 + $c2);
         }
      }
   }
   return \@mean_gt;
}

############################################

__PACKAGE__->meta->make_immutable;

1;
