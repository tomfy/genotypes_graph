package Genotype;
use Moose;
use namespace::autoclean;
use Carp;
use List::Util qw ( min max sum );

use constant MISSING_DATA => 3;
#Readonly my  $BIG_NUMBER => 1_000_000_000;

# a class for genotypes
# each has an identifier (integer),
# and a genotype represented as an array ref
# containing numbers 0,1,2,3
# 0 = AA, 1 = Aa, 2 = aa, 3 = missing.

has id => (
           isa => 'Int',
           is => 'ro',
           required => 1,
          );

has genotype => (
                 isa => 'ArrayRef', # e.g. [0,1,1,1,0,2,1,1,2,3,1,1,0,1,0,3,2];
                 is => 'ro',
                 required => 1,
                );

has pedigree => (
                 isa => 'Str',
                 is => 'ro',
                 default => sub { undef },
                );

around BUILDARGS => sub {
   my $orig = shift;
   my $class = shift;

   if (@_ == 1  and  (!ref $_[0]) ) { # argument should be string with >id generation pedigree \n sequence.
      my $arg = $_[0];
    #  print "XXXXX: $arg \n";
      my ($idline, $sequence) = split(/\n/, $arg);
      my ($id, $generation, $pedigree) = (undef, undef, undef);
      if ($idline =~ /^>(\S+)\s+(\d+)\s+(\S+)\s*$/x) {
         ($id, $generation, $pedigree) = ($1, $2, $3);
      } elsif ($idline =~ /^ > (\S+) /x) {
         $id = $1;
      } else {
         die "in Genotype BUILDARGS argument $arg does not have expected >id format.\n";
      }
      $sequence =~ s/\s+//gx;
      my @sequence_array = split(//, $sequence);
      return {id => $id, genotype => \@sequence_array, pedigree => $pedigree};
   }
   # otherwise, don't modify arguments, should be array ref with id, genotype, and optionally pedigree
};


sub distance{ # calculate distance between this genotype obj. and another
   my $self = shift;
   my $other_genotype = shift;
   my $this_gt = $self->genotype();            # array ref of 0,1,2,3
   my $other_gt = $other_genotype->genotype(); # array ref of 0,1,2,3

   my $distance = 0;
   my $count_both = 0; # count of snps with data present in both sequences
   my $count_missing = 0; # count of snps with data absent in one or both sequences
   if (scalar @$this_gt == scalar @$other_gt) {
      while (my($i, $c1) = each @{$this_gt}) {
         my $c2 = $other_gt->[$i];
         if ( ($c1 == MISSING_DATA)  or ($c2 == MISSING_DATA) ) {
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


############################################

__PACKAGE__->meta->make_immutable;

1;
