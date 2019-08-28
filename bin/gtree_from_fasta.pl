#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw (min max sum);
use Getopt::Long;
use lib '/home/tomfy/Orthologger/lib/';
use TomfyMisc qw ' fasta2seqon1line ';
use constant MISSING_DATA => '3';
use Math::GSL::RNG  qw( :all );

no warnings 'recursion';

use File::Basename 'dirname';
use Cwd 'abs_path';
my ( $bindir, $libdir );
BEGIN {     # this has to go in Begin block so happens at compile time
   $bindir =
     dirname( abs_path(__FILE__) ) ; # the directory containing this script
   $libdir = $bindir . '/../lib';
   $libdir = abs_path($libdir); # collapses the bin/../lib to just lib
}
use lib $libdir;

use GenotypeTree;
use GenotypeTreeNode;
use Genotype;

{                               ###########
   my $input_filename = undef;  # input fasta file name.
   my $show_all = 0;    # false -> output info on only the leaf nodes.
   my $compact = 1;
   my $newick_out = 1;
   my $rng_type = $gsl_rng_mt19937;
   my $rng_seed = undef;
   my $p_missing = -1; # prob. of missing data at each snp.

   GetOptions(
              'input_filename=s' => \$input_filename,
              'show_all!' => \$show_all,
              'compact_tree!' => \$compact,
              'newick_out!' => \$newick_out,
              'rng_type=s' => \$rng_type,
              'seed=i' => \$rng_seed,
              'p_missing=f' => \$p_missing,
             );

   # ############## set up the random number generator ##########################################
   if ($rng_type eq 'sys') {
      $rng_type = $gsl_rng_default;
   } elsif ($rng_type eq 'mt') {
      $rng_type = $gsl_rng_mt19937;
   } elsif ($rng_type eq 'lux') {
      $rng_type = $gsl_rng_ranlxd2;
   }
   print STDERR "RNG type: $rng_type \n";
   my $the_rng = (defined $rng_seed)? Math::GSL::RNG->new($rng_type, $rng_seed) : Math::GSL::RNG->new($rng_type) ; # RNG
   # ############################################################################################



   my $input_filename_stem = $input_filename;
   $input_filename_stem =~ s/\.\S+$//; # remove the part after the (last) '.'
   my $input_string = '';
   {
      if (defined $input_filename) {
         if (-f $input_filename) {
            open my $fhin, "<", $input_filename or die "open $input_filename for reading failed.\n";
            while (my $line = <$fhin>) {
               #   print STDERR $line;
               next if($line =~ /^\s*#/);
               $input_string .= $line;
               #  print STDERR $fasta_as_read;
            }
            close $fhin;
         } else {
            die "file $input_filename does not exist.\n";
         }
      } else {
         die "Must specify input file with -input_filename .\n";
      }
   }

   my @genotype_objects = ();
   my $sequence_length = undef;
   if ($input_filename =~ /\.fasta$/) {
      my @fasta_lines = split ("\n", TomfyMisc::fasta2seqon1line($input_string));
      #    $input_string = undef;
      while (@fasta_lines) {
         my $line = shift @fasta_lines;
         #    print $line;
         next if($line =~ /^\s*#/); # skip comment lines
         if ($line =~ /^>(\S+)\s+(\S+)\s+(\S+)/) {
            my ($id, $generation, $pedigree) = ($1, $2, $3);
            my $sequence = shift @fasta_lines;
            $sequence =~ s/\s+//g;
            $sequence_length = length $sequence if(!defined $sequence_length);
            die "Sequence lengths must all be the same.\n" if(length $sequence != $sequence_length);
            $sequence = lose_data($sequence, $the_rng, $p_missing);
            print STDERR "id: $id  seq: $sequence \n";
        #    print STDERR "AAA: $sequence \n";
            #  my @sequence_array = split('', $sequence);
            push @genotype_objects, Genotype->new(
                                                  { id => $id,
                                                    generation => $generation,
                                                    pedigree => $pedigree,
                                                    sequence => $sequence }
                                                 );
         }
      }
   } else {
      die "Input file must be fasta format. \n";
   }
   my $gtree = GenotypeTree->new( { depth => $sequence_length } );
   for my $gobj (@genotype_objects) {
      #     $gobj->lose_data($the_rng, 0.5);
      if (!$compact) {
         $gtree->add_genotype($gobj);
      } else {
         $gtree->add_genotype_compact($gobj);
      }
   }
   # print "\n", $gtree->as_string(!$show_all), "\n";
print "root ids:  ", $gtree->root()->id_as_string(), "\n";
   print $gtree->as_newick(), "\n" if($newick_out);

   for my $gobj (@genotype_objects) {
      $gtree->search($gobj);
   }

}                               # end main

######################

sub lose_data{
   my $sequence = shift;
   my $rng = shift;
   my $p_missing = shift;

#   print "$sequence \n";
   my $L = length $sequence;
   for my $i (0..$L-1) {
      if (gsl_rng_uniform($rng->raw()) < $p_missing) {
         substr($sequence, $i, 1, MISSING_DATA);
      }
   }
#  print "$sequence \n";
   return $sequence;
}
