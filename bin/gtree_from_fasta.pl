#!/usr/bin/perl
use strict;
use warnings;
use Readonly;
use List::Util qw (min max sum);
use Getopt::Long;
use lib '/home/tomfy/Orthologger/lib/';
use TomfyMisc qw ' fasta2seqon1line ';

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
   GetOptions(
              'input_filename=s' => \$input_filename,
             );

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

      while (@fasta_lines) {
         my $line = shift @fasta_lines;
         #    print $line;
         next if($line =~ /^\s*#/); # skip comment lines
         if ($line =~ /^>(\S+)\s+(\S+)\s+(\S+)/) {
            my ($id, $generation, $pedigree) = ($1, $2, $3);
            my $sequence = shift @fasta_lines;
            $sequence =~ s/\s+//g;
            print "XXX:  $id  $generation  $pedigree \n", "$sequence\n";
            $sequence_length = length $sequence if(!defined $sequence_length);
            die "Sequence lengths must all be the same.\n" if(length $sequence != $sequence_length);
            my @sequence_array = split('', $sequence);
            push @genotype_objects, Genotype->new(
                                                  { id => $id,
                                                    generation => $generation,
                                                    pedigree => $pedigree,
                                                    sequence => \@sequence_array }
                                                 );
         }

      }
   } else {
      die "Input file must be fasta format. \n";
   }

   my $gtree = GenotypeTree->new( { depth => $sequence_length } );
   for my $gobj (@genotype_objects) {
      $gtree->add_genotype($gobj);
   }

   print "\n\n[", $gtree->as_string(1), "]\n";

}                               # end main
