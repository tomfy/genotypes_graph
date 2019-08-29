#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw (min max sum);

use File::Basename 'dirname';
use Cwd 'abs_path';
my ( $bindir, $libdir );
BEGIN {     # this has to go in Begin block so happens at compile time
  $bindir =
    dirname( abs_path(__FILE__) ) ; # the directory containing this script
  $libdir = $bindir . '/../lib';
  $libdir = abs_path($libdir);	# collapses the bin/../lib to just lib
}
use lib $libdir;

use GenotypeGraph;
use GenotypeGraphNode;
use Genotype;
use Getopt::Long;
use lib '/home/tomfy/Orthologger/lib/';
use TomfyMisc qw ' fasta2seqon1line ';

{                               ###########

   my $n_nearest_to_keep =  5; # for each genotype make this many directed edges in graph, to the $n_nearest_to_keep closest other nodes.
   my $input_filename = undef; # input fasta file name.
   my $output_distance_matrix = 1; # whether to output a distance matrix ( .dmatrix filename ending)
   my $output_graph = 1; # whether to output the graph (.grph filename ending0
   my $multiplier = 10000; # controls # significant digits. 10000 -> 0.6492361... is output as 6492
   my $show_sequence = 0; # if true, will output the sequence at the end of line.
   GetOptions(
              'input_filename=s' => \$input_filename,
              'nearest=i' => \$n_nearest_to_keep, # e.g. '*.newick'
              'distance_matrix_out!' => \$output_distance_matrix,
              'graph_out!' => \$output_graph,
              'multiplier=i' => \$multiplier,
              'sequence_out!' => \$show_sequence,
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

   my $genotype_graph;

   if ($input_filename =~ /\.fasta$/) {
      my $fasta_string = TomfyMisc::fasta2seqon1line($input_string);
      $genotype_graph = GenotypeGraph->new(
                                           { fasta => $fasta_string,
                                             n_edges => $n_nearest_to_keep }
                                          );
   } elsif ($input_filename =~ /\.grph/) {
      $genotype_graph = GenotypeGraph->new(
                                           { idnnd => $input_string,
                                             n_edges => $n_nearest_to_keep }
                                          );
   }

   if ($output_graph) {
      my $graph_out_filename = $input_filename_stem . '.grph';
      open my $fhout, ">", $graph_out_filename or die "Couldn't open $graph_out_filename for writing.\n";
      print $fhout $genotype_graph->as_string($show_sequence);
      close $fhout;
   }
   if ($output_distance_matrix) {
      my $distance_out_filename = $input_filename_stem . '.dmatrix';
      open my $fhout, ">", $distance_out_filename or die "Couldn't open $distance_out_filename for writing.\n";
      print $fhout $genotype_graph->distance_matrix_as_string($multiplier);
      close $fhout;
   }

}                               # end main

################ subroutines ############################

# find all the distances and store in a hash.
# keys are node ids, and values are hash refs of node id : distance pairs;
# i.e. $idA__idB_distance{$idA}->{$idB} is the distance between nodes with ids $idA and $idB
# sub store_distances{
#    my $id_seq = shift;

#    my @ids = keys %$id_seq;

#    my %idA__idB_distance = ();  # hash of hash refs
#    while (my ($i1, $id1) = each @ids) {
#       #   print STDERR "$i1 $id1 \n";
#       my $s1 = $id_seq->{$id1};
#       for (my $i2 = $i1+1; $i2 < scalar @ids; $i2++) {
#          my $id2 = $ids[$i2];
#          #    print STDERR "    $i2 $id2 \n";
#          my $s2 = $id_seq->{$id2};
#          my ($d_ij, $both_snp_count) = distance($s1, $s2);
#          #   print "$id1  $id2  $d_ij \n";
#          my $d = ($both_snp_count > 0)? $d_ij/$both_snp_count : BIG_NUMBER; #
#          if (!exists $idA__idB_distance{$id1}) {
#             $idA__idB_distance{$id1} = {$id2 => $d_ij};
#          } else {
#             $idA__idB_distance{$id1}->{$id2} = $d_ij;
#          }
#          if (!exists $idA__idB_distance{$id2}) {
#             $idA__idB_distance{$id2} = {$id1 => $d_ij};
#          } else {
#             $idA__idB_distance{$id2}->{$id1} = $d_ij;
#          }
#       }
#    }
#    return \%idA__idB_distance;
# }

# sub prune_to_n_closest{
#    my $idA__idB_dist = shift;
#    my $n_keep = shift;

# }


# sub distance{
#    my $seq1 = shift;
#    my $seq2 = shift;

#    die "Sequence lengths are different - bye.\n" if(length $seq1 != length $seq2);

#    my $distance = 0;
#    my $count_both = 0; # count of snps with data present in both sequences
#    my $count_missing = 0; # count of snps with data absent in one or both sequences
#    for (my $i=0; $i < length $seq1; $i++) {
#       my ($c1, $c2) = (substr($seq1, $i, 1), substr($seq2, $i, 1));
#       if ($c1 eq '-'  or $c2 eq '-') {
#          $count_missing++;
#       } else {
#          $count_both++;
#          if ($c2 != $c1) {
#             if ($c1 == 0) {
#                $distance += $c2 - $c1;
#             } elsif ($c1 == 1) {
#                $distance += 1;
#             } elsif ($c1 == 2) {
#                $distance += $c1 - $c2;
#             } else {
#                die "c1 has unhandled value: $c1 \n";
#             }
#          }
#       }
#    }
#    return ($distance, $count_both);
# }
