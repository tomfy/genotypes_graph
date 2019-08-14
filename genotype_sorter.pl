#!/usr/bin/perl
use strict;
use warnings;
use Readonly;
use List::Util qw (min max sum);
use GenotypeGraph;
use GenotypeGraphNode;
use Genotype;
use Getopt::Long;
use lib '/home/tomfy/Orthologger/lib/';
use TomfyMisc qw ' fasta2seqon1line ';

# Readonly my  $BIG_NUMBER => 1_000_000_000;
#use constant BIG_NUMBER => 1000000000;

{                               ###########

   my $n_nearest_to_keep =  5; # for each genotype make this many directed edges in graph, to the $n_nearest_to_keep closest other nodes.
   my $fasta_input_filename = undef;
   GetOptions(
              'fasta=s' => \$fasta_input_filename,
              'nearest=i' => \$n_nearest_to_keep, # e.g. '*.newick'
             );

   my $fasta_as_read = '';
   {
      if (defined $fasta_input_filename) {
         if (-f $fasta_input_filename) {
            open my $fhin, "<", $fasta_input_filename or die "open $fasta_input_filename failed.\n";
            while (my $line = <$fhin>) {
            #   print STDERR $line;
               next if($line =~ /^\s*#/);
               $fasta_as_read .= $line;
             #  print STDERR $fasta_as_read;
            }
            close $fhin;
         } else {
            die "file $fasta_input_filename does not exist.\n";
         }
      } else {                  # read from stdin.
         while (my $line = <>) {
            next if($line =~ /^\s*#/);
            $fasta_as_read .= $line;
         }
      }

   }
# print $fasta_as_read, "\n";
   my $fasta_string = TomfyMisc::fasta2seqon1line($fasta_as_read);

   # ############# find all the distances and store in a hash.  ####################
   # keys are node ids, and values are hash refs of node id / distance pairs;
   # i.e. $idA__idB_distance{$idA}->{$idB} is the distance between nodes with ids $idA and $idB
   #   my %idA__idB_distance = %{ store_distances( \%id_seq ) };
   # my $genotype_graph = GenotypeGraph->new(
   #                                         'id_sequence' => \$id_seq,
   #                                         'id12dist' => \%idA__idB_distance, 
   #                                         'n_edges' => $n_nearest_to_keep);
   #


my $genotype_graph = GenotypeGraph->new(
                                           { fasta => $fasta_string,
                                             n_edges => $n_nearest_to_keep }
                                          );

   print $genotype_graph->as_string();

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
