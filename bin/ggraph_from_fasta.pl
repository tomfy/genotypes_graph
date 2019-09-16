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
   my $other_fasta = undef;
   my $output_distance_matrix = 1; # whether to output a distance matrix ( .dmatrix filename ending)
   my $output_graph = 1; # whether to output the graph (.grph filename ending0
   my $multiplier = 10000; # controls # significant digits. 10000 -> 0.6492361... is output as 6492
   my $show_sequence = 0; # if true, will output the sequence at the end of line.
   my $n_extras = 0; # number of extra 'neighbors' to give each node, in addition to the $n_nearest_to_keep nearest nodes.
   my $seed = 1234579;
   my $search_pq_size = 20;
   my $error_prob = 0;

   GetOptions(
              'input_filename|fasta1|f1=s' => \$input_filename,
              'fasta2|f2=s' => \$other_fasta,
              'nearest=i' => \$n_nearest_to_keep, # e.g. '*.newick'
              'distance_matrix_out!' => \$output_distance_matrix,
              'graph_out!' => \$output_graph,
              'multiplier=i' => \$multiplier,
              'sequence_out!' => \$show_sequence,
	      'extras=i' => \$n_extras,
	      'seed=i' => \$seed,
	      'pq_size=i' => \$search_pq_size,
	      'error_prob=f' => \$error_prob,
             );

   if ($seed > 0) {
      srand($seed);
   } else {
      srand();
   }
   my $input_filename_stem = $input_filename;
   $input_filename_stem =~ s/\.\S+$//; # remove the part after the (last) '.'
   my $input_string = '';
   if (defined $input_filename) {
      $input_string = file_to_string($input_filename);
   } else {
      die "Must specify input file with -input_filename .\n";
   }

   my $genotype_graph;

   if ($input_filename =~ /\.fasta$/) {
      my $fasta_string = TomfyMisc::fasta2seqon1line($input_string);
      $genotype_graph = GenotypeGraph->new(
                                           { fasta => $fasta_string,
                                             n_near => $n_nearest_to_keep,
					     n_extras => $n_extras,
					   }
                                          );
   } elsif ($input_filename =~ /\.grph/) {
      $genotype_graph = GenotypeGraph->new(
                                           { idnnd => $input_string,
                                             n_near => $n_nearest_to_keep,
					     n_extras => $n_extras,
					   }
                                          );
   }
   my $other_fasta_string = (defined $other_fasta)?  TomfyMisc::fasta2seqon1line(file_to_string($other_fasta)) : $input_string;
   my $other_id_gobjs = GenotypeGraph::fasta_string_to_gobjs($other_fasta_string);

   #  for(my $i = 0; $i < 100; $i++){
   #  for my $genotypegraph_node_to_clone (values %{$genotype_graph->nodes()}){ # -> genotype();
   for my $g_to_search_for (values %$other_id_gobjs) {
      #  my $g_to_search_for = $genotypegraph_node_to_clone->genotype()->clone(id => $genotypegraph_node_to_clone->id() + 1000000);
      $g_to_search_for->{id} += 100000;
      $g_to_search_for->add_noise($error_prob);
      $genotype_graph->search_for_best_match($g_to_search_for, $search_pq_size);
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


sub file_to_string{
   my $filename = shift;
   my $input_string = '';
  if (-f $filename) {
      open my $fhin, "<", $filename or die "open $filename for reading failed.\n";
      while (my $line = <$fhin>) {
         #   print STDERR $line;
         next if($line =~ /^\s*#/);
         $input_string .= $line;
         #  print STDERR $fasta_as_read;
      }
      close $fhin;
   } else {
      die "file $filename does not exist.\n";
   }
   return $input_string;
}
