#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw (min max sum);
use Time::HiRes qw( gettimeofday );
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

  my $input_filename = undef;	# input fasta file name.
  my $other_fasta = undef;
  my $error_prob = 0;

  my $output_distance_matrix = 1; # whether to output a distance matrix ( .dmatrix filename ending)
  my $output_graph = 1; # whether to output the graph (.graph filename ending0
  my $multiplier = 10000; # controls # significant digits. 10000 -> 0.6492361... is output as 6492

  my $n_nearest_to_keep = 20; # for each genotype make this many directed edges in graph, to the $n_nearest_to_keep closest other nodes
  my $n_nearest_for_search = 5;
  my $n_extras = 0; # number of extra 'neighbors' to give each node, in addition to the $n_nearest_to_keep nearest nodes.

  my $do_search = 1; # default is to do search. -nosearch to skip the search.
  my $n_independent_searches = 2;
  my $search_pq_size = 20;
  my $n_futile_rounds = 1;
  my $seed = 1234579;
  my $graph_search_outfilename = 'gr_search_out';
  my $exhaustive_search_outfilename = 'exh_search_out';
  my $do_exhaustive_search = 0;


  GetOptions(
	     'input_filename|fasta1|f1|stem=s' => \$input_filename,
	     'fasta2|f2=s' => \$other_fasta,
	     'error_prob=f' => \$error_prob,

	     'distance_matrix_out|dmatrix!' => \$output_distance_matrix,
	     'graph_out!' => \$output_graph,
	     'multiplier=i' => \$multiplier,

	     'search!' => \$do_search,
	     'keep=i' => \$n_nearest_to_keep, # e.g. '*.newick'
	     'neighbors|nearest=i' => \$n_nearest_for_search,
	     'extras=i' => \$n_extras,
	     'starts=i' => \$n_independent_searches,
	     'pq_size=i' => \$search_pq_size,
	     'rounds=i' => \$n_futile_rounds,
	     'seed=i' => \$seed, # rng seed - but results are not reproducible even with same seed (due to hashes?)
	     'exhaustive_search!' => \$do_exhaustive_search,
	    );

  die "Input file for constructing graph must be specified, is undefined. Bye.\n" if(!defined $input_filename);

  if ($seed > 0) {
    srand($seed);
  } else {
    srand();
  }
  my $input_filename_stem = $input_filename;
  $input_filename_stem =~ s/\.(fasta|graph|dmatrix)$//; # remove the part after the last '.'
  my $fasta1_filename = $input_filename_stem . '.fasta'; # get an array or hash of Genotype objs from this fasta file. 
  my $fasta1_string = (-f $fasta1_filename)?
    file_to_string($fasta1_filename)
    :  print STDERR "# Specified fasta file: $fasta1_filename,  does not exist.\n",
    "Will construct graph from .graph or .dmatrix file\n",
    "No sequence info; no search possible.\n";;

  my $genotype_graph;

  my $t0 = gettimeofday();
  my $input_string;
  if ( -f ($input_filename = $input_filename_stem . '.graph') ) { # .graph file exists - construct graph from it.
    print STDERR "# Constructing graph from .graph file: $input_filename\n";
    print "# Constructing graph from .graph file: $input_filename\n";
    $output_distance_matrix = 0; # turn off writing distance_matrix file.
    $input_string =  file_to_string($input_filename);
    $genotype_graph = GenotypeGraph->new(
					 { fasta_string => $fasta1_string,
					   graph_string => $input_string}
					);

  } elsif ( -f ($input_filename = $input_filename_stem . '.dmatrix') ) {
    print STDERR "# Constructing graph from .dmatrix file: $input_filename\n";
    print "# Constructing graph from .dmatrix file: $input_filename\n";

    $input_string =  file_to_string($input_filename);
    #   print "AAA: [$input_string]\n";
    $genotype_graph = GenotypeGraph->new(
					 {  fasta_string => $fasta1_string,
					    dmatrix_string => $input_string,
					    n_keep => $n_nearest_to_keep,
					    n_near => $n_nearest_for_search,
					    n_extras => $n_extras,
					 }
					);

  } elsif ( -f ($input_filename = $input_filename_stem . '.fasta')  ) {
    print STDERR "# Constructing graph from .fasta file: $input_filename\n";
    print "# Constructing graph from .fasta file: $input_filename\n";
    my $fasta1_string = TomfyMisc::fasta2seqon1line(file_to_string($input_filename));
    # print "[$fasta_string]\n";
    $genotype_graph = GenotypeGraph->new(
					 { fasta_string => $fasta1_string,
					   n_keep => $n_nearest_to_keep,
					   n_near => $n_nearest_for_search,
					   n_extras => $n_extras,
					 }
					);
  }
  my $t1 = gettimeofday();
  printf( STDERR "Done constructing graph. Time to construct: %10.3g\n", $t1-$t0);
  #die "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n";
  my $t1point5 = undef;
  if ($do_search) {
    # ###########  Search  ####################
    if (defined $other_fasta  and  -f $other_fasta) {
      my $other_fasta_string = TomfyMisc::fasta2seqon1line(file_to_string($other_fasta));
      my $other_id_gobjs = GenotypeGraph::fasta_string_to_gobjs($other_fasta_string);
      print STDERR "#  searching for matches to genotypes in file: $other_fasta.\n";
      open my $fhgrout, ">", $graph_search_outfilename or die;
      print $fhgrout "# searching for matches to genotypes in file: $other_fasta.\n";
      print $fhgrout "# graph parameters: n_keep: $n_nearest_to_keep  n_near: $n_nearest_for_search  n_extras: $n_extras \n";
      print $fhgrout "# search parameters: n independent searches $n_independent_searches  pq size: $search_pq_size  n futile rounds: $n_futile_rounds \n";
      #  for(my $i = 0; $i < 100; $i++){
      #  for my $genotypegraph_node_to_clone (values %{$genotype_graph->nodes()}){ # -> genotype();
      print $fhgrout "# number of genotypes searched for: ", scalar keys %$other_id_gobjs, "\n";

      for my $g_to_search_for (values %$other_id_gobjs) {
        #  my $g_to_search_for = $genotypegraph_node_to_clone->genotype()->clone(id => $genotypegraph_node_to_clone->id() + 1000000);
        $g_to_search_for->{id} += 100000;
        $g_to_search_for->add_noise($error_prob);
        print $fhgrout $genotype_graph->search_for_best_match($g_to_search_for, $n_independent_searches, $search_pq_size, $n_futile_rounds), "\n";
      }
      close $fhgrout;
      $t1point5 = gettimeofday();
printf( STDERR "Done searching graph. Time to conduct graph search: %10.3g\n", $t1point5-$t1);

      if($do_exhaustive_search){
      open my $fhexout, ">", $exhaustive_search_outfilename or die;
      for my $g_to_search_for (values %$other_id_gobjs) {
        #  my $g_to_search_for = $genotypegraph_node_to_clone->genotype()->clone(id => $genotypegraph_node_to_clone->id() + 1000000);
	#      $g_to_search_for->{id} += 100000;
	#      $g_to_search_for->add_noise($error_prob);
        print $fhexout $genotype_graph->exhaustive_search($g_to_search_for, 4), "\n";
      }
    }
    }else{
       die "Search requested but no file of sequences to search for specified or file doesn't exist.\n";
    }
  }
  my $t2 = gettimeofday();
  printf( STDERR "Done searching exhaustively. Time to conduct exhaustive search: %10.3g\n", $t2-$t1point5); 
  # ############## end of search  ##############

  if ($output_graph) {
    my $graph_out_filename = $input_filename_stem . '.graph';
    open my $fhout, ">", $graph_out_filename or die "Couldn't open $graph_out_filename for writing.\n";
    print $fhout "# n_keep: $n_nearest_to_keep   n_near: $n_nearest_for_search   n_extras: $n_extras \n";
    print $fhout $genotype_graph->as_string(0);
    close $fhout;
 }
  if ($output_distance_matrix) {
     my $distance_out_filename = $input_filename_stem . '.dmatrix';
     open my $fhout, ">", $distance_out_filename or die "Couldn't open $distance_out_filename for writing.\n";
     print $fhout $genotype_graph->distance_matrix_as_string($multiplier);
     close $fhout;
  }
  my $t3 = gettimeofday();
  my ($tex, $tgr) = (-1, -1);
  if (defined $t1point5) {
     $tgr = $t1point5 - $t1;
     $tex = $t2 - $t1point5;
  }
  printf("times: construct: %12.3f  gsearch: %12.3f  exsearch: %12.3f  output: %12.3f  total: %12.3f\n", $t1-$t0, $tgr, $tex, $t3-$t2, $t3-$t0);
}                               # end main


sub file_to_string{
  my $filename = shift;
  my $input_string = '';
  if (-f $filename) {
    open my $fhin, "<", $filename or die "open $filename for reading failed.\n";
    while (my $line = <$fhin>) {
      #   print STDERR $line;
      #     next if($line =~ /^\s*#/);
      $input_string .= $line;
      #  print STDERR $fasta_as_read;
    }
    close $fhin;
  } else {
    die "file $filename does not exist.\n";
  }
  return $input_string;
}
