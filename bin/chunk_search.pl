#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw (min max sum);
use Time::HiRes qw( gettimeofday );
use Inline 'C';
use constant MISSING_DATA => 'X';
use constant BIG_NUMBER => 1_000_000_000;

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

#use GenotypeGraph;
#use GenotypeGraphNode;
use Genotype;
use ChunkSet;
use Getopt::Long;
use lib '/home/tomfy/Orthologger/lib/';
use TomfyMisc qw ' fasta2seqon1line ';

{                               ###########

  # my @xs = (0..100);
  # my $randxs = randomize_array(\@xs);
  # print "randomized ints: ", join(' ', @$randxs), "\n";
  # exit;
  
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
  my $chunk_size = 6;


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
	     'chunk_size=i' => \$chunk_size,
	    );

  die "Input file for constructing graph must be specified, is undefined. Bye.\n" if(!defined $input_filename);

  if ($seed > 0) {
    srand($seed);
  } else {
    srand();
  }

  # my %id_sequence = (); # hash to hold id:sequence pairs
  my $input_filename_stem = $input_filename;
  $input_filename_stem =~ s/\.(fasta|graph|dmatrix)$//; # remove the part after the last '.'
  my $fasta1_filename = $input_filename_stem . '.fasta'; # get an array or hash of Genotype objs from this fasta file. 
  my $fasta1_string = (-f $fasta1_filename)?
    TomfyMisc::fasta2seqon1line(file_to_string($fasta1_filename))
      :  print STDERR "# Specified fasta file: $fasta1_filename,  does not exist.\n",
      "Will construct graph from .graph or .dmatrix file\n",
      "No sequence info; no search possible.\n";;


  ############## construct data structure  ##################
  my $t0 = gettimeofday();
  my @chunk_hashes = ();
  my $id_sequence = fasta_string_to_hash($fasta1_string);
  my $ids_stored = 0;
  while ( my ($id, $seq) = each %$id_sequence) {
    #    print STDERR "$id  $seq \n";
    #   my %chunkseq_ids = ();	# keys: chunk_sequences,
    my $chunk_index = 0;
    while (length $seq >= $chunk_size) {
      my $chunk_sequence = substr($seq, 0, $chunk_size, '');
      #     print STDERR "  $chunk_index  $chunk_sequence\n";
      #   print "AAA:  ", (defined $chunk_hashes[$chunk_index])? ref $chunk_hashes[$chunk_index] : 'undef', "\n";
      $chunk_hashes[$chunk_index] //= {};
      if (defined $chunk_hashes[$chunk_index]->{$chunk_sequence}) {
	push @{$chunk_hashes[$chunk_index]->{$chunk_sequence}}, $id;
	#	print STDERR "XXX: ", join(' ', @{$chunk_hashes[$chunk_index]->{$chunk_sequence}}), "\n";
      } else {
	$chunk_hashes[$chunk_index]->{$chunk_sequence} = [$id];
	#	print STDERR "YYY: $chunk_index  $chunk_sequence  ids: ", join(' ', @{$chunk_hashes[$chunk_index]->{$chunk_sequence}}), "\n";
      }
      $chunk_index++;
    }
    $ids_stored++;
    print STDERR "# $ids_stored \n" if($ids_stored % 100  == 0);
  }
    my $t1 = gettimeofday();
  print STDERR "# done storing sequences \n";
  # while(my ($chi, $seqids) = each @chunk_hashes){
  #   print "# chunk index:  $chi \n";
  #   while(my ($s, $ids) = each %$seqids){
  #     print " seq, n ids:   $s  ", scalar @$ids, "\n"; # join(',', @$ids), "\n";
  #   }
  # }
#  exit;
  
  # read in another fasta file and search for each of its sequences ...
  my $t2 = gettimeofday();
  my $t3;
  if (defined $other_fasta  and  -f $other_fasta) {
    my $other_fasta_string = TomfyMisc::fasta2seqon1line(file_to_string($other_fasta));
    my $other_id_sequence = fasta_string_to_hash($other_fasta_string);
    my %count_distrib = ();
    my $n_searches_done = 0;
    while (my  ($other_id, $other_seq) = each %$other_id_sequence) {
      my %matchid_count = (); # keys: ids in db; values: counts of chunks matching search seq.
      my $chunk_index = 0;
      my $chunk_start = 0;
      my $other_seq_length = length $other_seq;
      while ($chunk_start <= $other_seq_length - $chunk_size) {
	my $other_chunk_sequence = substr($other_seq, $chunk_start, $chunk_size);
	my $the_matching_ids = $chunk_hashes[$chunk_index]->{$other_chunk_sequence};
	for (@$the_matching_ids) {
	  $matchid_count{$_}++;
	}
	$chunk_index++;
	$chunk_start += $chunk_size;
      }
      for (values %matchid_count) {
	$count_distrib{$_}++;
      }
      print "$other_id     ";
      my @best_matchids = sort {$matchid_count{$b} <=> $matchid_count{$a}} keys %matchid_count;
      for (@best_matchids[0..24]) {
#	print "A: ", $id_sequence->{$_}, "\n";
#	print "B: $other_seq\n";
	printf("%4d %4d %8.5f   ", $_,  $matchid_count{$_}, distance_C($id_sequence->{$_}, $other_seq));
      }
      print "\n";
      $n_searches_done++;
      print STDERR "# n searches done: $n_searches_done \n" if($n_searches_done % 100  == 0);
    }
    $t3 = gettimeofday();
    for (0..100) {
      print STDERR "$_  ", $count_distrib{$_} // 0 , "\n";
    }
  }

  print "# construct data structure: ", $t1-$t0, "   search: ", $t3-$t2, "\n";

}				# end main


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

sub fasta_string_to_hash{
  my $fasta_sequence = shift;
  my $id_seq = {};

  my @lines = split("\n", $fasta_sequence);
  while (@lines) {
    my $idline = shift @lines;
    next if($idline =~ /^\s*#/);
    if ($idline =~ /^>(\S+)/) {
      my $id = $1;
      my $seq = shift @lines;
      $seq =~ s/\s+$//;
      $id_seq->{$id} = $seq;
    }
  }
  return $id_seq;
}

sub scrambled_chunks{
  my $sequence = shift;
  my $chunk_size = shift;
  my $indices = shift;	# array ref of indices: [13,2,45,5,17,6, ...] 

  my @chunk_seqs = ();
  my @chars = split('', $sequence);
  while (scalar @$indices >= $chunk_size) {
    my @chunk_indices = splice( @$indices, 0, $chunk_size);
    my $chunk_sequence = join('', @chars[@chunk_indices]);
    push @chunk_seqs, $chunk_sequence;
  }
  return \@chunk_seqs;
}

sub randomize_array{
  my $xs = shift;
  my $n_to_do = scalar @$xs;
  my $jswap = -1; # where randomly chosen elem will be swapped to.
  while($n_to_do > 1){
    my $rand_index = int(rand($n_to_do));
    my $tmp = $xs->[$jswap];
    $xs->[$jswap] = $xs->[$rand_index];
    $xs->[$rand_index] = $tmp;
    $n_to_do--;
    $jswap--;
  }
  return $xs;
}



# sub distance{ # calculate distance between this genotype obj. and another
#    my $this_gt = shift;    # string
#    my $other_gt = shift;   # string

#    my $Dcalc = 1;

#    if ($Dcalc == 1) {		       # use inline::C function: (much faster)
#       my ($d, $agmr, $hgmr, $ogmr);
#       distances_C(0, BIG_NUMBER, $this_gt, $other_gt, $d, $agmr, $hgmr, $ogmr);
# #      printf("%10.6f %10.6f, %10.6f, %10.6f \n", $d, $agmr, $hgmr, $ogmr);
#       return ($d, $agmr, $hgmr, $ogmr);
#     }elsif($Dcalc == 2){
#       my ($d, $agmr, $hgmr, $ogmr) = (distance_C($this_gt, $other_gt), 5, 5, 5);
#       return ($d, $agmr, $hgmr, $ogmr);
#    } else {			# pure perl. much slower.
#       my $distance = 0;
#       my $count_both = 0; # count of snps with data present in both sequences
#       my $count_missing = 0; # count of snps with data absent in one or both sequences
#       if (length $this_gt == length $other_gt) {
#          for my $i (0 .. (length $this_gt) - 1) {
#             my $c1 = substr($this_gt, $i, 1);
#             my $c2 = substr($other_gt, $i, 1); # $other_gt->[$i];
#             if ( ($c1 eq MISSING_DATA)  or ($c2 eq MISSING_DATA) ) {
#                $count_missing++;
#             } else {
#                $count_both++;
#                if ($c2 != $c1) {
#                   if ($c1 == 0) {
#                      $distance += $c2 - $c1;
#                   } elsif ($c1 == 1) {
#                      $distance += 1;
#                   } elsif ($c1 == 2) {
#                      $distance += $c1 - $c2;
#                   } else {
#                      die "c1 has unhandled value: $c1 \n";
#                   }
#                }
#             }
#          }
#       }
#       return ($count_both > 0)? $distance/$count_both : 1000;
#    }
# }




__DATA__
############################################
########### inline C stuff #################
__C__

  /* distance between two genotype strings  */
/* 0-1, 1-2 -> d=1; 0-2 -> d=2  */
double distance_C(char* str1, char* str2) {
  int i = 0;
  char c1;
  char c2;
  int dist = 0;
  int count = 0;
  while(c1 = str1[i]) {
    c2 = str2[i];

    if (c1 == '0') {
      if(c2 == '0'){
        count++;
      }else if(c2 == '1'){
        dist++;
        count++;
      }else if(c2 == '2'){
        dist += 2;
        count++;
      } // else c2 is missing data
} else if(c1 == '1'){

  if (c2 == '0') {
  dist++;
count++;
} else if(c2 == '1'){
          count++;
       } else if(c2 == '2'){
                 dist += 1;
                 count++;
              }

} else if(c1 == '2'){

   if (c2 == '0') {
       dist += 2;
       count++;
    } else if(c2 == '1'){
              dist++;
              count++;
           } else if(c2 == '2'){
                     count++;
                  }

}

i++;
}
//  printf("%d  %d \n", dist, count);
if (count > 0) {
   //  printf("%g\n", 1.0*dist/count);
   return 1.0*dist/count;
} else
  return 1000;
}



double distances_C(int start, int n_homozyg, char* str1, char* str2, SV* dist, SV* agmr, SV* hgmr, SV* ogmr) {
     // calculates 3 kinds of distance: dist, agmr and hgmr
     // start at snp with index start, go until end or until n_homozyg pairs with both homozygous.
       int i = start;
     char c1;
     char c2;
    
       int count_02 = 0; // 0<->2 (both directions)
         int count_0022 = 0; // 0<->0 or 2<->2
           int count_11 = 0; // 1<->1
             int count_1x = 0; // 1<->0 or 1<->2
       int hcount = 0; // count only snps homozygous in both
         while (c1 = str1[i]) {
            c2 = str2[i];

            if (c1 == '0') {

               if (c2 == '0') {
                   count_0022++;
                } else if(c2 == '1'){
                          count_1x++;
                       } else if(c2 == '2'){
                                 count_02++;
                              }
               // else c2 is missing data
            } else if(c1 == '1'){

               if (c2 == '0') {
                  count_1x++;
                } else if(c2 == '1'){
                          count_11++;
                       } else if(c2 == '2'){
                                count_1x++;
                              }

            } else if(c1 == '2'){

               if (c2 == '0') {
                  count_02++;
                } else if(c2 == '1'){
                          count_1x++;
                       } else if(c2 == '2'){
                                 count_0022++;
                              }
            }
            i++;
//	     printf("%d %d %d \n", start, n_homozyg, count_0022 + count_02);
            if((count_0022 + count_02) >= n_homozyg){ break; }  // stop after n_homozyg homozygous-only pairs
         } // end of while

         double homozyg_denom = count_0022 + count_02; // both genotypes of pair homozygous
         double other_denom = count_11 + count_1x; // 0 or 1 genotypes of pair are homozygous ( no missing data)
         double all_denom = homozyg_denom + other_denom; // all except if there is missing data

//	 printf("%8d %8d %8d %8d   %8.1f\n", count_0022, count_02, count_11, count_1x, all_denom);
     double homozyg_dist = (homozyg_denom > 0)? count_02/homozyg_denom : 1000;
     double other_dist = (other_denom > 0)? count_1x/other_denom : 1000;
     double all_dist = (all_denom > 0)? (count_1x + count_02)/all_denom : 1000;
     double d_dist = (all_denom > 0)? (count_1x + 2*count_02)/all_denom : 1000;

     sv_setnv(hgmr, homozyg_dist);
     sv_setnv(ogmr, other_dist);
     sv_setnv(agmr, all_dist);
     sv_setnv(dist, d_dist);

  }


/* number of mismatches between to strings up to some max  */
  /* returns 0,1,2,...,max_mismatches or -1 if there are more mismatches */
void count_mismatches_up_to_limit_C(char* str1, char* str2, int max_mismatches, SV* n_ok, SV* n_mm) {
  int i = 0;
  char c1;
  char c2;
  int mismatch_count = 0;
  while(c1 = str1[i]) {
     c2 = str2[i];
     if (c1 != c2) {
       mismatch_count++;
     }
     if(mismatch_count > max_mismatches){
       mismatch_count = -1;
       break;
     }
     i++;
/* i is now the number of characters which have been compared without exceeding the mismatch limit. */
  }
 /* Inline_Stack_Vars;
Inline_Stack_Reset;
Inline_Stack_Push(sv_2mortal(newSVuv(i)));
Inline_Stack_Push(sv_2mortal(newSVuv(mismatch_count)));
Inline_Stack_Done; */

  sv_setiv(n_ok, i);
  sv_setiv(n_mm, mismatch_count);

}

