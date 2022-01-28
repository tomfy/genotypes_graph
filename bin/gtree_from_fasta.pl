#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw (min max sum shuffle);
use Getopt::Long;
use Time::HiRes qw( gettimeofday );
use lib '/home/tomfy/Orthologger/lib/';
# use TomfyMisc qw ' fasta2seqon1line ';
use constant MISSING_DATA => 'X';
use Math::GSL::RNG  qw( :all );
use Inline C => '../lib/inlinec.c';

no warnings 'recursion';

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

use GenotypeTree;
use GenotypeTreeNode;
use Genotype;

{                               ###########
  my $input_filename = undef;	# input fasta file name.
  my $newick_out = 0;
  my $rng_type = $gsl_rng_mt19937;
  my $rng_seed = undef;
  my $p_missing = 0; # replace genotypes with MISSING_DATA with this probability
  my $query_fasta = undef;	# query genotype sets fasta file.
  my $max_mismatches = 1;
  my $error_prob = 0;	    # add errors to data with this probability
  my $chunk_size = 100000000;
  my $n_chunks = 1;
  my $sort_markers = 0;
  my $min_usable_pairs_to_store = 2;
  my $check_prob = 0; # do 'exhaustive' (i.e. treeless) calculation for this fraction of accession pairs.
  my $gmr_prob = 0; # get agmr, hgmr for this fraction of pairs. for comparison with tree-based fom.

  GetOptions(
	     'fasta1|f1=s' => \$input_filename,
	     'fasta2|f2=s' => \$query_fasta, # queries
	     'sort_markers!' => \$sort_markers,
	     'chunk_size=i' => \$chunk_size,
	     'n_chunks=i' => \$n_chunks,
	     'max_mismatch=i' => \$max_mismatches,
	     'min_to_store=i' => \$min_usable_pairs_to_store,
	     'check_prob=f' => \$check_prob,
	     'gmr_prob=f' => \$gmr_prob,

	     'newick_out!' => \$newick_out,
	     'rng_type=s' => \$rng_type,
	     'seed=i' => \$rng_seed,
	     'p_missing=f' => \$p_missing,
	     'error_prob=f' => \$error_prob,
	    );

  # ############## set up the random number generator ###############################
  if ($rng_type eq 'sys') {
    $rng_type = $gsl_rng_default;
  } elsif ($rng_type eq 'mt') {
    $rng_type = $gsl_rng_mt19937;
  } elsif ($rng_type eq 'lux') {
    $rng_type = $gsl_rng_ranlxd2;
  }
  my $the_rng = (defined $rng_seed)? Math::GSL::RNG->new($rng_type, $rng_seed) : Math::GSL::RNG->new($rng_type) ; # RNG
  # ################ RNG is set up ##################################################

  print "# genotype input files: $input_filename  $query_fasta \n";
  print "# p_missing: $p_missing   error_prob: $error_prob   RNG type: $rng_type\n";
  print "# max_mismatch: $max_mismatches \n";
  print "# chunk_size: ", (defined $chunk_size)? $chunk_size : 'undef',   "  n_chunks $n_chunks \n";

  my ($t_start, $t_readin, $t_exhaustive, $t_compare);

  ### #######  read in genotypes  ###############
  $t_start = gettimeofday();

  my $input_filename_stem = $input_filename;
  $input_filename_stem =~ s/\.\S+$//; # remove the part after the (last) '.'

  my ($sequence_length, $id_gtsobj1, $mindex_okcount1) = read_genotypes_from_fasta($input_filename, $the_rng, $p_missing, $error_prob);

  my ($query_sequence_length, $qid_gtsobj, $qmindex_okcount) = read_genotypes_from_fasta($query_fasta, $the_rng, $p_missing, $error_prob);

  my ($n_db, $n_query) = (scalar keys %$id_gtsobj1, scalar keys %$qid_gtsobj);
  my $n_accessions1 = scalar keys %$id_gtsobj1;
  my $n_accessions2 = scalar keys %$qid_gtsobj;
  my $mindex_okcount = sum_hashes($mindex_okcount1, $qmindex_okcount);

  my @array_0123etc = ();

  for (0..$sequence_length-1) {
    $array_0123etc[$_] = $_;
  }
  @array_0123etc = shuffle(@array_0123etc);
  my @sorted_match_indices = ($sort_markers)?
    sort {$mindex_okcount->{$b} <=> $mindex_okcount->{$a} }
    keys %$mindex_okcount	# sort by okcount (high to low)
    :
    @array_0123etc;

  my $block_size = 1000;
  for (my $i = 0; $i < scalar @sorted_match_indices; $i += $block_size) {
    my $iend = min($i+$block_size-1, $#sorted_match_indices);
    print STDERR "block $i; i, iend: $i $iend \n";
    @sorted_match_indices[$i..$iend] = shuffle(@sorted_match_indices[$i..$iend]);
  }
  my @chunk_index_arrayrefs = ();

  my $n_chunks_used = $n_chunks;
  my $i_start = 0;
  for my $i_chunk (1..$n_chunks) {
    my $i_end = $i_start + $chunk_size - 1;
    # print STDERR "$i_end  ", scalar @sorted_match_indices, "\n";
    if ($i_end >= $#sorted_match_indices) {
      $i_end = $#sorted_match_indices;
      my @chunk_indices = @sorted_match_indices[$i_start .. $i_end];
      push @chunk_index_arrayrefs, \@chunk_indices;
      $n_chunks_used = $i_chunk;
      last;
    } else {
      my @chunk_indices = @sorted_match_indices[$i_start .. $i_end];
      push @chunk_index_arrayrefs, \@chunk_indices;
    }
    $i_start += $chunk_size;
  }
  print "# n_chunks_used: $n_chunks_used \n";
  $t_readin = gettimeofday() - $t_start;
  # ################ done reading in genotypes and constructing genotype objects #########

  #  print STDERR "time to read input, construct genotype objects: ", $t1-$t0, " sec.\n";

  ### do exhaustively for fraction $check_prob (if > 0). For comparison: ###
  my $n_compare_exh = 0;
  $t_start = gettimeofday();
  my %Ex_qid_matchidinfosum = ();
  if ($check_prob > 0) {
    while (my ($qid, $q_gtsobj) = each %$qid_gtsobj) {
      my $ex_mid_matchinfosum = {};
      my $mid_matchinfo = {};
      (my $count, $mid_matchinfo) = exhaustive_search_1query($id_gtsobj1, $q_gtsobj, \@chunk_index_arrayrefs, $max_mismatches, $the_rng, $check_prob, $ex_mid_matchinfosum, $min_usable_pairs_to_store);
      $n_compare_exh += $count;
      $Ex_qid_matchidinfosum{$qid} = $ex_mid_matchinfosum;
    }
  }
  $t_exhaustive = gettimeofday() - $t_start;
  #### end of exhaustive search ####

 
  ########### for each chunk, construct tree, search tree to find matches to queries,  #########

  my ($t_treeconstruct, $t_treesearch, $t_check) = (0, 0, 0);

  #############  construct tree for all chunks ##############
  $t_start = gettimeofday();
  my @trees = ();
  while (my($i_chunk, $chunk_indices) = each @chunk_index_arrayrefs) {
    my $gtree = GenotypeTree->new( { depth => scalar @$chunk_indices } );
    while (my ($qid, $gobj1) = each %$id_gtsobj1) {
      $gtree->add_genotype($gobj1, $chunk_indices);
    }
    print $gtree->as_newick(), "\n\n" if($newick_out);
    push @trees, $gtree;
    print STDERR "# chunk $i_chunk tree constructed.\n" if($i_chunk%10 == 0);
  }
  print STDERR "# trees for all $n_chunks_used chunks constructed.\n";
  $t_treeconstruct = gettimeofday - $t_start;
  ############### trees constructed ###############

  ######## Do tree searches ########################
  $t_start = gettimeofday();
  #my %tree_qid_results = ();
  my %Tree_qid_matchidinfosum = (); 
  my $n_compare_tree = 0;
#  print "XXXXXXX\n";
  while (my ($qid, $q_gtsobj) = each %$qid_gtsobj) {
      my $tr_mid_matchinfosum = {};
      # print "qid: $qid \n";
      while (my ($i_chunk, $chunk_indices) = each @chunk_index_arrayrefs) {
	  #print "$i_chunk \n";
	########## search in tree for matches to sequences from $query_fasta file ####
      my $gtree = $trees[$i_chunk];
      my %qid_matchidinfo = ();
      $gtree->search($q_gtsobj, $chunk_indices, $max_mismatches, $min_usable_pairs_to_store, $tr_mid_matchinfosum);
	  $n_compare_tree += $n_db;
	  #print "ZZZZ:  [", scalar keys %$tr_mid_matchinfosum, "]\n";
    }

#    while(my($k, $v) = each %{$tr_mid_matchinfosum}){
#	print STDERR "qid: $qid mid: $k   ", join(", ", @$v), "\n";
#    }
    
    $Tree_qid_matchidinfosum{$qid} = $tr_mid_matchinfosum;
  }
  $t_treesearch += gettimeofday() - $t_start;
  ############ done with tree search for matches to this query  ######

  ######## compare tree results with exhaustive results: #####################
  $t_start = gettimeofday();
  my %ngoodbad_count = ();
  my ($total_trex_agree_count, $total_trex_disagree_count) = (0, 0); # number of queries for which exhaustive and tree searches disagree.
  if ($check_prob > 0) {
    while (my ($qid, $q_gtobj) = each %$qid_gtsobj) {
      my ($trex_agree_count, $trex_disagree_count) = (0, 0);
      while (my($mid, $tv) = each %{$Tree_qid_matchidinfosum{$qid}}) {
	my $xv = $Ex_qid_matchidinfosum{$qid}->{$mid} // undef;
	next if(!defined $xv); # since may only be doing exhaustive for a subset, only check those with exhaustive result present.
	if (join(";", @$tv) eq join(";", @$xv)) {
	  $trex_agree_count++;
	} else {
	  print "match id: $mid  exhaustive: ", join(" ", @$xv), "   tree: ", join(" ", @$tv), " \n";
	  $trex_disagree_count++;
	}
      }
      $ngoodbad_count{"$trex_agree_count    $trex_disagree_count"}++;
      $total_trex_disagree_count += $trex_disagree_count;
      $total_trex_agree_count += $trex_agree_count;
    }
  }

  if ($total_trex_disagree_count > 0) {
    my $trex_disagree_string = sprintf("# all chunks  ngood nbad  count \n");
    while (my($ngb, $count) = each %ngoodbad_count) {
      my ($n_agree, $n_disagree) = split(" ", $ngb);
      $trex_disagree_string .= sprintf("#       %12s  %3i \n",  $ngb, $count) if($n_disagree > 0);
    }
    print $trex_disagree_string;
  }
  $t_check = gettimeofday - $t_start;

  
  # report results
  my %quad_count = ();	       #key: '01' <-> est_agmr < a0, agmr > a0
  my $n_agmrs = 0;
  #$t_start = gettimeofday();

  my $t_agmr = 0;
  while (my($qid, $mid_sums) = each %Tree_qid_matchidinfosum) {
    my %mid_estagmr = ();
    my $q_sequence = $qid_gtsobj->{$qid}->sequence();
    while (my($mid, $sums) = each %$mid_sums) {
      my $est_agmr = est_agmr($n_chunks_used, $max_mismatches, $min_usable_pairs_to_store, @$sums); #($denom > 0)? $total_mms/$denom : 2;
      $mid_estagmr{$mid} = $est_agmr;
      my $est_agmr_a = 0; #est_agmr_a($n_chunks_used, $max_mismatches, $min_usable_pairs_to_store, @$sums); #($denom > 0)? $total_mms/$denom : 2;
    }
     $t_start = gettimeofday();
    while (my($mid, $sums) = each %$mid_sums) {
     my $est_agmr = $mid_estagmr{$mid};
      next if($est_agmr > 0.5);
      next if(gsl_rng_uniform($the_rng->raw()) >= $gmr_prob);

      my $m_sequence = $id_gtsobj1->{$mid}->sequence();
      my ($agmr_n, $agmr_d, $hgmr_n, $hgmr_d) = (0, 0, 0, 0);
      #    agmr_hgmr_C($q_sequence, $m_sequence, $agmr_n, $agmr_d, $hgmr_n, $hgmr_d);
      agmr_C($q_sequence, $m_sequence, $agmr_n, $agmr_d);
      $n_agmrs++;
      my $agmr = ($agmr_d > 0)? $agmr_n/$agmr_d : -1;
      # agmr_perl($q_sequence, $m_sequence);

      my $hgmr = ($hgmr_d > 0)? $hgmr_n/$hgmr_d : -1;

      my $quad = ($est_agmr < 0.15)? '0' : '1';
      $quad .= ($agmr < 0.15)? '0' : '1';
      $quad_count{$quad}++;

      print  "$qid  $mid   ", join("  ", @$sums), "   $est_agmr  $agmr  $hgmr \n";
    }
    $t_agmr += gettimeofday() - $t_start;
      # 	my @sorted_matchids = sort {$mid_estagmr{$a} <=> $mid_estagmr{$b}} keys %mid_estagmr;
      # print STDERR "$qid    ";
      # for (@sorted_matchids[0..4]){
      # 	print STDERR "$_ ", $mid_estagmr{$_}, "  ";
      # }print STDERR "\n";
    
  }
 # $t_compare = gettimeofday() - $t_start;
  while (my($k, $v) = each %quad_count) {
    print STDERR "$k  $v \n";
  }
  print "n_compare_exh: $n_compare_exh  \n";
  print "n_compare_tree: $n_compare_tree  n_chunks_used: $n_chunks_used \n";
  my $timing_string = '';
  $timing_string .= sprintf("#      input   n gobjs: db %10i  query: %5i   time %10.3f sec.\n", $n_db, $n_query, $t_readin);
  $timing_string .= sprintf("#   treeless   n compare:  %10i                 time %10.3f  %10.4f\n", $n_compare_exh, $t_exhaustive, ($n_compare_exh>0)? 1000*$t_exhaustive/($n_compare_exh/$n_chunks_used) : -1);
  $timing_string .= sprintf("# treeconstr   n constr:   %10i                 time %10.3f \n", scalar @trees, $t_treeconstruct);
  $timing_string .= sprintf("# treesearch   n compare:  %10i                 time %10.3f  %10.4f\n", $n_compare_tree, $t_treesearch, 1000*$t_treesearch/($n_compare_tree/$n_chunks_used));
  $timing_string .= sprintf("#  true agmr   n agmrs:    %10i                 time %10.3f  %10.4f\n", $n_agmrs, $t_agmr, 1000*(($n_agmrs > 0)? $t_agmr/$n_agmrs : -1.0));
  $timing_string .= sprintf("#      check   %5i of %5i agree.\n",  $total_trex_agree_count, $total_trex_agree_count+$total_trex_disagree_count);
  print $timing_string, "\n";
}				# end main


####  ############### subroutines ################## ####
sub est_agmr{
  my $n_chunks_used = shift;
  my $max_mismatches = shift;
  my $min_usable_pairs_to_store = shift;
  my $total_stored_mms = shift;
  my $n_usable_compared = shift;
  my $n_stored_chunks = shift;
  my $n_unstored_chunks = $n_chunks_used - $n_stored_chunks;
  my $total_mms = $total_stored_mms + $n_unstored_chunks*$max_mismatches;
  my $denom = $n_usable_compared + 0.5*$n_unstored_chunks*($min_usable_pairs_to_store); #-1);
  my $est_agmr = ($denom > 0)? $total_mms/$denom : 2;
  return $est_agmr;
}

sub est_agmr_a{
  my $n_chunks_used = shift;	# L
  my $max_mismatches = shift;
  my $min_usable_pairs_to_store = shift; # k0
  my $total_stored_mms = shift;
  my $n_usable_compared = shift;
  my $n_stored_chunks = shift;				     # L>
  my $n_unstored_chunks = $n_chunks_used - $n_stored_chunks; # L<
  my $total_mms = $total_stored_mms + $n_unstored_chunks*$max_mismatches;
  my ($L, $k0, $Lge) = ($n_chunks_used, $min_usable_pairs_to_store, $n_stored_chunks);
  my $denom = $n_usable_compared;
  if ($L > $Lge) {
    my $xi = $Lge/$L;
    $denom += ($L-$Lge)*(1/(1 - $xi**($k0-1)) - ($k0-1)*$xi/(1 - $xi));
  }
  my $est_agmr = ($denom > 0)? $total_mms/$denom : 2;
  return $est_agmr;
}

sub notX_indices{
  my $sequence = shift;
  my @notxindices = ();
  for (my $i=0; $i< length $sequence; $i++) {
    push @notxindices, $i if(substr($sequence, $i, 1) ne MISSING_DATA);
  }
  return \@notxindices;
}

sub neitherXnor1_indices{
  my $sequence = shift;
  my @notx1indices = ();
  for (my $i=0; $i< length $sequence; $i++) {
    my $gt = substr($sequence, $i, 1);
    next if($gt eq '1');
    push @notx1indices, $i if($gt ne MISSING_DATA);
  }
  return \@notx1indices;
}

sub lose_data{ # with prob. $p_missing replace genotypes with MISSING_DATA
  my $sequence = shift;
  my $rng = shift;
  my $p_missing = shift;
  if ($p_missing > 0) {
    my $L = length $sequence;
    for my $i (0..$L-1) {
      if (gsl_rng_uniform($rng->raw()) < $p_missing) {
	substr($sequence, $i, 1, MISSING_DATA);
      }
    }
  }
  return $sequence;
}

sub add_errors_to_data{ # with prob. $error_prob replace genotypes with other (wrong) gts.
  my $sequence = shift;
  my $rng = shift;
  my $p_error = shift;
  if ($p_error > 0) {
    my $pesqrd = $p_error*$p_error;
    my $L = length $sequence;
    for my $i (0..$L-1) {
      my $random = gsl_rng_uniform($rng->raw());
      if ($random < $p_error) {
	my $gt = substr($sequence, $i, 1);
	if ($gt eq '0') {
	  if ($random < $pesqrd) {
	    $gt = '2';
	  } else {
	    $gt = '1';
	  }
	} elsif ($gt eq '2') {
	  if ($random < $pesqrd) {
	    $gt = '0';
	  } else {
	    $gt = '1';
	  }
	} elsif ($gt eq '1') {
	  if ($random < 0.5*$p_error) {
	    $gt = '0';
	  } else {
	    $gt = '2';
	  }
	}
	substr($sequence, $i, 1, $gt);
      }
    }
  }
  return $sequence;
}


sub count_mismatches{
  my $gobj1 = shift;
  my $gobj2 = shift;
  my $gs1 = $gobj1->sequence();
  my $gs2 = $gobj2->sequence();
  my $L1 = length $gs1;
  die if(length $gs2 != $L1);
  my $n_mismatches = 0;
  for my $i (0..$L1-1) {
    my ($c1, $c2) = ( substr($gs1, $i, 1), substr($gs2, $i, 1) );
    if (($c1 ne MISSING_DATA)  and  ($c2 ne MISSING_DATA)  and  ($c1 ne $c2)) {
      $n_mismatches++;
    }
  }
  return $n_mismatches;
}

sub count_mismatches_up_to_limit_perl{
  my $max_mismatches = shift;
  my $gobj1 = shift;
  my $gobj2 = shift;
  my $gs1 = $gobj1->sequence();
  my $gs2 = $gobj2->sequence();
  my $L1 = length $gs1;
  die if(length $gs2 != $L1);
  my $n_mismatches = 0;
  for my $i (0..$L1-1) {
    my ($c1, $c2) = ( substr($gs1, $i, 1), substr($gs2, $i, 1) );
    if (($c1 ne MISSING_DATA)  and  ($c2 ne MISSING_DATA)  and  ($c1 ne $c2)) {
      $n_mismatches++;
    }
    if ($n_mismatches > $max_mismatches) {
      $n_mismatches = -1;
      last;
    }
  }
  return $n_mismatches;
}

sub read_genotypes_from_fasta{ # read fasta with genotype sequences. return an array ref of genotype objects.
  my $input_filename = shift;
  my $the_rng = shift;
  my $p_missing = shift // 0;
  my $error_prob = shift // 0;
  my $input_string = '';
  if (-f $input_filename) {
    open my $fhin, "<", $input_filename or die "open $input_filename for reading failed.\n";
    while (my $line = <$fhin>) {
      next if($line =~ /^\s*#/);
      $input_string .= $line;
    }
    close $fhin;
  } else {
    die "file $input_filename does not exist.\n";
  }

  #my @marker_ok_counts = ();
  my %mindex_okcount = ();

  #  print STDERR "input file name: [$input_filename]\n";
  my %count_errors = ();
  my %id_gtsobj = ();
  my $sequence_length = undef;
  # my @fasta_lines = split ("\n", TomfyMisc::fasta2seqon1line($input_string));
  my @fasta_lines = split("\n", $input_string);
  while (@fasta_lines) {
    my $line = shift @fasta_lines;
    next if($line =~ /^\s*#/);	# skip comment lines
    #  if ($line =~ /^>(\S+)\s+(\S+)\s+(\S+)/) {
    my ($id, $generation, $pedigree);
    $id = $1 if($line =~ /^>(\S+)/);
    $generation= $1 if($line =~ /^>\S+\s+(\S+)/);
    $pedigree = $1 if($line =~ /^>\S+\s+\S+\s+(\S+)/);
    my $sequence = shift @fasta_lines;
    #   print "the sequence $sequence \n";
    $sequence =~ s/\s+//g;
    $sequence_length = length $sequence if(!defined $sequence_length);
    die "Sequence lengths must all be the same.\n" if(length $sequence != $sequence_length);
    #   print STDERR "before: $sequence\n";
    my $sequence_w_errors = add_errors_to_data($sequence, $the_rng, $error_prob);
    if ($error_prob > 0) {
      for (my $i = 0; $i < length $sequence; $i++) {
	my $key = substr($sequence, $i, 1) . substr($sequence_w_errors, $i, 1);
	$count_errors{$key}++;
      }
    }
    $sequence = lose_data($sequence_w_errors, $the_rng, $p_missing);

    for (my $i = 0; $i < length $sequence; $i++) {
      $mindex_okcount{$i}++ if(substr($sequence, $i, 1) ne 'X');
    }

    #  print STDERR "after:  $sequence\n";
    $id_gtsobj{$id} = Genotype->new(
				    { id => $id,
				      generation => $generation,
				      pedigree => $pedigree,
				      sequence => $sequence }
				   );
    # }
  }
  # while (my($k, $v) = each %count_errors) {
  #   print STDERR "$k  $v \n";
  # }
  #   my @sorted_marker_indices = sort {$mindex_okcount{$b} <=> $mindex_okcount{$a}} keys %mindex_okcount;
  # for my $j (0..49){ # (@sorted_marker_indices){
  #   my $i = $sorted_marker_indices[$j];
  #     print STDERR "i, okcount: $i ", $mindex_okcount{$i}, "\n";
  #   }
  return ($sequence_length, \%id_gtsobj, \%mindex_okcount);
}

sub sum_hashes{
  my $hr1 = shift;
  my $hr2 = shift;

  while (my ($k, $v) = each %$hr2) {
    $hr1->{$k} += $v;
  }
  return $hr1;
}

sub exhaustive_search_1query{
  my $db_gts_objs = shift;	# array ref of Genotypes objects
  my $q_gts_obj = shift;
  my $chunk_index_arrayrefs = shift;
  my $max_mismatches = shift;
  my $rng = shift;
  my $prob = shift;
  my $mid_matchinfosum = shift;
  my $min_to_store = shift;
 
  my $count = 0;
  my %mid_matchinfo = ();
  while (my($i, $db_gts_obj) = each %$db_gts_objs) {
    # my $sequence1 = $db_gt_obj->get_chunk($index_arrayref);
    my $mid =  $db_gts_obj->id();
    if (gsl_rng_uniform($rng->raw()) < $prob) {
      for my $index_arrayref (@$chunk_index_arrayrefs) {
	my $q_sequence = $q_gts_obj->get_chunk($index_arrayref); # $chunk_index_arrayrefs[0]);
	my ($n_pairs_compared, $n_usable_pairs_compared, $mismatch_count) = (0, 0, 0);
	#GenotypeTreeNode::
	count_mismatches_C(
			   $db_gts_obj->get_chunk($index_arrayref),
			   $q_sequence,
			   # $sequence1,
			   # $sequence2,
			   $max_mismatches, 0,
			   $n_pairs_compared, $n_usable_pairs_compared, $mismatch_count);

	if ($n_usable_pairs_compared >= $min_to_store) {
	  $mid_matchinfo{$mid} = [$mismatch_count, $n_usable_pairs_compared];
	  if (1) {
	    if (!exists $mid_matchinfosum->{$mid}) {
	      $mid_matchinfosum->{$mid} = [$mismatch_count, $n_usable_pairs_compared, 1];
	    } else {
	      #      print STDERR "ref mid_matchinfosum:  ", ref $mid_matchinfosum, "\n";
	      #            print STDERR "ref mid_matchinfosum->[mid]:  ", ref $mid_matchinfosum->{$mid}, "\n";

	      $mid_matchinfosum->{$mid}->[0] += $mismatch_count;
	      $mid_matchinfosum->{$mid}->[1] += $n_usable_pairs_compared;
	      $mid_matchinfosum->{$mid}->[2]++;
	    }
	  }
	}
	$count++;
      }
    }
  }
  return ($count, \%mid_matchinfo);
}

sub agmr_perl{
  my $gts1 = shift;
  my $gts2 = shift;
  my ($L1, $L2) = (length $gts1, length $gts2);
  die "lengths $L1, $L2 should be equal!" if($L1 != $L2);
  my ($numer, $denom) = (0, 0);
  for (my $i = 0; $i < $L1; $i++) {
    my $c1 = substr($gts1, $i, 1);
    if ($c1 ne 'X') {
      my $c2 = substr($gts2, $i, 1);
      if ($c2 ne 'X') {
	$denom++;
	$numer++ if($c1 ne $c2);
      }
    }
  }
  #   print "PERL: $numer $denom \n";
  return ($denom>0)? $numer/$denom : -1;
}

sub hgmr_perl{
  my $gts1 = shift;
  my $gts2 = shift;
  my ($L1, $L2) = (length $gts1, length $gts2);
  die "lengths $L1, $L2 should be equal!" if($L1 != $L2);
  my ($numer, $denom) = (0, 0);
  for (my $i = 0; $i < $L1; $i++) {
    my $c1 = substr($gts1, $i, 1);
    if ($c1 ne 1  and  $c1 ne 'X') {
      my $c2 = substr($gts2, $i, 1);
      if ($c2 ne 1  and  $c2 ne 'X') {
	#	if($c1 ne 1  and $c2 ne 1){
	$denom++;
	$numer++ if($c1 ne $c2);
	#}
      }
    }
  }
  #  print "PERL: hgmr  $numer $denom \n";
  return ($denom>0)? $numer/$denom : -1;
}

# ############################################
# ########### inline C stuff #################
# ######### has been moved to ################
# ######### ../lib/inlinec.c  ################
