#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw (min max sum);
use Getopt::Long;
use Time::HiRes qw( gettimeofday );
use lib '/home/tomfy/Orthologger/lib/';
use TomfyMisc qw ' fasta2seqon1line ';
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
  my $show_all = 0;	# false -> output info on only the leaf nodes.
  my $compact = 1;
  my $newick_out = 0;
  my $rng_type = $gsl_rng_mt19937;
  my $rng_seed = undef;
  my $p_missing = 0;		# prob. of missing data at each snp.
  my $algorithm = 'both';	# 'quick' or 'slow' or 'both'
  my $fasta2 = undef;		# query genotype sets fasta file.
  my $max_mismatches = 0;
  my $use_inline_C = 1;
  my $error_prob = 0;	    # add errors to data with this probability
  my $chunk_size = 100000000;
  my $n_chunks = 1;
  my $sort_markers = 0;
  my $min_usable_pairs_to_store = 2;

  GetOptions(
	     'input_filename|fasta1|f1=s' => \$input_filename,
	     'fasta2|f2=s' => \$fasta2, # queries
	     'show_all!' => \$show_all,
	     'compact_tree!' => \$compact,
	     'newick_out!' => \$newick_out,
	     'rng_type=s' => \$rng_type,
	     'seed=i' => \$rng_seed,
	     'p_missing=f' => \$p_missing,
	     'algorithm=s' => \$algorithm,
	     'max_mismatch=i' => \$max_mismatches,
	     'C|useC!' =>  \$use_inline_C,
	     'error_prob=f' => \$error_prob,
	     'chunk_size=i' => \$chunk_size,
	     'n_chunks=i' => \$n_chunks,
	     'sort_markers!' => \$sort_markers,
	     'min_to_store=i' => \$min_usable_pairs_to_store,
	    );

  # $max_mismatches--; # this is so the command line parameter is the number of mismatches 

  #  print STDERR "input file name: $input_filename\n";

  #my %exhaustive_eqpairs = ();	# 
  my $ok = '?';
  # ############## set up the random number generator ##########################################
  if ($rng_type eq 'sys') {
    $rng_type = $gsl_rng_default;
  } elsif ($rng_type eq 'mt') {
    $rng_type = $gsl_rng_mt19937;
  } elsif ($rng_type eq 'lux') {
    $rng_type = $gsl_rng_ranlxd2;
  }
  #  print STDERR "RNG type: $rng_type \n";
  my $the_rng = (defined $rng_seed)? Math::GSL::RNG->new($rng_type, $rng_seed) : Math::GSL::RNG->new($rng_type) ; # RNG

  print "# genotype input files: $input_filename  $fasta2 \n";
  print "# p_missing: $p_missing   error_prob: $error_prob \n";
  print "# max_mismatch: $max_mismatches \n";
  print "# chunk_size: ", (defined $chunk_size)? $chunk_size : 'undef',   "  n_chunks $n_chunks\n";

  ### #######  read in genotypes  ###############
  my ($t0, $t1, $t2, $t3, $t4, $t5);
  my $t_output = 0;
  $t0 = gettimeofday();

  my @array_0123etc = ();

  my $input_filename_stem = $input_filename;
  $input_filename_stem =~ s/\.\S+$//; # remove the part after the (last) '.'

  my ($sequence_length, $id_gtsobj1, $mindex_okcount1) = read_genotypes_from_fasta($input_filename, $the_rng, $p_missing, $error_prob);
  for (0..$sequence_length-1) {
    $array_0123etc[$_] = $_;
  }
  my ($sequence_length2, $id_gtsobj2, $mindex_okcount2) = read_genotypes_from_fasta($fasta2, $the_rng, $p_missing, $error_prob);
  my $n_accessions1 = scalar keys %$id_gtsobj1;
  my $n_accessions2 = scalar keys %$id_gtsobj2;
  my $mindex_okcount = sum_hashes($mindex_okcount1, $mindex_okcount2);
  my @sorted_match_indices = ($sort_markers)?
    sort {$mindex_okcount->{$b} <=> $mindex_okcount->{$a} }
    keys %$mindex_okcount	# sort by okcount (high to low)
    :
    @array_0123etc;
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

  

  $t1 = gettimeofday();
  #  print STDERR "time to read input, construct genotype objects: ", $t1-$t0, " sec.\n";

  ### do exhaustive comparison: ###
  my %exhaust_chunk_results = ();
  if ($algorithm ne 'quick') {
    my $N2 = scalar keys %$id_gtsobj2;
    while (my($i_chunk, $chunk_indices) = each @chunk_index_arrayrefs) {
      my %qid_matchidinfo = ();
      #   $for my $j (0 .. ($N2-1)) {
      #	$my $g2 = $go->[$j];
      while (my ($qid, $q_gtsobj) = each %$id_gtsobj2) {
	my $mid_matchinfo = {};
	$mid_matchinfo = exhaustive_search_1query($id_gtsobj1, $q_gtsobj, $chunk_indices, $max_mismatches);
	# while(my($mid, $minfo) = each %$mid_matchinfo){
	#   print STDERR "X $i_chunk  $qid $mid [$minfo]\n";
	# }
	$qid_matchidinfo{$q_gtsobj->id()} = $mid_matchinfo;
      }
      $exhaust_chunk_results{$i_chunk} = \%qid_matchidinfo;
    }
  }
  $t2 = gettimeofday();
  #### end of exhaustive search ####

 
  #### construct tree ####
  my $n_nogood_total = 0; # number of queries for which exhaustive and tree searches disagree.
  my ($t_treeconstruct, $t_treesearch, $t_check) = (0, 0, 0);
  if ($algorithm ne 'slow') {	# search using trees
    my %qid_mid_sums;
    while (my($i_chunk, $chunk_indices) = each @chunk_index_arrayrefs) {
      print STDERR "# doing chunk $i_chunk\n";
      my $ta = gettimeofday();
      my $gtree = GenotypeTree->new( { depth => scalar @$chunk_indices } );
      # for my $gobj (@$gobjects1) {
      while (my ($qid, $gobj1) = each %$id_gtsobj1) {
	$gtree->add_genotype($gobj1, $chunk_indices);
      }
      $t_treeconstruct += gettimeofday - $ta;
      my %ngoodbad_count = ();
      #### search in tree for sequences from $fasta2 file ####

      # for my $gobj2 (@$gobjects2) {
      while (my ($id2, $gobj2) = each %$id_gtsobj2) {
	# my $id2 = $gobj2->id();
	$qid_mid_sums{$id2} = {} if(! exists $qid_mid_sums{$id2});
	my $tb = gettimeofday();
	my $quick_mid_matchinfo = $gtree->search($gobj2, $chunk_indices, $max_mismatches, $min_usable_pairs_to_store);
	# matchinfo is  $max_bad_count  $bad_count  ($n_gts_above+$n_ok_characters)   $n_good_to_exceed_max_bad_count";
	while (my($mid, $minfo) = each %$quick_mid_matchinfo) {
	  my ($n_mms, $n_ok_before, $N) = split(" ", $minfo);
	  $qid_mid_sums{$id2}->{$mid} = [0, 0, 0] if(! exists $qid_mid_sums{$id2}->{$mid});
	  $qid_mid_sums{$id2}->{$mid}->[0] += $n_mms; # number of mismatches found, (sum over chunks)
	  $qid_mid_sums{$id2}->{$mid}->[1] += $N; # number of usable gt pairs (i.e. with neither being missing data) compared (sum over chunks)
	  $qid_mid_sums{$id2}->{$mid}->[2]++; # number of chunks with at least $min_to_store usable pairs compared.

	  #  print STDERR "T $i_chunk  $id2 $mid [$minfo]\n";
	}

	my $tc = gettimeofday();
	$t_treesearch += $tc - $tb;

	####### Check whether match info can predict best agmrs ##############
	#	while(my($qid, $mid_sums) = each %qid_mid_sums){
	#	  print STDERR
	#	}

	######## compare tree results with exhaustive results: #####################
	if ($algorithm eq 'both') {
	#  while(my($x, $y) = each $exhaust_chunk_results->{0}
	  my ($ok_count, $nogood_count) = (0, 0);
	  while (my($mid, $v) = each %$quick_mid_matchinfo) {
	    my $xv = $exhaust_chunk_results{$i_chunk}->{$gobj2->id()}->{$mid} // 'undef';
	    if ($xv eq $v) {
	      $ok_count++;
	    } else {
	      print "match id: $mid  exhaustive: $xv  tree: $v \n";
	      $nogood_count++;
	    }
	  }

	  $ngoodbad_count{"$ok_count    $nogood_count"}++;
	  $n_nogood_total += $nogood_count;
	}
	$t_check += gettimeofday - $tc;
      }

      print "# chunk  ngood nbad  count \n" if($i_chunk == 0);
      while (my($ngb, $count) = each %ngoodbad_count) {
	print "#   $i_chunk     $ngb    $count \n";
      }

      print $gtree->as_newick(), "\n\n" if($newick_out);

    }				# end of loop over chunks

    # report results
    my %mid_agmr = ();
    my %mid_agmr_fom = ();
    my %mid_hgmr = ();
    my %mid_hgmr_fom = ();
    my $td = gettimeofday();
    while (my($qid, $mid_sums) = each %qid_mid_sums) {
      #print STDERR "$qid \n";
      my $q_sequence = $id_gtsobj2->{$qid}->sequence();
      my $neitherXnor1_indices = neitherXnor1_indices($q_sequence);
      while (my($mid, $sums) = each %$mid_sums) {
	my $m_sequence = $id_gtsobj1->{$mid}->sequence();
#	print STDERR "qid: $qid   ", substr($q_sequence, 0, 30), "\n";
#	print STDERR "mid:   $mid   ", substr($m_sequence, 0, 30), "\n";
	#	my $agmr = agmr($q_sequence, $m_sequence);
#	my $perl_hgmr = hgmr($q_sequence, $m_sequence);
	
	my ($agmr_n, $agmr_d, $hgmr_n, $hgmr_d) = (0, 0, 0, 0);
#	agmr_C($q_sequence, $m_sequence, $agmr_n, $agmr_d);
#	hgmr_C($q_sequence, $m_sequence, $hgmr_n, $hgmr_d);
#        hgmr_Cx($q_sequence, $neitherXnor1_indices, $m_sequence, $hgmr_n, $hgmr_d);
	agmr_hgmr_C($q_sequence, $m_sequence, $agmr_n, $agmr_d, $hgmr_n, $hgmr_d);
	my $agmr = ($agmr_d > 0)? $agmr_n/$agmr_d : -1;
	my $hgmr = ($hgmr_d > 0)? $hgmr_n/$hgmr_d : -1;
	
#	my ($hgmr1_n, $hgmr1_d) = (0, 0);
#	hgmr_C($q_sequence, $m_sequence, $hgmr1_n, $hgmr1_d);
#	my $hgmr1 = ($hgmr1_d > 0)? $hgmr1_n/$hgmr1_d : -1;
	# query_id  match_id   n_mismatches  n_usable_pairs_compared  n_chunks_info_stored
	my ($total_stored_mms, $n_usable_compared, $n_stored_chunks) = @$sums;
	my $n_unstored_chunks = $n_chunks_used - $n_stored_chunks;
	my $total_mms = $total_stored_mms + $n_unstored_chunks*$max_mismatches;
	my $denom = $n_usable_compared + $n_unstored_chunks*($min_usable_pairs_to_store-1);
	my $fom = ($denom > 0)? $total_mms/$denom : 2; #figure of merit: small fom should imply small agmr.
	$mid_agmr{$mid} = $agmr;
	$mid_agmr_fom{$mid} = $fom;
	print STDERR "   $fom $agmr $hgmr \n"
	# printf( STDERR "  %s %2i %3i %7.3f %1i  %2i %3i  %8.4f  %8.5f\n",
	# 	$mid, $total_stored_mms, $n_usable_compared, (($n_usable_compared)? $total_stored_mms/$n_usable_compared : -1),
	# 	$n_stored_chunks, $total_mms, $denom, $fom, $agmr);
      }
      my @sorted_mids_agmr = sort {$mid_agmr{$a} <=> $mid_agmr{$b}} keys %mid_agmr;
      my @sorted_mids_fom = sort {$mid_agmr_fom{$a} <=> $mid_agmr_fom{$b}} keys %mid_agmr_fom;
      #   print STDERR "AAAA: ", join(",", @sorted_mids_agmr), "\n";
      #   print STDERR "FFFF: ", join(",", @sorted_mids_fom), "\n";
      
      # for (@sorted_mids_agmr[0..8]) {
      # 	printf( STDERR "%24s %8.6f  ", $_, $mid_agmr{$_});
      # }
      # print STDERR "\n";
      # for (@sorted_mids_fom[0..8]) {
      # 	printf( STDERR "%24s %8.6f  ", $_, $mid_agmr_fom{$_});
      # }
      # print STDERR "\n";
      
    }
    $t_output += gettimeofday() - $td;
  }				### end of tree search ###

  print STDERR "#     input,gobjs      N*M       treeconstr    search        check       output    OK?\n";
  printf(STDERR "# %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f   %3s\n",
	 $t1-$t0, $t2-$t1, $t_treeconstruct, $t_treesearch, $t_check, $t_output, ($algorithm eq 'both')? (($n_nogood_total == 0)? 'Y' : 'N') : '-' );
}				# end main


####  ############### subroutines ################## ####
sub notX_indices{
  my $sequence = shift;
  my @notxindices = ();
    for(my $i=0; $i< length $sequence; $i++){
      push @notxindices, $i if(substr($sequence, $i, 1) ne MISSING_DATA);
    }
  return \@notxindices;
}

sub neitherXnor1_indices{
  my $sequence = shift;
  my @notx1indices = ();
  for(my $i=0; $i< length $sequence; $i++){
    my $gt = substr($sequence, $i, 1);
    next if($gt eq '1');
      push @notx1indices, $i if($gt ne MISSING_DATA);
    }
  return \@notx1indices;
}

sub lose_data{ # with prob. $p_missing replace snp genotypes with MISSING_DATA
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

sub add_errors_to_data{ # with prob. $error_prob replace snp genotypes with other (wrong) gts.
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
  my @fasta_lines = split ("\n", TomfyMisc::fasta2seqon1line($input_string));
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
  my $index_arrayref = shift;
  my $max_mismatches = shift;
  my $q_sequence = $q_gts_obj->get_chunk($index_arrayref); # $chunk_index_arrayrefs[0]);
  my %mid_matchinfo = ();
  while (my($i, $db_gts_obj) = each %$db_gts_objs) {
    # my $sequence1 = $db_gt_obj->get_chunk($index_arrayref);
    my ($n_pairs_compared, $n_usable_pairs_compared, $mismatch_count) = (0, 0, 0);
    #GenotypeTreeNode::
    count_oks_and_mismatches_up_to_limit_C(
							     $db_gts_obj->get_chunk($index_arrayref),
							     $q_sequence,
							     # $sequence1,
							     # $sequence2,
							     $max_mismatches, 0,
							     $n_pairs_compared, $n_usable_pairs_compared, $mismatch_count);

    $mid_matchinfo{$db_gts_obj->id()} = "$mismatch_count  $n_pairs_compared  $n_usable_pairs_compared";
  }
  return \%mid_matchinfo;
}

sub agmr{
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

sub hgmr{
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

__DATA__

# ############################################
# ########### inline C stuff #################
__C__

  #define MISSING_DATA 'X'
  void agmr_C(char* str1, char* str2, SV* numer, SV* denom){
    char c1;
    char c2;
    int i=0;
    int n=0;
    int d=0;
    while(c1 = str1[i]){
      if(c1 != MISSING_DATA){
	c2 = str2[i];
	if(c2 != MISSING_DATA){
	  d++;
	  if(c1 != c2){
	    n++;
	  }
	}
      }
      i++;
    }
 //   printf("C:    %i %i \n", n, d);
    sv_setiv(numer, n);
    sv_setiv(denom, d);
  }

void hgmr_Cx(char* str1, AV* qindices, char* str2, SV* numer, SV* denom){ // slower!
  /* qindices is array of the indices for which str1 (query) is homozygous */

    char c1;
    char c2;
//    int i=0;
    int n=0;
  int d=0;
  int L = av_len(qindices);
  for(int i=0; i<L; i++){
  //  v = ;
      int idx = SvNV( * av_fetch( qindices, i, 0) );
      c1 = str1[idx];
      c2 = str2[idx];
      if((c2 != '1') && (c2 != MISSING_DATA)){
	d++;
	if(c1 != c2){
	  n++;
	}
      }
    }
 //   printf("C:    %i %i \n", n, d);
    sv_setiv(numer, n);
    sv_setiv(denom, d);
}

void hgmr_C(char* str1, char* str2,
	    SV* numer, SV* denom){
    char c1;
    char c2;
    int i=0;
    int n=0;
    int d=0;
    while(c1 = str1[i]){
      if((c1 != '1') && (c1 != MISSING_DATA)){
	c2 = str2[i];
	if((c2 != '1') && (c2 != MISSING_DATA)){
	  d++;
	  if(c1 != c2){
	    n++;
	  }
	}
      }
      i++;
    }
    sv_setiv(numer, n);
    sv_setiv(denom, d);
  }

void agmr_hgmr_C(char* str1, char* str2, SV* a_numer, SV* a_denom, SV* h_numer, SV* h_denom){
    char c1;
    char c2;
    int i=0;
    int n_a=0;
    int r_a=0;
    int n_h=0;
    int r_h=0;
    while(c1 = str1[i]){
      if(c1 != MISSING_DATA){
	c2 = str2[i];
	if(c2 != MISSING_DATA){

	  if(c1 != c2){
	    n_a++;
	     if((c1 != '1')  &&  (c2 != '1')){
	      n_h++;
	    }
	  }else{
	    r_a++;
	     if((c1 != '1')  &&  (c2 != '1')){
	      r_h++;
	    }
	  }

	}
      }
      i++;
    }
 //   printf("C:    %i %i \n", n, d);
    sv_setiv(a_numer, n_a);
    sv_setiv(a_denom, r_a+n_a);
    sv_setiv(h_numer, n_h);
    sv_setiv(h_denom, r_h+n_h);
  }

