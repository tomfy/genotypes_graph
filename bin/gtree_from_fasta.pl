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
use Inline 'C';

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
  my $p_missing = -1;		# prob. of missing data at each snp.
  my $algorithm = 'both';	# 'quick' or 'slow' or 'both'
  my $fasta2 = undef;
  my $max_mismatches = 0;
  my $use_inline_C = 1;
  my $error_prob = 0;	    # add errors to data with this probability
  my $chunk_size = undef;
  my $n_chunks = 1;
  my $sort_markers = 0;

  GetOptions(
	     'input_filename|fasta1|f1=s' => \$input_filename,
	     'fasta2|f2=s' => \$fasta2,
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
	    );

  #  print STDERR "input file name: $input_filename\n";

  my %exhaustive_eqpairs = ();	# 
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
  print "# chunk_size: $chunk_size  n_chunks $n_chunks\n";

  ### #######  read in genotypes  ###############
  my ($t0, $t1, $t2, $t3, $t4, $t5);
  $t0 = gettimeofday();

  my @array_0123etc = ();

  my $input_filename_stem = $input_filename;
  $input_filename_stem =~ s/\.\S+$//; # remove the part after the (last) '.'

  my ($sequence_length, $gobjects1, $mindex_okcount1) = read_genotypes_from_fasta($input_filename, $the_rng, $p_missing, $error_prob);
  for(0..$sequence_length-1){ $array_0123etc[$_] = $_; }
  my ($sequence_length2, $gobjects2, $mindex_okcount2) = read_genotypes_from_fasta($fasta2, $the_rng, $p_missing, $error_prob);
  my $n_accessions1 = scalar @$gobjects1;
  my $n_accessions2 = scalar @$gobjects2;
  my $mindex_okcount = sum_hashes($mindex_okcount1, $mindex_okcount2);
  my @sorted_match_indices = ($sort_markers)?
    sort {$mindex_okcount->{$b} <=> $mindex_okcount->{$a} }
    keys %$mindex_okcount # sort by okcount (high to low)
    :
    @array_0123etc;
  my @chunk_index_arrayrefs = ();
  
  my $n_chunks_used = $n_chunks;
  my $i_start = 0;
  for my $i_chunk (1..$n_chunks) {
    my $i_end = $i_start + $chunk_size - 1;
    if ($i_end >= scalar @sorted_match_indices) {
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
  my $ex_str = '';
  if ($algorithm ne 'quick') {
    my $N = scalar @$gobjects1;
    while (my($i, $g1) = each @$gobjects1) {
      for my $j (0 .. (@$gobjects2-1)) {
	my $g2 = $gobjects2->[$j];
	my $n_mismatches;
	# if(0){
	#   my ($n_ok, $n_mm); # = $g1->count_mismatches_up_to_limit($g2, $max_mismatches); # inline C
	#   #count_mismatches_up_to_limit($g2, $max_mismatches); # inline C
	#   count_oks_and_mismatches_up_to_limit_C($g1->sequence(), $g2->sequence(), $max_mismatches, $n_ok, $n_mm);
	#   #  print "$i $j $n_ok $n_mm \n";
	#   $n_mismatches = $n_mm;
	# }else{
	$n_mismatches = count_mismatches_up_to_limit_C($g1->sequence(), $g2->sequence(), $max_mismatches);
	#	}
	if ($n_mismatches >= 0) {
	  $exhaustive_eqpairs{$g2->id()} .= $g1->id() . ',';
	} else {		# nuthin
	}
      }
    }
    $t2 = gettimeofday();
    #    print STDERR "time to do N*M genotype-genotype comparisons: ", $t2-$t1, " sec.\n";

    my @sorted_ids = sort  keys %exhaustive_eqpairs;
    for my $anid (@sorted_ids) {
      my $s = $exhaustive_eqpairs{$anid};
      $s =~ s/,\s*$//;
      my @eeks = split(',', $s);
      @eeks = sort  @eeks;
      $ex_str .=  "$anid  " . join(",", @eeks) . "\n" if(scalar @eeks > 0);
    }
  }
  $t2 = gettimeofday();
  #### end of exhaustive ####

  my $gtree = undef;
  my ($t2_3, $t3_4, $t2_4) = ('---', '---', '---');
  #### construct tree ####
  my $quick_str = '';
  if ($algorithm ne 'slow') {
    #   print $sequence_length // 'XXX', " \n";
    while (my($i_chunk, $chunk_indices) = each @chunk_index_arrayrefs) {
      print STDERR "\n", "# doing chunk $i_chunk\n";
      my $gtree = GenotypeTree->new( { depth => $sequence_length } );
      for my $gobj (@$gobjects1) {
	$gtree->add_genotype($gobj, $chunk_indices);
      }
      $t3 = gettimeofday();
      $t2_3 = $t3-$t2;
      #### search in tree for sequences from $fasta2 file ####
      for my $gobj2 (@$gobjects2) {
	$quick_str .= $gobj2->id() . "  " . $gtree->search($gobj2, $max_mismatches, $chunk_indices) . "\n";
      }
      #print "exhaust str: [\n$ex_str]\n";
      # print "quick str:   [\n$quick_str]\n";
      $t4 = gettimeofday();
      $t3_4 = $t4 - $t3;
      $t2_4 = $t4 - $t2;
      #     print "time to search genotype tree: ", $t4-$t3, " sec.\n";

      #    print "root ids:  ", $gtree->root()->id_as_string(), "\n";
      print $gtree->as_newick(), "\n\n" if($newick_out);
  

      #   print STDERR "time to do search tree for N genotypes: ", $t4-$t3, " sec.\n";

    } # end of loop over chunks
  } 
  if ($algorithm eq 'both') {
    # check that quick gives same matches as exhaustive
    my %exedges = ();
    my %qedges = ();
    my @exs = split("\n", $ex_str);
    my @qs = split("\n", $quick_str);

    for my $line (@exs) {
      my ($id1, $id2str) = split(' ', $line);
      my @id2s = split(',', $id2str);
      for my $id2 (@id2s) {
	my $edge_id = ($id1 lt $id2)? $id1 . " " . $id2 : $id2 . " " . $id1;
	$exedges{$edge_id} = 1;
      }
    }
    for my $line (@qs) {
      my ($id1, $id2str) = split(' ', $line);
      my @id2s = split(',', $id2str);
      for my $id2 (@id2s) {
	my $edge_id = ($id1 lt $id2)? $id1 . " " . $id2 : $id2 . " " . $id1;
	$qedges{$edge_id} = 1;
      }
    }

    my $bad_edge = undef;
    for (keys %exedges) {
      if (!exists $qedges{$_}) {
	#	$ok = 'not equal';
	$bad_edge = $_;
	last;
      }
    }
    if (!defined $bad_edge) {
      for (keys %qedges) {
	if (!exists $exedges{$_}) {
	  #	$ok = 'not equal';
	  $bad_edge = $_;
	  last;
	}
      }
    }
    $ok = (defined $bad_edge)? 'N' : 'Y';
  }
  $t5 = gettimeofday();

  if ($algorithm eq 'both') {
    print STDERR "#     input,gobjs      N*M       treeconstr    search     constr&search   check        total    agree? \n";
    printf(STDERR "# %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f      %s\n",
	   $t1-$t0, $t2-$t1, $t2_3, $t3_4, $t2_4, $t5-$t4, $t5-$t0, $ok);
  } elsif ($algorithm eq 'slow') {
    print STDERR "#     input,gobjs      N*M         total  \n";
    printf(STDERR "# %12.3f %12.3f %12.3f \n", $t1-$t0, $t2-$t1, $t5-$t0);
  } else {			#  quick only 
    print STDERR "#     input,gobjs   treeconstr    search     constr&search   total\n";
    printf(STDERR "# %12.3f %12.3f %12.3f %12.3f %12.3f \n",
	   $t1-$t0, $t2_3, $t3_4, $t2_4, $t5-$t0);
  }
  #  my $n_search = @$gobjects2;
  #  print "# search_recursive calls: ", $gtree->count_search_recursive_calls(), 
  #    "  $n_search  ", $gtree->count_search_recursive_calls()/$n_search, "\n" if(defined $gtree);
} # end main


####  ############### subroutines ################## ####

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
  my @genotype_objects = ();
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
    push @genotype_objects, Genotype->new(
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
  return ($sequence_length, \@genotype_objects, \%mindex_okcount);
}

sub sum_hashes{
  my $hr1 = shift;
  my $hr2 = shift;

  while (my ($k, $v) = each %$hr2) {
    $hr1->{$k} += $v;
  }
  return $hr1;
}


__DATA__

# ############################################
# ########### inline C stuff #################
__C__

  #define MISSING_DATA 'X'

  /* number of mismatches between two strings up to some max  */
 /* returns 0,1,2,...,max_mismatches or -1 if there are more mismatches */
int count_mismatches_up_to_limit_C(char* str1, char* str2, int max_mismatches) {
  int i = 0;
  char c1;
  char c2;
  int mismatch_count = 0;
  while(c1 = str1[i]) {
    if(c1 != MISSING_DATA){
      c2 = str2[i];
      if(c2 != MISSING_DATA){
        if (c1 != c2) {
          mismatch_count++;
        }
        if(mismatch_count > max_mismatches){
          return -1; // indicates number of mismatches > limit
}
}
}
  i++;
}
  return mismatch_count;
}

  // lines to use in alternative way of get multiple values back from Inline::C function
  /* Inline_Stack_Vars;
Inline_Stack_Reset;
Inline_Stack_Push(sv_2mortal(newSViv(i)));
Inline_Stack_Push(sv_2mortal(newSViv(mismatch_count)));
Inline_Stack_Done; */

