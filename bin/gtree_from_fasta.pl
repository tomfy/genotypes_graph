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
#use Inline 'C';

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

  ### #######  read in genotypes  ###############
  my ($t0, $t1, $t2, $t3, $t4, $t5);
  $t0 = gettimeofday();

  my $input_filename_stem = $input_filename;
  $input_filename_stem =~ s/\.\S+$//; # remove the part after the (last) '.'

  my ($sequence_length, $gobjects) = read_genotypes_from_fasta($input_filename, $the_rng, $p_missing);
  my ($sequence_length2, $gobjects2) = read_genotypes_from_fasta($fasta2, $the_rng, $p_missing);

  print "# input fasta 1: $input_filename\n";
  print "# input fasta 2: $fasta2\n";
  $t1 = gettimeofday();
  #  print STDERR "time to read input, construct genotype objects: ", $t1-$t0, " sec.\n";

  ### do exhaustive comparison: ###
  my $ex_str = '';
  if ($algorithm ne 'quick') {
    my $N = scalar @$gobjects;
    while (my($i, $g1) = each @$gobjects) {
      for my $j (0 .. (@$gobjects2-1)) {
	my $g2 = $gobjects2->[$j];
	my $n_mismatches = $g1->count_mismatches_up_to_limit($g2, $max_mismatches); # inline C
	if ($n_mismatches >= 0) {
	  $exhaustive_eqpairs{$g2->id()} .= $g1->id() . ',';
	} else {		# nuthin
	}
      }
    }
    $t2 = gettimeofday();
    #    print STDERR "time to do N*M genotype-genotype comparisons: ", $t2-$t1, " sec.\n";

    my @sorted_ids = sort { $a <=> $b } keys %exhaustive_eqpairs;
    for my $anid (@sorted_ids) {
      my $s = $exhaustive_eqpairs{$anid};
      $s =~ s/,\s*$//;
      my @eeks = split(',', $s);
      @eeks = sort { $a <=> $b } @eeks;
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
    $gtree = GenotypeTree->new( { depth => $sequence_length } );
    for my $gobj (@$gobjects) {
      $gtree->add_genotype($gobj);
    }
    $t3 = gettimeofday();
    $t2_3 = $t3-$t2;
    #### search in tree for sequences from $fasta2 file ####
    for my $gobj2 (@$gobjects2) {
      $quick_str .= $gobj2->id() . "  " . $gtree->search($gobj2, $max_mismatches) . "\n";
    }
    #  print "exhaust str: [\n$ex_str]\n";
    #  print "quick str:   [\n$quick_str]\n";
    $t4 = gettimeofday();
    $t3_4 = $t4 - $t3;
    $t2_4 = $t4 - $t2;
    #     print "time to search genotype tree: ", $t4-$t3, " sec.\n";

    #    print "root ids:  ", $gtree->root()->id_as_string(), "\n";
    print $gtree->as_newick(), "\n\n" if($newick_out);
  

    #   print STDERR "time to do search tree for N genotypes: ", $t4-$t3, " sec.\n";

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
	my $edge_id = ($id1 < $id2)? $id1 . " " . $id2 : $id2 . " " . $id1;
	$exedges{$edge_id} = 1;
      }
    }
    for my $line (@qs) {
      my ($id1, $id2str) = split(' ', $line);
      my @id2s = split(',', $id2str);
      for my $id2 (@id2s) {
	my $edge_id = ($id1 < $id2)? $id1 . " " . $id2 : $id2 . " " . $id1;
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
  my $n_search = @$gobjects2;
  print "# search_recursive calls: ", $gtree->count_search_recursive_calls(), 
    "  $n_search  ", $gtree->count_search_recursive_calls()/$n_search, "\n" if(defined $gtree);
}                               # end main


####  ############### subroutines ################## ####

sub lose_data{ # with prob. $p_missing replace snp genotypes with MISSING_DATA
  my $sequence = shift;
  my $rng = shift;
  my $p_missing = shift;
  my $L = length $sequence;
  for my $i (0..$L-1) {
    if (gsl_rng_uniform($rng->raw()) < $p_missing) {
      substr($sequence, $i, 1, MISSING_DATA);
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

  #  print STDERR "input file name: [$input_filename]\n";
  my @genotype_objects = ();
  my $sequence_length = undef;
  my @fasta_lines = split ("\n", TomfyMisc::fasta2seqon1line($input_string));
  while (@fasta_lines) {
    my $line = shift @fasta_lines;
    next if($line =~ /^\s*#/);	# skip comment lines
    if ($line =~ /^>(\S+)\s+(\S+)\s+(\S+)/) {
      my ($id, $generation, $pedigree) = ($1, $2, $3);
      my $sequence = shift @fasta_lines;
      #   print "the sequence $sequence \n";
      $sequence =~ s/\s+//g;
      $sequence_length = length $sequence if(!defined $sequence_length);
      die "Sequence lengths must all be the same.\n" if(length $sequence != $sequence_length);
      $sequence = lose_data($sequence, $the_rng, $p_missing);
      push @genotype_objects, Genotype->new(
					    { id => $id,
					      generation => $generation,
					      pedigree => $pedigree,
					      sequence => $sequence }
					   );
    }
  }
  return ($sequence_length, \@genotype_objects);
}


# __DATA__

# ############################################
# ########### inline C stuff #################
# __C__

#   /* number of mismatches between to strings up to some max  */
# /* returns 0,1,2,...,max_mismatches or -1 if there are more mismatches */
# double count_mismatches_up_to_limit_C(int max_mismatches, char* str1, char* str2) {
#   int i = 0;
#   char c1;
#   char c2;
#   int mismatch_count = 0;
#   while(c1 = str1[i]) {
#      c2 = str2[i];
#      if (c1 != c2) {
#        mismatch_count++;
#      }
#      if(mismatch_count > max_mismatches){
#        return -1; // indicates number of mismatches > limit
# }
#   i++;
# }
#   return mismatch_count;
# }

