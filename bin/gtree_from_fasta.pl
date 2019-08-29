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
  my $newick_out = 1;
  my $rng_type = $gsl_rng_mt19937;
  my $rng_seed = undef;
  my $p_missing = -1;		# prob. of missing data at each snp.
  my $algorithm = 'both';	# 'quick' or 'slow' or 'both'

  GetOptions(
	     'input_filename=s' => \$input_filename,
	     'show_all!' => \$show_all,
	     'compact_tree!' => \$compact,
	     'newick_out!' => \$newick_out,
	     'rng_type=s' => \$rng_type,
	     'seed=i' => \$rng_seed,
	     'p_missing=f' => \$p_missing,
	     'algorithm=s' => \$algorithm,
	    );

  my %exhaustive_eqpairs = ();	# 

  # ############## set up the random number generator ##########################################
  if ($rng_type eq 'sys') {
    $rng_type = $gsl_rng_default;
  } elsif ($rng_type eq 'mt') {
    $rng_type = $gsl_rng_mt19937;
  } elsif ($rng_type eq 'lux') {
    $rng_type = $gsl_rng_ranlxd2;
  }
  print STDERR "RNG type: $rng_type \n";
  my $the_rng = (defined $rng_seed)? Math::GSL::RNG->new($rng_type, $rng_seed) : Math::GSL::RNG->new($rng_type) ; # RNG
  # ############################################################################################
my ($t0, $t1, $t2, $t3, $t4);
  $t0 = gettimeofday();

  my $input_filename_stem = $input_filename;
  $input_filename_stem =~ s/\.\S+$//; # remove the part after the (last) '.'
  my $input_string = '';
  {
    if (defined $input_filename) {
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
      next if($line =~ /^\s*#/); # skip comment lines
      if ($line =~ /^>(\S+)\s+(\S+)\s+(\S+)/) {
	my ($id, $generation, $pedigree) = ($1, $2, $3);
	my $sequence = shift @fasta_lines;
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
  } else {
    die "Input file must be fasta format. \n";
  }
  $t1 = gettimeofday();
  print STDERR "time to read input, construct genotype objects: ", $t1-$t0, " sec.\n";

  ### do exhaustive comparison: ###
  my $ex_str = '';
  if ($algorithm ne 'quick') {
    my $N = scalar @genotype_objects;
    while (my($i, $g1) = each @genotype_objects) {
      for my $j ($i .. $N-1) {
	my $g2 = $genotype_objects[$j];
	if (compare_two_genotype_objects($g1, $g2) eq 'equal') {
	  $exhaustive_eqpairs{$g1->id()} .= $g2->id() . ',';
	  $exhaustive_eqpairs{$g2->id()} .= $g1->id() . ',' if($j != $i);
	} else {		# nuthin
	}
      }
    }
    $t2 = gettimeofday();
    print STDERR "time to do N choose 2 genotype-genotype comparisons: ", $t2-$t1, " sec.\n";

    my @sorted_ids = sort { $a <=> $b } keys %exhaustive_eqpairs;
    for my $anid (@sorted_ids) {
      my $s = $exhaustive_eqpairs{$anid};
      $s =~ s/,\s*$//;
      my @eeks = split(',', $s);
      @eeks = sort { $a <=> $b } @eeks;
      $ex_str .=  "$anid  " . join(",", @eeks) . "\n";
    }
    print "[$ex_str]\n";
  }
  $t2 = gettimeofday();
  #### end of exhaustive ####
  
  my $quick_str = '';
  if ($algorithm ne 'slow') {
    my $gtree = GenotypeTree->new( { depth => $sequence_length } );
    for my $gobj (@genotype_objects) {
      #     $gobj->lose_data($the_rng, 0.5);
      if (!$compact) {
	$gtree->add_genotype($gobj);
      } else {
	$gtree->add_genotype_compact($gobj);
      }
    }
    $t3 = gettimeofday();
    print STDERR "time to do construct genotype tree: ", $t3-$t2, " sec.\n";
  
    print "root ids:  ", $gtree->root()->id_as_string(), "\n";
    print $gtree->as_newick(), "\n\n" if($newick_out);
    $t3 = gettimeofday();
 
    print "Now search for genotypes in tree: \n";
 
    for my $gobj (@genotype_objects) {
      $quick_str .= $gtree->search($gobj);
    }
    $t4 = gettimeofday();
    print STDERR "time to do search tree for N genotypes: ", $t4-$t3, " sec.\n";
    print "[$quick_str]\n";
  }
  print STDERR "algorithm: [$algorithm] \n";
  if ($algorithm eq 'both') {
    if ($quick_str eq $ex_str) {
      print STDERR "Strings are equal.\n";
    } else {
      print STDERR "Strings are not equal.\n";

      my @exs = split("\n", $ex_str);
      my @qs = split("\n", $quick_str);
      while(my($i, $q) = each @qs){
	my $x = $exs[$i];
	print STDERR "$q  $x \n";
	exit if($q ne $x);
      }
    }
  }
  print STDERR "    input/gobjs   Nchoose2   treeconstruct  treesearch      total \n";
  printf(STDERR "%12.3f %12.3f %12.3f %12.3f %12.3f \n",$t1-$t0, $t2-$t1, $t3-$t2, $t4-$t3, $t4-$t0);
}                               # end main



######################

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

sub compare_two_genotype_objects{
  my $gobj1 = shift;
  my $gobj2 = shift;
  my $gs1 = $gobj1->sequence();
  my $gs2 = $gobj2->sequence();
  my $L1 = length $gs1;
  die if(length $gs2 != $L1);
  my $n_equal_characters = 0;
  for my $i (0..$L1-1) {
    my ($c1, $c2) = ( substr($gs1, $i, 1), substr($gs2, $i, 1) );
    last if(($c1 ne MISSING_DATA)  and  ($c2 ne MISSING_DATA)  and  ($c1 ne $c2));
    $n_equal_characters++;
  }
  return ($n_equal_characters == $L1)? 'equal' : 'unequal';
}
