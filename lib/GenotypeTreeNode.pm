package GenotypeTreeNode;
#use Moose;
use Mouse;
use namespace::autoclean;
use Carp::Assert;
use List::Util qw ( min max sum );
use constant MISSING_DATA => 'X';
use Inline 'C';

no warnings 'recursion';

use constant DEBUG => 0;

# a class for nodes of a tree.
has ids => (
            isa => 'ArrayRef',
            is => 'ro',
            default => sub { [] },
           );

has count => (
	      #      traits => ['Counter'],
              isa => 'Int',
              is => 'rw',
              default => 0,
              # handles => {
              #             inc_counter => 'inc',
              #            },
             );

has depth => ( # 0 for root, counts not number of nodes below root, but number of snps, so all leaves have depth == length of genotype sequences.
              isa => 'Int',
              is => 'rw',
              required => 1,
             );

has genotype => ( # the one or more snp genotype associated with this node. String is 'R' for the root, otherwise a sequence of 0, 1, 2, and MISSING_DATA.
                 isa => 'Str',
                 is => 'rw',
                 required => 1,
                );

has parent => (
               isa => 'Maybe[Object]',
               is => 'rw',
               default => undef,
              );

has children => (
                 isa => 'HashRef', # keys are 0, 1, 2, or '-', values are GenotypeTreeNode objects
                 is => 'rw',
                 default => sub { {} },
                );

has tree => (
	     isa => 'Object',
	     is => 'rw',
	     required => 0,
	    );

# sub add_child{                  # 
#   my $self = shift;
#   my $id = shift;
#   my $g = shift;		# string of 0, 1, 2, or '-'
#   my $d = $self->depth();	# depth of this node
#   my $child_node;

#   if (!defined $self->children()->{$g}) {
#     $child_node = GenotypeTreeNode->new( {tree => $self->tree(),
# 					  parent => $self, depth => $d+1, genotype => $g, ids => [$id] } );
#     $self->children()->{$g} = $child_node;
#   } else {
#     $child_node = $self->children()->{$g};
#     push @{$child_node->ids()}, $id;
#   }
#   #  $child_node->inc_counter();
#   $child_node->{count}++;
#   return $child_node;
# }

sub add_genotype_recursive{
  my $self = shift;
  my $id2 = shift;
  my $gt2s = shift;
  my $gt1s = $self->genotype();
  my $d = $self->depth();
  my $parent_depth = $self->parent()->depth();
  my $n_equal_characters = 0;

  for my $i (0 .. (length $gt1s) - 1) {
    last if(substr($gt1s, $i, 1) ne substr($gt2s, $i, 1));
    $n_equal_characters++;
  }

  my ($l1, $l2) = (length $gt1s, length $gt2s);

  if ($n_equal_characters < length $gt1s) { # $gt1s and beginning of $gt2s are not the same. Make new node between self and parent.
    # $self gets two new children, child1 gets $self's children, and child2 is a leaf (no children)
    my $g_common = substr($gt1s, 0, $n_equal_characters, '');
    my $g2_common = substr( $gt2s, 0, $n_equal_characters, '');
    assert ($g2_common eq $g_common) if DEBUG;

    #   my @new_child1_ids = @{ $self->ids() };
    my $new_child1 = GenotypeTreeNode->new( {tree => $self->tree(),
					     parent => $self,
					     depth => $d,
					     genotype => $gt1s,
					     #      ids => \@new_child1_ids,
					     ids => [@{$self->ids()}], 
					     children => $self->children() } );

    my $new_child2 = GenotypeTreeNode->new( {tree => $self->tree(),
					     parent => $self,
					     depth => $parent_depth + $n_equal_characters + length $gt2s,
					     genotype => $gt2s,
					     ids => [$id2] } );

    while (my($g, $ch) = each %{$self->children()}) { # these need to make $new_child1 their parent!
      $ch->parent($new_child1);
    }
    $self->genotype($g_common);
    $self->depth($parent_depth + $n_equal_characters);
    $self->add_id($id2);

    $self->children( { substr($gt1s, 0, 1) => $new_child1, substr($gt2s, 0, 1) => $new_child2 } );

  } else {		# $gt1s and beginning of $gt2 are the same, so
    substr($gt2s, 0, $n_equal_characters, '');
    if ((length $gt2s) > 0) {	# theres a bit of $gt2s left
      my $gt2head = substr($gt2s, 0, 1);
      $self->add_id($id2);
      if (exists $self->children()->{$gt2head}) {
	$self->children()->{$gt2head}->add_genotype_recursive($id2, $gt2s);
      } else {
	my $new_child = 
	  GenotypeTreeNode->new( {tree => $self->tree(),
				  parent => $self,
				  depth => $d + length $gt2s, 
				  genotype => $gt2s,
				  ids => [$id2] } );
	$self->children()->{$gt2head} = $new_child;
      }
    } else {			# this must be a leaf node
      push @{$self->ids()}, $id2;
    }
  }
}

sub search_recursive{
  my $self = shift;
  my $id2 = shift;
  my $gt2s = shift;		# searching for this one.
  my $max_bad_count = shift // die; # the max number of mismatching characters allowed along this branch
  my $matchid_matchinfo = shift; # keys query id, values: string: [id_string n_to_spare] (n_to_spare is $max_bad_count - $bad_count)
  my $n_gts_above = shift;
  my $n_good_gts_above = shift;
   my $gt1s = $self->genotype(); # the genotype along the branch above this node.
  my ($L1, $L2) = (length $gt1s, length $gt2s); # $L1: gts on this branch (above node); $L2 gts left in query.
  $self->{tree}->{count_search_recursive_calls}++;

  assert ($L2 >= $L1) if DEBUG;
  assert ($L2 == $L1  or (keys %{$self->children()} > 0)) if DEBUG;

  print STDERR "\n", "newick: ", $self->newick_recursive(), "\n";

  my $matching_id_string = '';
  my $bad_count = 0; # the number of mismatching characters so far along this branch.
  my $n_ok_characters = 0; # the number of characters before the bad_count gets too big.
  my $n_good_pairs_to_limit = 0; # number of gt pairs (neither being MISSING_DATA) to exceed max_bad
  if (0) {  #  perl
    for my $i (0 .. $L1 - 1) {
      my ($g1, $g2) = (substr($gt1s, $i, 1), substr($gt2s, $i, 1));
      if ($max_bad_count == 0) {
	last if($g1 ne MISSING_DATA  and  $g2 ne MISSING_DATA  and  $g2 ne $g1 );
	$n_ok_characters++;
      } else {			#
	$bad_count++ if ($g1 ne MISSING_DATA  and  $g2 ne MISSING_DATA  and  $g2 ne $g1 );
	if ($bad_count > $max_bad_count) {
	  last;
	} else {
	  $n_ok_characters++;
	}
      }
    }
  } else {  #  using Inline::C
    # $n_ok_characters in the following means the number of characters before $bad_count reaches $max_bad_count; it includes missing data
    # n_ok_data_chars
    count_oks_and_mismatches_up_to_limit_C($gt1s, $gt2s, $max_bad_count, # $n_gts_above, $n_good_gts_above,
					   $n_ok_characters, $n_good_pairs_to_limit, $bad_count);
    print STDERR "$id2   $L1 $L2   $n_gts_above $n_ok_characters ", $n_ok_characters + $n_gts_above,
      "   $n_good_gts_above $n_good_pairs_to_limit ", $n_good_pairs_to_limit + $n_good_gts_above,
      "   $max_bad_count  $bad_count ", $max_bad_count - $bad_count, "   ", join(',', @{$self->ids()}), "\n";
  }

  if ($n_ok_characters < $L1) { # this branch is ruled out.
    print STDERR "# n_ok_characters ( $n_ok_characters ) < L ( $L1 ). branch ruled out   ids:  ", join(',', @{$self->ids()}), "\n";
  } else {
    assert ($n_ok_characters == $L1) if DEBUG; # ok so far, now examine children if any.
    if ($L2 == $L1) {	 # moment of truth - this must be a leaf node.
      print STDERR "$id2  :  ", join(',', @{$self->ids()}) . "   $max_bad_count $bad_count  ", $max_bad_count - $bad_count, "\n";
      # these are identical genotypes (or rather differ by fewer than $max_bad_count!
      my $this_node_id_string = join(',', @{$self->ids()});
      my $n_to_spare = $max_bad_count - $bad_count;
   #   print "qid:  $id2   matching ids:  ", join(',', @{$self->ids()}) , "\n";
      for my $mid (@{$self->ids()}){
#	print "mid: $mid \n";
	$matchid_matchinfo->{$mid} = $n_to_spare;
      }
      $matching_id_string .= $this_node_id_string . ',';
    } elsif ($L2 > $L1) {	# Ok so far, but must look further ...
      substr($gt2s, 0, $n_ok_characters, '');
      my $g2head = substr($gt2s, 0, 1);
      print STDERR "query gt2s: $gt2s  g2head: $g2head \n";
      if ($g2head eq MISSING_DATA  or  $bad_count < $max_bad_count) { # must check all subtrees
	while ( my ($gh, $child) = each %{ $self->children() }) {
	  print STDERR "  br1. $gh ", join(',', @{$child->ids()}), "\n";
	  $matching_id_string .= $child->search_recursive($id2, $gt2s, $max_bad_count - $bad_count, $matchid_matchinfo,
							  $n_gts_above+$n_ok_characters, $n_good_gts_above+$n_good_pairs_to_limit);
	}
      } else { # $g2head is 0,1, or 2; and $bad_count == $max_bad_count
	for my $gh ($g2head, MISSING_DATA) {
	  my $child = $self->children()->{$gh} // undef;
	  if (defined $child) {
	    print STDERR " br2. $gh ", join(',', @{$child->ids()}), "\n";
	    if ($gh eq $g2head  or  $gh eq MISSING_DATA) {
	      $matching_id_string .= $child->search_recursive($id2, $gt2s, $max_bad_count - $bad_count, $matchid_matchinfo,
							      $n_gts_above+$n_ok_characters, $n_good_gts_above+$n_good_pairs_to_limit);
	    }
	  }else{
	    print STDERR "bad_count == max_bad_count, and no child with head eq gh $gh  ($g2head MISSING_DATA) \n";
	  }
	}
      }
    }
  }
  return $matching_id_string;
}


sub check_node{
  my $self = shift;
  my $error = 0;
  my @check_ids = ();
  if (scalar %{$self->children()} > 0) {
    while (my($g, $chobj) = each %{$self->children()}) {
      $error += 1 if($chobj->parent() ne $self);
      $error += 10 if($self->depth() + length $chobj->genotype() != $chobj->depth());
      push @check_ids,  @{$chobj->ids()};
    }
    my $ids_str = join(',', sort  @{$self->ids()});
    my $check_ids_str = join(',', sort  @check_ids);
    #    print "ids strs in check_node: [$ids_str]  [$check_ids_str] \n";
    $error += 100 if ($check_ids_str ne $ids_str);
  }
  return $error;
}

sub newick_recursive{
  my $self = shift;
  my $newick_str = '';

  assert($self->check_node() == 0) if DEBUG;

  if (keys %{ $self->children() } ) { # not a leaf node.
    $newick_str .= '(';
    for my $g (0,1,2,MISSING_DATA) {
      $newick_str .= ( exists $self->children()->{$g} )? $self->children()->{$g}->newick_recursive() . ',' : '';
    }
    $newick_str =~ s/,$//;
    $newick_str .=
      #')' . $self->genotype() . ':1';
      ')' .
      #	  join('-', @{$self->ids()}) . '_' .
      $self->depth() . '_' . $self->genotype() . ':1';
  } else {			# leaf
    my $gstring = $self->genotype();
    substr($gstring, 4, -4, '...') if(length $gstring > 50); # if long, only show beginning and end
    my $depth = $self->depth();
    #	print "XXX ", join('-', @{$self->ids()}), "  $depth\n"; # if($depth != 200);
    $newick_str .= join('-', @{$self->ids()}) . '_' . $self->depth() . '_' . $gstring . ':1';
  }
  return $newick_str;
}

sub as_string{
  my $self = shift;
  my $leaves_only = shift;
  my $string = '';
  my @child_gs = keys  %{$self->children()};
  if ($leaves_only) {
    if (scalar @child_gs == 0) {
      $string .= 'd: ' . $self->depth() . '  ';
      $string .= 'g: ' . $self->genotype() . '  ';
      $string .=  $self->count() . ' ';
      $string .=  join(',', @{$self->ids()}) . "\n";
    }
  } else {
    $string .= 'depth: ' . $self->depth() . '  ';
    $string .= 'count: ' . $self->count() . '  ';
    $string .= 'ids: ' . join(',', @{$self->ids()}) . '  ';
    $string .= 'genotype: ' . $self->genotype() . '  ';
    $string .= 'children: ' . join(',', @child_gs) . "\n";
  }
  return $string;
}

sub as_string_recursive{
  my $self = shift;
  my $leaves_only = shift;
  my $string = $self->as_string($leaves_only);
  my @gs = sort keys %{$self->children()};
  #  print "child genotypes: ", join('  ', @gs), "\n"; 
  for my $g (@gs) {
    my $child_node = $self->children()->{$g};
    #   print "g: $g ", $child_node->genotype(), "  depths: ", $self->depth(), "  ", $child_node->depth(), "\n";
    $string .= $child_node->as_string_recursive($leaves_only);
  }
  return $string;
}

sub leaves_as_string{
  my $self = shift;
  my $string = '';
  my @child_gs = keys  %{$self->children()};
  if (scalar @child_gs == 0) {
    $string .= 'd: ' . $self->depth() . '  ';
    $string .= 'g: ' . $self->genotype() . '  ';
    $string .=  $self->count() . ' ';
    $string .=  join(',', @{$self->ids()}) . "\n";
  }
  return $string;
}

sub leaves_as_string_recursive{
  my $self = shift;
  my $string = $self->leaves_as_string();
  my @gs = sort keys %{$self->children()};
  #  print "child genotypes: ", join('  ', @gs), "\n"; 
  for my $g (@gs) {
    my $child_node = $self->children()->{$g};
    #   print "g: $g ", $child_node->genotype(), "  depths: ", $self->depth(), "  ", $child_node->depth(), "\n";
    $string .= $child_node-> leaves_as_string_recursive();
  }
  return $string;
}

sub add_id{
  my $self = shift;
  my $id_to_add = shift;
  push @{$self->ids()}, $id_to_add;
}

sub id_as_string{
  my $self = shift;
  return join(',', @{$self->ids()});
}

############################################

__PACKAGE__->meta->make_immutable;

1;

# ############################################
# ########### inline C stuff #################

__DATA__

__C__

#define MISSING_DATA 'X'

#  /* number of mismatches between two strings up to some max  */
# /* returns 0,1,2,...,max_mismatches or -1 if there are more mismatches */
  void count_oks_and_mismatches_up_to_limit_C(char* str1, char* str2, int max_mismatches, SV* n_before_limit, SV* n_good_to_limit, SV* n_mm){ //  SV* n_good_to_limit, SV* n_mm) {
  int i = 0; // n_bef_lim;
  int n_good_pair = 0; // n_good2lim; // count the number of good pairs (i.e. with neither being MISSING_DATA
  char c1;
  char c2;
  int mismatch_count = 0;
// printf("%d %d %d \n", max_mismatches, n_bef_lim, n_good2lim);
  while(c1 = str1[i]) {
    if(c1 != MISSING_DATA){
      c2 = str2[i];
      if(c2 != MISSING_DATA){
	n_good_pair++;
        if (c1 != c2) {
          mismatch_count++;
        }
        if(mismatch_count > max_mismatches){ // i not incremented
        //  mismatch_count = -1;
          break;
        }
      }
    }
    i++; //  i is now the number of characters which have been compared without exceeding the mismatch limit. 
  } // end of while
// printf("%s %s  %d %d %c\n", str1, str2, i, mismatch_count, MISSINGDATA);
							  sv_setiv(n_before_limit, i);
							  sv_setiv(n_good_to_limit, n_good_pair);
							sv_setiv(n_mm, mismatch_count);
}

   /* number of mismatches between two strings up to some max  */
/* returns 0,1,2,...,max_mismatches or -1 if there are more mismatches */
/*
    void count_homozygous_oks_and_mismatches_up_to_limit_C(char* str1, char* str2, int max_mismatches, SV* n_ok, SV* n_mm) {
  int i = 0;
  int n_both_homozyg = 0;
  char c1;
  char c2;
  int mismatch_count = 0;
int zero2_20_count = 0;
int both_homozygous_count = 0;
  while(c1 = str1[i]) {
    if(c1 != MISSING_DATA  and  c1 != '1'){
      c2 = str2[i];
      if(c2 != MISSING_DATA  and  c2 != '1'){
        n_both_homozyg++;
        if (c1 != c2) {
          mismatch_count++;
        }
        if(mismatch_count > max_mismatches){
          mismatch_count = -1;
          break;
        }
      }
    }
    i++;
//  i is now the number of characters which have been compared without exceeding the mismatch limit. 
  }
 printf("%s %s  %d %d %c\n", str1, str2, i, mismatch_count, MISSINGDATA);
  sv_setiv(n_ok, i);
  // sv_setiv(n_n
  sv_setiv(n_mm, mismatch_count);
}

*/
