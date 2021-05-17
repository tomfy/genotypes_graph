package GenotypeTreeNode;
#use Moose;
use Mouse;
use namespace::autoclean;
use Carp::Assert;
use List::Util qw ( min max sum );
use constant MISSING_DATA => 'X';
use Inline C => '../lib/inlinec.c';

no warnings 'recursion';

# use constant DEBUG => 0;

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
  # my $id2 = shift; # not used
  my $gt2s = shift;		# searching for this one.
  my $max_bad_count = shift // die; # the max number of mismatching characters allowed along this branch
  my $min_to_store = shift; #
  my $bad_count = shift;       # number of mismatches upstream of this
  my $matchid_matchinfo = shift; # key match id, values: string: [] (n_to_spare is $max_bad_count - $bad_count)
  my $n_gts_above = shift;
  my $n_good_gts_above = shift;
  my $gt1s = $self->genotype(); # the genotype along the branch above this node.
  my ($L1, $L2) = (length $gt1s, length $gt2s); # $L1: gts on this branch (above node); $L2 gts left in query.
  $self->{tree}->{count_search_recursive_calls}++;
  #my $min_to_store = 5; # 3*($max_bad_count+1);

  assert ($L2 >= $L1) if DEBUG;
  assert ($L2 == $L1  or (keys %{$self->children()} > 0)) if DEBUG;

  my $n_ok_characters = 0; # the number of characters before the bad_count gets too big.
  my $n_good_pairs_to_limit = 0; # number of gt pairs (neither being MISSING_DATA) to exceed max_bad
  #  using Inline::C
  # $n_ok_characters in the following means the number of characters before $bad_count reaches $max_bad_count; it includes missing data

  if (0) { # using perl count mismatches between strings
    ($n_ok_characters, $n_good_pairs_to_limit, $bad_count) =
      count_oks_and_mismatches_up_to_limit_perl($gt1s, $gt2s, $max_bad_count, $bad_count);
  } else { # using C count mismatches between strings
    count_oks_and_mismatches_up_to_limit_C($gt1s, $gt2s, $max_bad_count, $bad_count,
					   $n_ok_characters, $n_good_pairs_to_limit, $bad_count);
  }
  # print STDERR "qid: $id2   L1,L2: $L1 $L2",
  #  "   $n_gts_above $n_ok_characters ", $n_ok_characters + $n_gts_above,
  #  "   $n_good_gts_above $n_good_pairs_to_limit ", $n_good_pairs_to_limit + $n_good_gts_above,
  #  "   max_bad_count, bad_count: $max_bad_count  $bad_count ", $max_bad_count - $bad_count, "   ", join(',', @{$self->ids()}), "\n";

  if ($n_ok_characters < $L1) { # this branch is ruled out.
    if ( (my $n_good_to_exceed_max_bad_count = $n_good_gts_above + $n_good_pairs_to_limit) >= $min_to_store) {
      for my $mid (@{$self->ids()}) {
	$matchid_matchinfo->{$mid} = "$bad_count  " . ($n_gts_above+$n_ok_characters) . "  $n_good_to_exceed_max_bad_count";
      }
    }
  } else {
    assert ($n_ok_characters == $L1) if DEBUG; # ok so far, now examine children if any.
    if ($L2 == $L1) {	# moment of truth - this must be a leaf node.
      # these a differ by fewer than $max_bad_count!
      if ( (my $n_good_to_exceed_max_bad_count = $n_good_gts_above + $n_good_pairs_to_limit) >= $min_to_store) {
	for my $mid (@{$self->ids()}) {
	  $matchid_matchinfo->{$mid} = "$bad_count  " . ($n_gts_above+$n_ok_characters) . "  $n_good_to_exceed_max_bad_count";
	}
      }
    } elsif ($L2 > $L1) {	# Ok so far, but must look further ...
      substr($gt2s, 0, $n_ok_characters, ''); 
      my $g2head = substr($gt2s, 0, 1);
      if (1  or  $g2head eq MISSING_DATA  or  $bad_count < $max_bad_count) { # must check all subtrees
	while ( my ($gh, $child) = each %{ $self->children() }) {
	  $child->search_recursive($gt2s, $max_bad_count, $min_to_store, $bad_count, $matchid_matchinfo,
				   $n_gts_above+$n_ok_characters, $n_good_gts_above+$n_good_pairs_to_limit);
	}
      } else { # $g2head is 0,1, or 2; and $bad_count == $max_bad_count
	for my $gh ($g2head, MISSING_DATA) {
	  my $child = $self->children()->{$gh} // undef;
	  if (defined $child) {
	    if ($gh eq $g2head  or  $gh eq MISSING_DATA) {
	      $child->search_recursive($gt2s, $max_bad_count, $min_to_store, $bad_count, $matchid_matchinfo,
				       $n_gts_above+$n_ok_characters, $n_good_gts_above+$n_good_pairs_to_limit);
	    }
	  } else {	    # in this case (bad_count == max_bad_count
	    #print STDERR "bad_count == max_bad_count, and no child with head eq gh $gh  ($g2head MISSING_DATA) \n";
	  }
	}
      }
    }
  }
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

sub count_oks_and_mismatches_up_to_limit_perl{
  my $str1 = shift;
  my $str2 = shift;
  my $max_mismatches = shift;
  my $mismatches_so_far = shift;
#  my $n_before_limit = shift;
#  my $n_good_to_limit = shift;
  my $n_good_pair = 0;
  my $mismatch_count = $mismatches_so_far;
  my $i = 0;
  for( ; $i< length $str1; $i++){
    my $c1 = substr($str1, $i, 1);
    if($c1 ne MISSING_DATA){
      my $c2 = substr($str2, $i, 1);
      if($c2 ne MISSING_DATA){
	$n_good_pair++;
	if($c1 ne $c2){
	  $mismatch_count++;
	}
	if($mismatch_count >= $max_mismatches){
	  last;
	}
      }
    }
    
  }
  return ($i, $n_good_pair, $mismatch_count);
}
############################################

__PACKAGE__->meta->make_immutable;

1;

# ############################################
# ########### inline C stuff #################
# ###### is now in ../lib/inlinec.c ##########
# ############################################
