package GenotypeTree;
#use Moose;
use Mouse;
use namespace::autoclean;
use Carp;
use List::Util qw ( min max sum );
use constant MISSING_DATA => 'X';

has root => (
             isa => 'Object',   # should be a GenotypeTreeNode Object
             is => 'ro',
             default => sub { GenotypeTreeNode->new( { depth => 0, genotype => 'R' } ); },
            );

has depth => (                  # equal to length of sequences
              isa => 'Int',
              is => 'ro',
              required => 1,
             );

has size => (		  # number of genotypes which have been added.
	     isa => 'Int',
	     is => 'rw',
	     default => 0,
	    );

has count_search_recursive_calls => (
                                     isa => 'Int',
                                     is => 'ro',
                                     default => 0,
				    );

# class to represent a tree of genotypes

sub BUILD{
  my $self = shift;
  $self->root()->tree($self);
}

sub add_genotype{
  my $self = shift;
  my $gobj = shift;
  my $chunk_indices = shift // undef;
  my $genotype_string = $gobj->get_chunk($chunk_indices); #
  my $id = $gobj->id();
  my $root = $self->root();
  $root->add_id($id);
  my $ghead = substr($genotype_string, 0, 1);
  if (exists $root->children()->{$ghead}) {
    $root->children()->{$ghead}->add_genotype_recursive($id, $genotype_string);
  } else {
    my $new_node = GenotypeTreeNode->new( {tree => $self,
					   parent => $root,
					   depth => length $genotype_string,
					   genotype => $genotype_string,
					   ids => [$id] } );
    $root->children()->{$ghead} = $new_node;
  }
  $self->size($self->size()+1);
}

sub search{
  my $self = shift;
  my $gobj = shift;
  my $chunk_indices = shift // undef;
  my $max_bad_count = shift // die;
  my $min_to_store = shift; # if number of pairs with no missing data when max_mismatch_count is reached is < this, do not store info
  my $mid_matchinfosum = shift;
  my $root = $self->root();
  my $genotype_string = $gobj->get_chunk($chunk_indices); # query genotype string

  my $ghead = substr($genotype_string, 0, 1);
  my $matchid_matchinfo = {};
  my $n_gts_above = 0;
  my $n_good_gts_above;

  while (my($gh, $child) = each %{$root->children()}) {	# search all children of root if 1st char is missing data.
    $child->search_recursive(
			     $genotype_string,
			     $max_bad_count, $min_to_store, 0,
			     0, 0, $mid_matchinfosum);
  }
}

sub as_newick{
  my $self = shift;
  return $self->root()->newick_recursive();
}

sub as_string{
  my $self = shift;
  my $leaves_only = shift // 0;
  my $string = ($leaves_only)? 
    $self->root()->leaves_as_string_recursive($leaves_only) : 
    $self->root()->as_string_recursive($leaves_only);
  return $string;
}



###########################################
__PACKAGE__->meta->make_immutable;

1;


# sub add_genotype_old{
#   my $self = shift;
#   my $gobj = shift ;
#   my $genotype = $gobj->sequence(); # ArrayRef, e.g. [0, 0, 1, 2, 1, 2, 2, 0, 0, 0, 1, 0, '-', 1, 0];
#   my $id = $gobj->id();
#   my $root = $self->root();
#   # $root->inc_counter();
#   $root->{count}++;
#   my @ids_array = push @{$root->ids()}, $id;
#   $root->ids( \@ids_array );
#   my $current_node = $self->root();
#   for my $g (@$genotype) {
#     my $next_node =  $current_node->add_child($id, $g);
#     $current_node = $next_node;
#   }
#   return $current_node->depth() . " " . join(',', $current_node->ids());
# }
