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

# class to represent a tree of genotypes

sub BUILD{
  my $self = shift;
  $self->root()->tree($self);
}

sub add_genotype{
  my $self = shift;
  my $gobj = shift ;
  my $genotype = $gobj->sequence(); # ArrayRef, e.g. [0, 0, 1, 2, 1, 2, 2, 0, 0, 0, 1, 0, '-', 1, 0];
  my $id = $gobj->id();
  my $root = $self->root();
  # $root->inc_counter();
  $root->{count}++;
  my @ids_array = push @{$root->ids()}, $id;
  $root->ids( \@ids_array );
  my $current_node = $self->root();
  for my $g (@$genotype) {
    my $next_node =  $current_node->add_child($id, $g);
    $current_node = $next_node;
  }
  return $current_node->depth() . " " . join(',', $current_node->ids());
}

sub add_genotype_compact{
  my $self = shift;
  my $gobj = shift;
  my $genotype_string = $gobj->sequence(); # entire (all snps) genotype as string.
  my $id = $gobj->id();
  my $root = $self->root();
  $root->add_id($id);
  my $ghead = substr($genotype_string, 0, 1);
  if (exists $root->children()->{$ghead}) {
    $root->children()->{$ghead}->add_genotype_compact($id, $genotype_string);
  } else {
    my $new_node = GenotypeTreeNode->new( {tree => $self,
					   parent => $root,
					   depth => length $genotype_string,
					   genotype => $genotype_string, 
					   ids => [$id] } );
    $root->children()->{$ghead} = $new_node;
  }
  $self->size($self->size()+1);
  #  print "done adding genotype.\n\n";
  # exit if($self->size() > 20);
}

sub search{
  my $self = shift;
  my $gobj = shift;
  my $root = $self->root();
  my $genotype_string = $gobj->sequence();
  my $ghead = substr($genotype_string, 0, 1);
  my $matching_ids = '';
  if ($ghead eq MISSING_DATA) {
    while (my($gh, $child) = each %{$root->children()}) { # search all children of root if 1st char is missing data.
      $matching_ids .= $child->search_recursive($gobj->id(), $gobj->sequence());
    }
  } else {
    while (my($gh, $child) = each %{$root->children()}) { # search all children of root if 1st char is missing data.
      if ($gh eq $ghead  or  $gh eq MISSING_DATA) {
	$matching_ids .= $child->search_recursive($gobj->id(), $gobj->sequence());
      }
    }
  }
  $matching_ids =~ s/,\s*$//;
  my @sorted_ids = sort {$a <=> $b} split(',', $matching_ids);

  return( scalar @sorted_ids > 0)?
    # $gobj->id() . "  " .
    join(",", @sorted_ids) : '';
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


