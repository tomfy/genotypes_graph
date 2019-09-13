package GenotypeGraphNode;
#use Moose;
use Mouse;
use namespace::autoclean;
use Carp;
use List::Util qw ( min max sum );


# a class for nodes of a (directed) graph.
# each node has an id number, and a hash of
# id/distance pairs holding the ids and distances to
# other nodes to which there is a directed edge.

has graph => (	 # the GenotypeGraph object to which this node belongs
              isa => 'Object',
              is => 'ro',
              required => 0,
	     );

has id => (
           isa => 'Int',
           is => 'ro',
           required => 1,
          );

# has neighbor_ids => (  # just use keys of neighbor_id_distance
# 		     isa => 'ArrayRef',
# 		     is => 'ro',
# 		     required => 0,
# 		    );

has neighbor_id_distance => (		# should be just the neighbors. 
                    isa => 'HashRef',
                    is => 'rw',
                    required => 1,
                   );

has genotype => (
                 isa => 'Object',
                 is => 'ro',
                 required => 1,
		);

sub as_string{
  my $self = shift;
  my $show_sequence = shift // 0;
  my $str = $self->id() . '  ' . $self->genotype()->generation() . '  ' . $self->genotype()->pedigree() . '   ';
  my @nn_ids =
    #    sort { $a <=> $b }  # uncomment to sort by id
    keys %{$self->neighbor_id_distance()};
    #    @{$self->neighbor_ids()};	# sort by id
  
  for my $idB (@nn_ids) {
    my $dist = $self->neighbor_id_distance->{$idB};
    $str .= sprintf("%4d %5.4f  ", $idB, $dist);
  }
  $str .= '   ' . join('', @{$self->genotype()->sequence()} ) if($show_sequence);
  return $str;
}

sub symmetrize_neighbors{
  my $self = shift;
  # make neighbors of self have self as neighbors
  my ($graph, $self_id) = ($self->graph(), $self->id());
  while(my($nid, $d) = each %{$self->neighbor_id_distance()}){ # distance from $self to $nid is $d
    my $neighbors_of_neighbor = $graph->nodes()->{$nid}->neighbor_id_distance(); # hashref with neighbor ids and dists of $nid
    if(!exists $neighbors_of_neighbor->{$self_id}){ # if $self is not already a neighbor of $nid
      $neighbors_of_neighbor->{$self_id} = $d;
    }
  }
}


############################################

__PACKAGE__->meta->make_immutable;

1;
