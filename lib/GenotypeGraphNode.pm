package GenotypeGraphNode;
use Moose;
use namespace::autoclean;
use Carp;
use List::Util qw ( min max sum );


# a class for nodes of a (directed) graph.
# each node has an id number, and a hash of
# id/distance pairs holding the ids and distances to
# other nodes to which there is a directed edge.

has id => (
           isa => 'Int',
           is => 'ro',
           required => 1,
          );

has nearest_neighbor_ids => (
                             isa => 'ArrayRef',
                             is => 'ro',
                             required => 1,
                            );

has id_distance => (
                    isa => 'HashRef',
                    is => 'rw',
                    required => 1,
                   );

has furthest_id_distance => (
                             isa => 'ArrayRef',
                             is => 'ro',
                             );

has genotype => (
                 isa => 'Object',
                 is => 'ro',
                 required => 1,
                 );


sub BUILD {
   my $self = shift;
   my %nearestid_distance = ();
   for my $idB (@{$self->nearest_neighbor_ids()}) {
 #     print "idB:  $idB \n";
      $nearestid_distance{$idB} = $self->id_distance()->{$idB};
   }
#    my $furthest_id = $self->nearest_neighbor_ids()->[-1];
#    my $furthest_distance = $self->id_distance()->{$furthest_id};
# print STDERR $self->{id}, "    ", $self->furthest_id_distance->[0], "   ",  $self->furthest_id_distance()->[1], " \n";
   $self->id_distance(\%nearestid_distance);
 #  $self->{furthest_id_distance} = [$furthest_id => $furthest_distance];
   return $self;
}

sub as_string{
   my $self = shift;
   my $show_sequence = shift // 0;
   my $str = $self->id() . '  ' . $self->genotype()->generation() . '  ' . $self->genotype()->pedigree() . '   ';
   for my $idB ( @{ $self->nearest_neighbor_ids() } ) {
      my $dist = $self->id_distance->{$idB};
      $str .= sprintf("%6d %5.4f  ", $idB, $dist);
   }
 #  $str .= '   ...   ' . join("  ", @{ $self->furthest_id_distance() });
   $str .= '   ' . join('', @{$self->genotype()->sequence()} ) if($show_sequence);
   return $str;
}






############################################

__PACKAGE__->meta->make_immutable;

1;
