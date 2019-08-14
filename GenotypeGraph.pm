package GenotypeGraph;
use Moose;
use namespace::autoclean;
use Carp;
use List::Util qw ( min max sum );
#use Readonly;

use constant BIG_NUMBER => 1_000_000_000;

# class to represent a graph, i.e. a set of node, which may be
# connected in pairs by edges. Here the nodes know which other nodes
# they are connected to. 

has nodes => (
              isa => 'HashRef', # keys are ids, values GenotypeGraphNode objects.
              is => 'ro',
              default => sub { {} },
             );

has n_edges => (
                isa => 'Int', # the number of outward edges to keep in the graph for each node.
                is => 'ro',
                default => BIG_NUMBER, # -> keep all edges.
               );

around BUILDARGS => sub {
   my $orig = shift;
   my $class = shift;

   my $args = shift; # hash ref, e.g. { fasta => <a fasta string for several sequences> , n_edges => 5 }
   if (defined $args->{fasta} ) { # construct from fasta string
      # print $args->{fasta}, "\n";
      my @fasta_lines = split("\n", $args->{fasta} );
      my %id_gobj = ();
      while (@fasta_lines) {
         my $id_line = shift @fasta_lines;
         my $sequence = shift @fasta_lines;
         my $gobj = Genotype->new( $id_line . "\n" . $sequence );
         $id_gobj{$gobj->id()} = $gobj;
      }

      my %idA__idB_distance = (); # hash of hash refs
      my @ids = keys %id_gobj;
      my $n_nodes = scalar @ids;
      my ($edge_count, $d_sum) = (0, 0);
      while (my ($i1, $id1) = each @ids) {
         #   print STDERR "$i1 $id1 \n";
         my $g1 = $id_gobj{$id1};
         for (my $i2 = $i1+1; $i2 < scalar @ids; $i2++) {
            my $id2 = $ids[$i2];
            #    print STDERR "    $i2 $id2 \n";
            my $g2 = $id_gobj{$id2};
            my ($d_ij, $both_snp_count) = $g1->distance($g2);
            #   print "$id1  $id2  $d_ij \n";
            my $d = ($both_snp_count > 0)? $d_ij/$both_snp_count : BIG_NUMBER; #
            $d_sum += $d;
            $edge_count++;
            if (!exists $idA__idB_distance{$id1}) {
               $idA__idB_distance{$id1} = {$id2 => $d};
            } else {
               $idA__idB_distance{$id1}->{$id2} = $d;
            }
            if (!exists $idA__idB_distance{$id2}) {
               $idA__idB_distance{$id2} = {$id1 => $d};
            } else {
               $idA__idB_distance{$id2}->{$id1} = $d;
            }
         }
      }

      my $id_node = {};
      my $n_edges = $args->{n_edges} // BIG_NUMBER;
      $n_edges = min($n_edges, $n_nodes-1);
      while (my($id1, $id2_dist) = each %idA__idB_distance) { # get nearest neighbors and construct nodes
         my @nearest_neighbor_ids = sort { $id2_dist->{$a} <=> $id2_dist->{$b} } keys %$id2_dist;
         my $furthest_id = $nearest_neighbor_ids[-1]; # last one = furthest
         my $furthest_d = $id2_dist->{$furthest_id};
         #   print STDERR "$id1 ...  $furthest_id  $furthest_d \n";
         @nearest_neighbor_ids = @nearest_neighbor_ids[0 .. $n_edges-1] if($n_edges < scalar keys %$id2_dist);
         my $a_node = GenotypeGraphNode->new( {
                                               id => $id1,
                                               genotype => $id_gobj{$id1}, # genotype object
                                               nearest_neighbor_ids => \@nearest_neighbor_ids,
                                               id_distance => $id2_dist,
                                               furthest_id_distance => [$furthest_id, $id2_dist->{$furthest_id}],
                                              } );
         $id_node->{$id1} = $a_node;
         # print "Node as string:  ", $a_node->as_string(), "\n";
      }  # end node construction loop
      return {
              nodes => $id_node,
              n_edges => $n_edges,
             };
   }
   # otherwise do nothing.
}; # end of BUILDARGS


sub get_node_by_id{
   my $self = shift;
   my $id = shift;
   return $self->nodes()->{$id};
}

sub as_string{
   my $self = shift;
   my $str = '';
   my @sorted_ids = sort { $a <=> $b } keys %{ $self->nodes() };
#   while (my ($node_id, $node_obj) = each %{ $self->nodes() } ) {
   for my $id (@sorted_ids){
      $str .= $self->get_node_by_id($id)->as_string() . "\n";
   }
   return $str;

}



# sub BUILD{
#    my $self = shift;
#    my $args = shift;         # a hash ref {'id12dist' => 
#    my $id_sequence = $args->{'id_sequence'};

#    #########################################################



#    #####################################################

#    my $id12dist = $args->{'id12dist'};
#    #   $self->n_edges($args->{'n_edges'});
#    my $n_edges = $self->n_edges();
#    while (my($id1, $id2_dist) = each %$id12dist) {
#       my @nearest_neighbor_ids = sort { $id2_dist->{$a} <=> $id2_dist->{$b} } keys %$id2_dist;
#       my $furthest_id = $nearest_neighbor_ids[-1]; # last one = furthest
#       my $furthest_d = $id2_dist->{$furthest_id};
#       #   print STDERR "$id1 ...  $furthest_id  $furthest_d \n";
#       @nearest_neighbor_ids = @nearest_neighbor_ids[0 .. $n_edges-1] if($n_edges > 0 and $n_edges < scalar keys %$id2_dist);
#       my $a_node = GenotypeGraphNode->new( {
#                                             id => $id1,
#                                             nearest_neighbor_ids => \@nearest_neighbor_ids,
#                                             id_distance => $id2_dist,
#                                             furthest_id_distance => [$furthest_id => $id2_dist->{$furthest_id}],
#                                            } );
#       $self->nodes()->{$id1} = $a_node;
#       #  print "Node as string:  ", $a_node->as_string(), "\n";
#    }
# }            
# end BUILD





###########################################
__PACKAGE__->meta->make_immutable;

1;


