package GenotypeGraph;
use Moose;
use namespace::autoclean;
use Carp;
use List::Util qw ( min max sum );
use Graph;

use constant BIG_NUMBER => 1_000_000_000;
#use constant MULTIPLIER => 1000;

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

has distances => (
                  isa => 'HashRef',
                  is => 'rw',
                  default => sub { {} },
                 );

has Graph_object => (
		     isa => 'Object',
		     is => 'rw',
		     default => sub { Graph->new() },
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

      # calculate and store distances (N choose 2 of them)
      my %idA__idB_distance = (); # hash of hash refs
      my @ids = sort { $a <=> $b } keys %id_gobj;
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

      # for each genotype, find the nearest neighbors and construct node.
      my $id_node = {};
      my $n_edges = $args->{n_edges} // BIG_NUMBER;
      $n_edges = min($n_edges, $n_nodes-1);

      while ( my ($i, $id1) = each @ids) {
	my $id2_dist = $idA__idB_distance{$id1};
	my $count = 0;
	if(1){
	   my $pq = List::PriorityQueue->new();
	  while(my($id2, $d) = each %$id2_dist){
	    $pq->insert($id2, $dist);
	    $count++;
	    if($count >= $n_edges){
	      $pq->worst(); # shift the 'worst' element, to keep size <= $n_edges.
	    }
	  }
	    
	}else{
         my @nearest_neighbor_ids = sort { $id2_dist->{$a} <=> $id2_dist->{$b} } keys %$id2_dist;
         my $furthest_id = $nearest_neighbor_ids[-1]; # last one = furthest
         my $furthest_d = $id2_dist->{$furthest_id};
         @nearest_neighbor_ids = @nearest_neighbor_ids[0 .. $n_edges-1] if($n_edges < scalar keys %$id2_dist);
         my $a_node = GenotypeGraphNode->new( {
                                               id => $id1,
                                               genotype => $id_gobj{$id1}, # genotype object
                                               nearest_neighbor_ids => \@nearest_neighbor_ids,
                                               id_distance => $id2_dist,
                                               furthest_id_distance => [$furthest_id, $id2_dist->{$furthest_id}],
                                              } );
         $id_node->{$id1} = $a_node;
       }
      }                         # end node construction loop
      return {
              nodes => $id_node,
              n_edges => $n_edges,
              distances => \%idA__idB_distance,
             };
   } elsif (defined $args->{idnnd}) {
      # construct from file with info on nearest neighbors to each node. typical line:
      # 290  2  ((()26,()24)218,(()52,()22)176)290    176  0.383   218  0.42   274  0.482   354  0.501   52  0.502   ...   328  0.757
      my $n_edges = BIG_NUMBER;
      my %idA__idB_distance = (); # hash of hash refs
      my $id_node = {};
      my @graph_string_lines = split("\n", $args->{idnnd});
    #  my $multiplier = 1;
      for my $line (@graph_string_lines){
         next if($line =~ /^\s*#/);
         # $multiplier = $1 if(/#\s+(\S+)/;
         $line =~ s/^\s*(\S+)\s+(\S+)\s+(\S+)\s+//;
         print "XXX $line \n";
         my ($id1, $generation, $pedigree) = ($1, $2, $3);
print "$id1, $generation, $pedigree \n";
         $line =~ s/\.{3}\s+(\S+)\s+(\S+)/.../;
         my ($furthest_id, $furthest_d) = ($1, $2);
         print "$furthest_id  $furthest_d \n";
         my $sequence = ($line =~ /\s*\.{3}\s+(\S+)\s*$/)? $1 : '---';
         $line =~ s/\s*\.{3}.*$//;
         my $id2_dist = ();
         my @nearest_neighbor_ids = ();
         while ($line) {
            $line =~ s/\s*(\S+)\s+(\S+)\s*//;
            my ($id2, $dist) = ($1, $2);
            $id2_dist->{$id2} = $dist;
            push @nearest_neighbor_ids, $id2;
         }
         $n_edges = min($n_edges, scalar @nearest_neighbor_ids);
         my $gobj = Genotype->new( '>' . "$id1 $generation $pedigree\n" . $sequence );
         print "gobj as string:  ", $gobj->id(), "  ", join('', @{$gobj->sequence()}), "\n";
         my $a_node = GenotypeGraphNode->new( {
                                               id => $id1,
                                               genotype => $gobj, # genotype object
                                               nearest_neighbor_ids => \@nearest_neighbor_ids,
                                               id_distance => $id2_dist,
                                               furthest_id_distance => [$furthest_id, $furthest_d],
                                              } );
         $id_node->{$id1} = $a_node;
      }
      return {
              nodes => $id_node,
              n_edges => $n_edges,
              distances => \%idA__idB_distance,
             };
   }                            # end of

  # otherwise do nothing.
};                              # end of BUILDARGS


sub get_node_by_id{
   my $self = shift;
   my $id = shift;
   return $self->nodes()->{$id};
}

sub as_string{
   my $self = shift;
   my $show_sequences = shift;
   my $str = '';
   my @sorted_ids = sort { $a <=> $b } keys %{ $self->nodes() };
   #   while (my ($node_id, $node_obj) = each %{ $self->nodes() } ) {
   for my $id (@sorted_ids) {
      $str .= $self->get_node_by_id($id)->as_string($show_sequences) . "\n";
   }
   return $str;
}

sub distance_matrix_as_string{
   my $self = shift;
   my $multiplier = shift;
   my $idA__idB_distance = $self->distances();
   my @ids = sort {$a <=> $b} keys %$idA__idB_distance;
   my $n_nodes = scalar @ids;
   my $d_matrix_string = '# ' . $multiplier ."\n";
   $d_matrix_string .= join(" ", @ids) . "\n";
   while ( my ($i, $id1) = each @ids) {
      my $id2_dist = $idA__idB_distance->{$id1};
      my @id2s = sort { $a <=> $b } keys %$id2_dist;
      $d_matrix_string .= sprintf("%2d  %s\n", $id1, join(" ", map (int($multiplier*$id2_dist->{$_} + 0.5), @id2s[$i..$n_nodes-2]) ) );
   }
   return $d_matrix_string;
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


