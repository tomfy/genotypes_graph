package GenotypeGraph;
use Moose;
use namespace::autoclean;
use Carp;
use List::Util qw ( min max sum );
# use Graph;
use MyPriorityQueue;

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

has n_near => (
                isa => 'Int', # the number of outward edges to keep in the graph for each node.
                is => 'ro',
                default => BIG_NUMBER, # -> keep all edges.
               );
has n_extras =>  (
                isa => 'Int', # the number randomly chosen nodes to add as neighbors in addition to near ones.
                is => 'ro',
                default => 0, # just keep the closest ones.
               );

has distances => ( # initially all N choose 2 distances (edge weights)
                  isa => 'HashRef',
                  is => 'rw',
                  default => sub { {} },
                 );

around BUILDARGS => sub {
  my $orig = shift;
  my $class = shift;

  my $args = shift; # hash ref, e.g. { fasta => <a fasta string for several sequences> , n_near => 5 }
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
	my $d = $g1->distance($g2);
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
    my $n_near = $args->{n_near} // BIG_NUMBER;
    $n_near = min($n_near, $n_nodes-1);
    my $n_extras = $args->{n_extras};

    while ( my ($i, $id1) = each @ids) {
      my $id2_dist = $idA__idB_distance{$id1};
      my $count = 0;
      my @neighbor_ids;
      if (1) {		  # using quickselect algorithm (a bit faster)
	@neighbor_ids = quickselect([keys %$id2_dist], $n_near, $id2_dist);
	@neighbor_ids = sort { $id2_dist->{$a} <=> $id2_dist->{$b} } @neighbor_ids;
      } else {
	@neighbor_ids = sort { $id2_dist->{$a} <=> $id2_dist->{$b} } keys %$id2_dist;
	@neighbor_ids = @neighbor_ids[0 .. $n_near-1] if($n_near < scalar keys %$id2_dist);
      }
      my %neighborid_dist = map(($_ => $id2_dist->{$_}), @neighbor_ids); # hash w ids, distances for just nearest $n_near
      get_extra_ids($id2_dist, \@neighbor_ids, \%neighborid_dist, $n_extras);

      $id_node->{$id1} = GenotypeGraphNode->new( {
					    id => $id1,
					    genotype => $id_gobj{$id1}, # genotype object
					    neighbor_ids => \@neighbor_ids,
					    neighbor_id_distance => \%neighborid_dist,
					   } );

    } # end node construction loop
    return {
	    nodes => $id_node,
	    n_near => $n_near,
	    n_extras => $n_extras,
	    distances => \%idA__idB_distance,
	   };
  } elsif (defined $args->{idnnd}) {
    # construct from file with info on nearest neighbors to each node. typical line:
    # 290  2  ((()26,()24)218,(()52,()22)176)290    176  0.383   218  0.42   274  0.482   354  0.501   52  0.502   ...   328  0.757
    my $n_near = BIG_NUMBER;
    my %idA__idB_distance = (); # hash of hash refs
    my $id_node = {};
    my @graph_string_lines = split("\n", $args->{idnnd});
    #  my $multiplier = 1;
    for my $line (@graph_string_lines) {
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
      my $neighborid_dist = ();
      my @neighbor_ids = ();
      while ($line) {
	$line =~ s/\s*(\S+)\s+(\S+)\s*//;
	my ($id2, $dist) = ($1, $2);
	$neighborid_dist->{$id2} = $dist;
	push @neighbor_ids, $id2;
      }
      $n_near = min($n_near, scalar @neighbor_ids);
      my $gobj = Genotype->new( '>' . "$id1 $generation $pedigree\n" . $sequence );
      print "gobj as string:  ", $gobj->id(), "  ", join('', @{$gobj->sequence()}), "\n";
      my $a_node = GenotypeGraphNode->new( {
					    id => $id1,
					    genotype => $gobj, # genotype object
					    neighbor_ids => \@neighbor_ids,
					    neighbor_id_distance => $neighborid_dist,
					    #           furthest_id_distance => [$furthest_id, $furthest_d],
					   } );
      $id_node->{$id1} = $a_node;
    }
    return {
	    nodes => $id_node,
	    n_near => $n_near,
	    n_extras => undef, 
	    distances => \%idA__idB_distance,
	   };
  }				# end of

  # otherwise do nothing.
};                              # end of BUILDARGS

sub BUILD {
  my $self = shift;
  my @nodes = values %{$self->nodes()};
  while (my($id, $node_obj) = each %{$self->nodes()} ) {
    #  $node_obj->graph($self); # can't do this because graph is read-only, but
    $node_obj->{graph} = $self; # can do it this way.
  }

}

sub get_node_by_id{
  my $self = shift;
  my $id = shift;
  return $self->nodes()->{$id};
}

sub search_for_best_match{
  my $self = shift;
  my $gobj = shift;		# Genotype object
  #  my $rand_id = {$self->nodes()};

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

####   ordinary subroutines  ###

sub get_extra_ids{
  my $id2_dist = shift;
  my $near_ids = shift;
  my $nearid_dist = shift;
  my $n_extras = shift;
  my $n_try = shift // 5;
  for (1..$n_extras) {
    for (1..$n_try) {
      my $extra_id = (keys %$id2_dist)[int(rand keys %$id2_dist)];
#	(keys %$id2_dist)[0 - $_];
      if (! exists $nearid_dist->{$extra_id}) {
	push @$near_ids, $extra_id;
	$nearid_dist->{$extra_id} = $id2_dist->{$extra_id};
	last;
      }
    }
  }
}

sub quickselect{		# get the nearest $k ids.
  my $id_list = shift;
  my $k = shift;	      # find this many.
  my $id_distance = shift;    # hash ref of all N-1 id:distance pairs.

  my $rand_index = int rand @$id_list;
  my $pivot_id = $id_list->[ int(rand(@$id_list)) ];
  my $pivot = $id_distance->{$pivot_id};

  my @lefts = (); my @rights = (); my @equals = ();
  if (1) { # use built-in grep, but run through array 2 (or sometimes 3) times; a bit slower.
    @lefts  = grep { $id_distance->{$_} < $pivot } @$id_list;
    @rights = grep { $id_distance->{$_} > $pivot } @$id_list;
    if (@lefts + @rights + 1 == scalar @$id_list) {
      push @equals, $pivot_id;
    } else {
      @equals = grep { $id_distance->{$_} == $pivot } @$id_list;
    }
  } else {			# a bit faster
    for my $an_id (@$id_list) { # store in separate arrays the ids with distance <, ==, and > the pivot 
      if ($id_distance->{$an_id} < $pivot) {
	push @lefts, $an_id;
      } elsif ($id_distance->{$an_id} > $pivot) {
	push @rights, $an_id;
      } else {
	push @equals, $an_id;
      }
    }
  }

  if ($k < @lefts) {  # kth will be in @lefts, but lefts has too many.
    return quickselect(\@lefts, $k, $id_distance);
  } elsif ($k > @lefts + @equals) { # kth will be in @rights
    return ((@lefts, @equals), quickselect(\@rights, $k - @lefts - @equals, $id_distance));
  } elsif ($k == @lefts) {	# done 
    return @lefts;
  } elsif ($k <= @lefts + @equals) { # just @lefts plus 1 or more from @equals
    push @lefts, @equals[0..$k-@lefts-1];
    return @lefts;
  } else {
    die "???\n";
  }
}

###########################################
__PACKAGE__->meta->make_immutable;

1;


