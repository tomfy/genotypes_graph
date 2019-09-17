package GenotypeGraph;
use strict;
use warnings;
use Moose;
#use Mouse;
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
	       default => 10000000, # BIG_NUMBER, # -> keep all edges.
	      );

has n_extras =>  (
		  isa => 'Int', # the number randomly chosen nodes to add as neighbors in addition to near ones.
		  is => 'rw',
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
#  die "fasta file not defined.\n" if (! defined $args->{fasta_string} ) ;

  my $fasta_str = $args->{fasta_string} // die '$args->{fasta_string} not defined.', "\n";
  my $id_gobj = fasta_string_to_gobjs($fasta_str); # $args->{fasta});
  my @ids = sort { $a <=> $b } keys %$id_gobj;
  my $n_nodes = scalar @ids;
  
  if (defined $args->{graph_string}) { # construct from .graph file
    # construct from file with info on nearest neighbors to each node. typical line:
    # 290  2  ((()26,()24)218,(()52,()22)176)290    176  0.383   218  0.42   274  0.482   354  0.501   52  0.502  328  0.757
    my $n_near = BIG_NUMBER; # just keep everything in the .graph file
    my %idA__idB_distance = (); # hash of hash refs
    my $id_node = {};
    my @graph_string_lines = split("\n", $args->{graph_string});
    for my $line (@graph_string_lines) {
      next if($line =~ /^\s*#/);
      $line =~ s/^\s*(\S+)\s+(\S+)\s+(\S+)\s+//;
      my ($id1, $generation, $pedigree) = ($1, $2, $3);
      $line =~ s/\.{3}\s+(\S+)\s+(\S+)/.../;
      my ($furthest_id, $furthest_d) = ($1, $2);
      my $sequence = ($line =~ /\s*\.{3}\s+(\S+)\s*$/)? $1 : '---';
      $line =~ s/\s*\.{3}.*$//;
      my $neighborid_dist = {};
      my @neighbor_ids = ();
      while ($line) {
	$line =~ s/\s*(\S+)\s+(\S+)\s*//;
	my ($id2, $dist) = ($1, $2);
	$neighborid_dist->{$id2} = $dist;
	push @neighbor_ids, $id2;
      }
      $n_near = min($n_near, scalar @neighbor_ids);
      my $gobj = $id_gobj->{$id1};
      my $a_node = GenotypeGraphNode->new( {
					    id => $id1,
					    genotype => $gobj, # genotype object
					    neighbor_id_distance => $neighborid_dist,
					   } );
      $id_node->{$id1} = $a_node;

    }
    return {
	    nodes => $id_node,
	    distances => \%idA__idB_distance,
	   };

  } elsif (defined $args->{dmatrix_string}) { # construct from .dmatrix file
    
    # read in and store distances from
    my %idA__idB_distance = (); # hash of hash refs
    print "AA: [", $args->{dmatrix_string}, "]\n";
    my @lines = split("\n", $args->{dmatrix_string});
    my $first_line = shift @lines;
    my $multiplier = ($first_line =~ /#\s+(\d+)/)? $1 : undef;
    my @dm_ids = split(" ", shift @lines);
#    print "n lines in dmatrix string: ", scalar @md
    while(my($i, $line) = each @lines){
      my @distances = split(/\s+/, $line); # except that the first col is id
      my $id1 = shift @distances;
      while(my($j, $d) = each @distances){
	$d /= (defined $multiplier)? $multiplier : 1;
	my $id2 = $dm_ids[$j + $i + 1];
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
    print join(' ', @ids), "\n";
    print join(' ', keys %idA__idB_distance), "\n";
    # for each genotype, find the nearest neighbors and construct node.
    my $id_node = {};
    my $n_near = $args->{n_near} // BIG_NUMBER;
    $n_near = min($n_near, $n_nodes-1);
    my $n_extras = $args->{n_extras};

    while ( my ($i, $id1) = each @ids) {
      my $id2_dist = $idA__idB_distance{$id1};
      print "$id1  ", join(', ', %$id2_dist), "\n";
      my $count = 0;
      my @neighbor_ids;
      my $m = 'pq';
      if ($m eq 'qsel') { # using quickselect algorithm (a bit faster)
	@neighbor_ids = quickselect([keys %$id2_dist], $n_near, $id2_dist);
	@neighbor_ids = sort { $id2_dist->{$a} <=> $id2_dist->{$b} } @neighbor_ids;
      } elsif ($m eq 'sort') { # sort the whole set of nodes (a bit slower)
	@neighbor_ids = sort { $id2_dist->{$a} <=> $id2_dist->{$b} } keys %$id2_dist;
	@neighbor_ids = @neighbor_ids[0 .. $n_near-1] if($n_near < scalar keys %$id2_dist);
      } else {		   # a bit faster than qsel for small $n_near 
	my $pq = MyPriorityQueue->new($n_near, $id2_dist);
	@neighbor_ids = @{ $pq->{queue} };
      }
      my %neighborid_dist = map(($_ => $id2_dist->{$_}), @neighbor_ids); # hash w ids, distances for just nearest $n_near
      get_extra_ids($id2_dist, \@neighbor_ids, \%neighborid_dist, $n_extras);

      $id_node->{$id1} = GenotypeGraphNode->new( {
						  id => $id1,
						  genotype => $id_gobj->{$id1}, # genotype object
						  #		    neighbor_ids => \@neighbor_ids,
						  neighbor_id_distance => \%neighborid_dist,
						 } );

    }				# end node construction loop

    return {
	    nodes => $id_node,
	    n_near => $n_near,
	    n_extras => $n_extras,
	    distances => \%idA__idB_distance,
	   };

    
  } else {			# construct from .fasta string


    
    # calculate and store distances (N choose 2 of them)
    my %idA__idB_distance = (); # hash of hash refs
#    my @ids = sort { $a <=> $b } keys %$id_gobj;
#    my $n_nodes = scalar @ids;
    my ($edge_count, $d_sum) = (0, 0);
    while (my ($i1, $id1) = each @ids) {
      #   print STDERR "$i1 $id1 \n";
      my $g1 = $id_gobj->{$id1};
      for (my $i2 = $i1+1; $i2 < scalar @ids; $i2++) {
	my $id2 = $ids[$i2];
	#    print STDERR "    $i2 $id2 \n";
	my $g2 = $id_gobj->{$id2};
	my $d = $g1->distance($g2);
	$d += 1e-6*($id2 + $id1); ##########
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
      my $m = 'pq';
      if ($m eq 'qsel') { # using quickselect algorithm (a bit faster)
	@neighbor_ids = quickselect([keys %$id2_dist], $n_near, $id2_dist);
	@neighbor_ids = sort { $id2_dist->{$a} <=> $id2_dist->{$b} } @neighbor_ids;
      } elsif ($m eq 'sort') { # sort the whole set of nodes (a bit slower)
	@neighbor_ids = sort { $id2_dist->{$a} <=> $id2_dist->{$b} } keys %$id2_dist;
	@neighbor_ids = @neighbor_ids[0 .. $n_near-1] if($n_near < scalar keys %$id2_dist);
      } else {		   # a bit faster than qsel for small $n_near 
	my $pq = MyPriorityQueue->new($n_near, $id2_dist);
	@neighbor_ids = @{ $pq->{queue} };
      }
      my %neighborid_dist = map(($_ => $id2_dist->{$_}), @neighbor_ids); # hash w ids, distances for just nearest $n_near
      get_extra_ids($id2_dist, \@neighbor_ids, \%neighborid_dist, $n_extras);

      $id_node->{$id1} = GenotypeGraphNode->new( {
						  id => $id1,
						  genotype => $id_gobj->{$id1}, # genotype object
						  #		    neighbor_ids => \@neighbor_ids,
						  neighbor_id_distance => \%neighborid_dist,
						 } );

    }				# end node construction loop

    return {
	    nodes => $id_node,
	    n_near => $n_near,
	    n_extras => $n_extras,
	    distances => \%idA__idB_distance,
	   };
  }				# end of

  # otherwise do nothing.
};				# end of BUILDARGS

sub BUILD {
  my $self = shift;
  #my @nodes = values %{$self->nodes()};
  for my $node_obj (values %{$self->nodes()} ) {
    $node_obj->{graph} = $self; # can do it this way.
    $node_obj->symmetrize_neighbors();
  }

  # if (1) {
  #   for my $node_obj (values %{$self->nodes()}) {
  #     $node_obj->symmetrize_neighbors();
  #   }
  # }
}

sub get_node_by_id{
  my $self = shift;
  my $id = shift;
  return $self->nodes()->{$id};
}

sub search_for_best_match{
  my $self = shift;
  my $gobj = shift;		#  Genotype objects to match
  my $pq_size_limit = shift // 10;
  my $independent_starts = shift // 2;
  print "#  ", $gobj->id(), "\n";
  my %initid_bestmatchiddist = ();
  for (1..$independent_starts) {
    my $init_node_id = (keys %{$self->nodes()})[ int(rand(keys  %{$self->nodes()})) ];

    my $pq = MyPriorityQueue->new($pq_size_limit); # for storing the best-so-far nodes;
    my $pq_a = MyPriorityQueue->new($pq_size_limit); # for storing the best-so-far nodes;
    my $pq_h = MyPriorityQueue->new($pq_size_limit); # for storing the best-so-far nodes;
    my $id_status = {$init_node_id => 0}; # 0: unchecked, 1: checked, 2: checked and neighbors checked

  
    my $count_d_calcs = 0;
    my $count_rounds = 0;
    my $count_futile_rounds = 0; # count the number of rounds since a better candidate has been found - use for deciding when to stop.
    my $active_ids = {$init_node_id => 1}; # neighbors of these need to be checked.
    while (1) {
      my $neighbor_ids = {};
      my $inserted_ids = {};
      my $inserted_ids_a = {};
      my $inserted_ids_h = {};

      for my $an_id (keys %$active_ids) {
	# check these ids. i.e. get distances, insert in pq, and keep track of which have been checked,
	# and which are in the pq after they have all been added (and some possibly bumped).

	my ($d, $a_dist, $h_dist) = $gobj->distance($self->nodes()->{$an_id}->genotype());
            
	$count_d_calcs++;

	my ($inserted, $bumped_id) = $pq->size_limited_insert($an_id, $d);
	$inserted_ids->{$an_id} = 1 if($inserted);
	delete $inserted_ids->{$bumped_id} if(defined $bumped_id); # so if an id is inserted, then bumped, it will not be in this hash.

	my ($inserted_a, $bumped_id_a) = $pq_a->size_limited_insert($an_id, $a_dist);
	$inserted_ids_a->{$an_id} = 1 if($inserted_a);
	delete $inserted_ids_a->{$bumped_id_a} if(defined $bumped_id_a); # so if an id is inserted, then bumped, it will not be in this hash.

	my ($inserted_h, $bumped_id_h) = $pq_h->size_limited_insert($an_id, $h_dist);
	$inserted_ids_h->{$an_id} = 1 if($inserted_h);
	delete $inserted_ids_h->{$bumped_id_h} if(defined $bumped_id_h); # so if an id is inserted, then bumped, it will not be in this hash.

	$id_status->{$an_id} = 1; # this one has been checked!

	print $gobj->id(), "  $count_rounds  $count_d_calcs  ", join(' ', map($_ // '-', $pq->peek_best())), " $d  ",
	  join(' ', map($_ // '-', $pq_a->peek_best())), " $a_dist  ",
	  join(' ', map($_ // '-', $pq_h->peek_best())), " $h_dist\n";
      }

      $count_rounds++;
      $count_futile_rounds = (keys %$inserted_ids > 0)? 0 : $count_futile_rounds+1;
      last if($count_futile_rounds > 1 and $count_rounds > 1);

      for my $an_id (keys %$inserted_ids) { # get the neighbors of these (only those which have not been checked yet)
	#  for my $a_neighbor_id (@{$self->nodes()->{$an_id}->neighbor_ids()}) {
	for my $a_neighbor_id (keys %{$self->nodes()->{$an_id}->neighbor_id_distance()}) {
	  $id_status->{$a_neighbor_id} //= 0;
	  $neighbor_ids->{$a_neighbor_id} = 1 if($id_status->{$a_neighbor_id} == 0); # skip any neighbors which have been checked already.
	}
      }
      for my $an_id (keys %$active_ids) { # these have been check and the set of their neighbors will need to be checked has been defined.
	$id_status->{$an_id} = 2;
      }
      $active_ids = $neighbor_ids; # neighbors in this round become active nodes for next round.

      # print "rounds:  $count_rounds  $count_futile_rounds distance calcs: $count_d_calcs    ";
      # for(my $i=0; $i <= 2; $i++){
      #   my ($an_id, $dist) = $pq->i_th_best($i);
      #   last if(!defined $an_id);
      #   print "$an_id  $dist    ";
      # }print "\n";

    }				# end of a round
    # print "rounds:  $count_rounds  $count_futile_rounds distance calcs: $count_d_calcs    ";
    # for(my $i=0; $i <= 2; $i++){
    #   my ($an_id, $dist) = $pq->i_th_best($i);
    #   last if(!defined $an_id);
    #   print "$an_id  $dist    ";
    # }print "\n";
    my ($best_id, $best_dist) = $pq->best();
    print STDERR "$count_d_calcs  $best_id  $best_dist   ";
    my ($nextbest_id, $nextbest_dist) = $pq->best();
    print STDERR "$nextbest_id  $nextbest_dist   ";
    $initid_bestmatchiddist{$gobj->id()} = $best_id . "_" . sprintf("%f6.4", $best_dist);
  }
  print STDERR "\n";
  return \%initid_bestmatchiddist
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

sub quickselect{                # get the nearest $k ids.
  my $id_list = shift;
  my $k = shift;	      # find this many.
  my $id_distance = shift;    # hash ref of all N-1 id:distance pairs.

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


sub fasta_string_to_gobjs{
  my $fasta_string = shift;
  my @fasta_lines = split("\n", $fasta_string );
  my %id_gobj = ();
  while (@fasta_lines) {
    my $id_line = shift @fasta_lines;
    my $sequence = shift @fasta_lines;
    my $gobj = Genotype->new( $id_line . "\n" . $sequence );
    $id_gobj{$gobj->id()} = $gobj;
  }
  return \%id_gobj;
}

sub store_idAidBdistance{
  my $idA__idB_distance = shift; # hashref of hashrefs
  my $idA = shift;
  my $idB = shift;
  my $d = shift;

  	if (!exists $idA__idB_distance->{$idA}) {
	  $idA__idB_distance->{$idA} = {$idB => $d};
	} else {
	  $idA__idB_distance->{$idA}->{$idB} = $d;
	}
	if (!exists $idA__idB_distance->{$idB}) {
	  $idA__idB_distance->{$idB} = {$idA => $d};
	} else {
	  $idA__idB_distance->{$idB}->{$idA} = $d;
	}

  return $idA__idB_distance;
}

###########################################
__PACKAGE__->meta->make_immutable;

1;


# __DATA__

# __C__

# void qselect(

