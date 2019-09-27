package GenotypeGraph;
use strict;
use warnings;
#use Moose;
use Mouse;
use namespace::autoclean;
use Carp;
use List::Util qw ( min max sum );
# use Graph;
use MyPriorityQueue;

use constant BIG_NUMBER => 1_000_000_000;

# class to represent a graph, i.e. a set of node, which may be
# connected in pairs by edges. Here the nodes know which other nodes
# they are connected to. 

has nodes => (
              isa => 'HashRef', # keys are ids, values GenotypeGraphNode objects.
              is => 'ro',
              default => sub { {} },
             );

has n_keep => (
	       isa => 'Int', # the number of outward edges to keep in the graph for each node.
	       is => 'ro',
	       default => BIG_NUMBER, # -> keep all edges.
	      );

has n_near => (
	       isa => 'Int', # the number of outward edges to keep in each node object of the graph.
	       is => 'ro',
	       default => BIG_NUMBER, # -> keep all edges.
	      );

has n_extras =>  (
		  isa => 'Int', # the number randomly chosen nodes to add as neighbors in addition to near ones.
		  is => 'rw',
		  default => 0, # just keep the closest ones.
		 );

has distances => ( #
                  isa => 'HashRef',
                  is => 'rw',
                  default => sub { {} },
                 );

around BUILDARGS => sub {
  my $orig = shift;
  my $class = shift;

  my $args = shift; # hash ref, e.g. { fasta => <a fasta string for several sequences> , n_near => 5 }

  my $fasta_str = $args->{fasta_string} // die '$args->{fasta_string} not defined.', "\n";
  my $id_gobj = fasta_string_to_gobjs($fasta_str); # $args->{fasta});
  my @ids = sort { $a <=> $b } keys %$id_gobj;
  my $n_keep = $args->{n_keep} // BIG_NUMBER; 
  my $n_near = $args->{n_near} // BIG_NUMBER;
  my $n_extras = $args->{n_extras} // 0;
  my $n_nodes = scalar @ids;

  if (defined $args->{graph_string}) { # construct from .graph file    ### FROM GRAPH ###
    # #### construct from .graph file, with info on nearest neighbors to each node. typical line:
    # #### 290  2  ((()26,()24)218,(()52,()22)176)290    176  0.383   218  0.42   274  0.482   354  0.501   52  0.502  328  0.757

    my $idA__idB_distance = {}; # hash of hash refs
    my $id_node = {};
    my @graph_string_lines = split("\n", $args->{graph_string});
    my $first_line = shift @graph_string_lines;
    my ($n_near, $n_extras) = ($first_line =~ /n_near:\s+(\d+)\s+n_extras:\s+(\d+)/)? ($1, $2) : (BIG_NUMBER, 0);
    for my $line (@graph_string_lines) {
      next if($line =~ /^\s*#/);
      $line =~ s/^\s*(\S+)\s+(\S+)\s+(\S+)\s+//;
      my ($id1, $generation, $pedigree) = ($1, $2, $3);
      $line =~ s/\.{3}\s+(\S+)\s+(\S+)/.../;
      my ($furthest_id, $furthest_d) = ($1, $2);
      my $sequence = ($line =~ /\s*\.{3}\s+(\S+)\s*$/)? $1 : '---';
      $line =~ s/\s*\.{3}.*$//;
      my $neighborid_dist = {};
      while ($line) {
	$line =~ s/\s*(\S+)\s+(\S+)\s*//;
	my ($id2, $dist) = ($1, $2);
	$neighborid_dist->{$id2} = $dist;
	store_idAidBdistance($idA__idB_distance, $id1, $id2, $dist);
      }
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
	    distances => $idA__idB_distance,
	    n_keep => $n_keep,
	    n_near => $n_near,
	    n_extras => $n_extras,
	   };

  } elsif (defined $args->{dmatrix_string}) { # construct from .dmatrix file  ### FROM DISTANCE MATRIX ###
    # #### construct graph from .dmatrix file
    # #### which has all  N choose 2 distances

    # read in and store distances from .dmatrix file.
    my $idA__idB_distance = {}; # hash ref of hash refs
    my @lines = split("\n", $args->{dmatrix_string});
    my $first_line = shift @lines;
    my $multiplier = ($first_line =~ /#\s+(\d+)/)? $1 : undef;
    my @dm_ids = map($_*1, split(" ", shift @lines));
    while (my($i, $line) = each @lines) {
      $line =~ s/^\s+//;
      my @distances = split(/\s+/, $line); # except that the first col is id
      my $id1 = shift @distances;
      while (my($j, $d) = each @distances) {
	if ($d =~ /\-+/){
	  # 
	}else {
	  $d /= (defined $multiplier)? $multiplier : 1;
	  my $id2 = $dm_ids[$j + $i + 1];
	  store_idAidBdistance($idA__idB_distance, $id1, $id2, $d);
	}
      }
    }

    # for each genotype, find the nearest neighbors and construct node.
    my $id_node = construct_nodes($id_gobj, $idA__idB_distance, $n_near, $n_extras);

    return {
	    nodes => $id_node,
	    n_near => $n_near,
	    n_extras => $n_extras,
	    distances => $idA__idB_distance,
	   };

    
  } else {		       # construct from .fasta string    ### FROM FASTA ###
    # #### construct graph from .fasta file
    # #### with N equal-length sequences. Calculate all N choose 2 distances
     
    # calculate and store distances (N choose 2 of them)
    if (0) { # store all
      my $idA__idB_distance = {}; # hash of hash refs
      while (my ($i1, $id1) = each @ids) {
	my $g1 = $id_gobj->{$id1};
	for (my $i2 = $i1+1; $i2 < scalar @ids; $i2++) {
	  my $id2 = $ids[$i2];
	  my $g2 = $id_gobj->{$id2};
	  my ($d, $a, $h, $o) = $g1->distance($g2);
	  store_idAidBdistance($idA__idB_distance, $id1, $id2, $d);
	}
      }

      # for each genotype, find the nearest neighbors and construct node.
      my $id_node = construct_nodes($id_gobj, $idA__idB_distance, $n_near, $n_extras);

      return {
	      nodes => $id_node,
	      n_near => $n_near,
	      n_extras => $n_extras,
	      distances => $idA__idB_distance,
	     };
    } else { # store $n_keep best in a PQ
      # my $idA__idB_distance = {}; # hash ref of pqs
      my %idA_idBdPQ = map(($_ => MyPriorityQueue->new($n_keep)) , @ids);
      while (my ($i1, $id1) = each @ids) {
	my $g1 = $id_gobj->{$id1};
	for (my $i2 = $i1+1; $i2 < scalar @ids; $i2++) {
	  my $id2 = $ids[$i2];
	  my $g2 = $id_gobj->{$id2};
	  my ($d, $a, $h, $o) = $g1->distance($g2);
#	  print STDERR "$id1  $id2  $d  \n";
	  $idA_idBdPQ{$id1}->size_limited_insert($id2, $d); # store in pq, keeping only top $n_keep distances
	  $idA_idBdPQ{$id2}->size_limited_insert($id1, $d); # store in pq, keeping only top $n_keep distances
	  #	  print STDERR "AAA n best: ", $idA_idBdPQ{$id1}->peek_n_best(4), "\n";
#	  print STDERR "after sli\n";
	}
      }

      # for each genotype, find the nearest neighbors and construct node.
      my $id_node = construct_nodes($id_gobj, \%idA_idBdPQ, $n_near, $n_extras);

      return {
	      nodes => $id_node,
	      n_keep => $n_keep,
	      n_near => $n_near,
	      n_extras => $n_extras,
	      distances => \%idA_idBdPQ,
	     };
    }
  }
  # otherwise do nothing, leave arguments unaltered.
};				# end of BUILDARGS

sub BUILD {
  my $self = shift;
  #my @nodes = values %{$self->nodes()};
  for my $node_obj (values %{$self->nodes()} ) {
    $node_obj->{graph} = $self; # can do it this way.
    $node_obj->symmetrize_neighbors();
  }
}

sub get_node_by_id{
  my $self = shift;
  my $id = shift;
  return $self->nodes()->{$id};
}

sub exhaustive_search{
  my $self = shift;
  my $gobj = shift;
  my $pq_size_limit = shift // 5;

  my $id_node = $self->nodes();

  my $pq =  MyPriorityQueue->new($pq_size_limit);
  while (my ($id, $node) = each %$id_node) {
    my ($d, $a, $h, $o) = $gobj->distance($node->genotype());
    $pq->size_limited_insert($id, $d);
  }
  my $search_out_string = $gobj->id() . "    ";
  my @bests = $pq->peek_n_best(5);
  for (@bests) {
    $search_out_string .= sprintf("%5s %7.5f   ", @$_);
  }
  #  $search_out_string .= "\n";
  return $search_out_string;
}

sub search_for_best_match{
  my $self = shift;
  my $gobj = shift;		#  Genotype objects to match
  my $independent_starts = shift // 2;
  my $pq_size_limit = shift // 10;
  my $n_futile_rounds = shift // 1;
  my $Ninit =  1;
  
  my $search_out_string = $gobj->id() . "    ";
  #  my %initid_bestmatchiddist = ();
  for (1..$independent_starts) {
    my $init_node_id = (keys %{$self->nodes()})[ int(rand(keys  %{$self->nodes()})) ];
     my @init_node_ids = ();
    for (1..$Ninit) {
      push @init_node_ids, (keys %{$self->nodes()})[ int(rand(keys  %{$self->nodes()})) ];
    }
    my $pq = MyPriorityQueue->new($pq_size_limit); # for storing the best-so-far nodes;
#    my $pq_a = MyPriorityQueue->new($pq_size_limit); # for storing the best-so-far nodes;
#    my $pq_h = MyPriorityQueue->new($pq_size_limit); # for storing the best-so-far nodes;
    my %id_status = map(($_ => 0), @init_node_ids); # 0: unchecked, 1: checked, 2: checked and neighbors checked

    my $count_d_calcs = 0;
    my $count_rounds = 0;
    my $count_futile_rounds = 0; # count the number of rounds since a better candidate has been found - use for deciding when to stop.
    my $active_ids = {map(($_ => 1), @init_node_ids)}; # neighbors of these need to be checked.
    while (1) {
      my $neighbor_ids = {};
      my $inserted_ids = {};
      my $inserted_ids_a = {};
      my $inserted_ids_h = {};

      for my $an_id (keys %$active_ids) {
	# check these ids. i.e. get distances, insert in pq, and keep track of which have been checked,
	# and which are in the pq after they have all been added (and some possibly bumped).

	my ($d, $a_dist, $h_dist, $o_dist) = $gobj->distance($self->nodes()->{$an_id}->genotype());
	$count_d_calcs++;

	my ($inserted, $bumped_id) = $pq->size_limited_insert($an_id, $d);
	$inserted_ids->{$an_id} = 1 if($inserted);
	delete $inserted_ids->{$bumped_id} if(defined $bumped_id); # so if an id is inserted, then bumped, it will not be in this hash.

	# if (0) { # pq's for the other kinds of distances - skip for now.
	#   my ($inserted_a, $bumped_id_a) = $pq_a->size_limited_insert($an_id, $a_dist);
	#   $inserted_ids_a->{$an_id} = 1 if($inserted_a);
	#   delete $inserted_ids_a->{$bumped_id_a} if(defined $bumped_id_a); # so if an id is inserted, then bumped, it will not be in this hash.

	#   my ($inserted_h, $bumped_id_h) = $pq_h->size_limited_insert($an_id, $h_dist);
	#   $inserted_ids_h->{$an_id} = 1 if($inserted_h);
	#   delete $inserted_ids_h->{$bumped_id_h} if(defined $bumped_id_h); # so if an id is inserted, then bumped, it will not be in this hash.
	# }

	
	$id_status{$an_id} = 1; # this one has been checked!

#		print $gobj->id(), "  $count_rounds  $count_d_calcs  ", join(' ', map($_ // '-', $pq->peek_best())), " $d  \n"; #,
	#	  join(' ', map($_ // '-', $pq_a->peek_best())), " $a_dist  ",
	#	  join(' ', map($_ // '-', $pq_h->peek_best())), " $h_dist\n";
      }

      $count_rounds++;
      $count_futile_rounds = (keys %$inserted_ids > 0)? 0 : $count_futile_rounds+1;
      last if($count_futile_rounds > $n_futile_rounds and $count_rounds > 1);

      for my $an_id (keys %$inserted_ids) { # get the neighbors of these (only those which have not been checked yet)
	#  for my $a_neighbor_id (@{$self->nodes()->{$an_id}->neighbor_ids()}) {

	# just get top few here
	for my $a_neighbor_id (keys %{$self->nodes()->{$an_id}->neighbor_id_distance()}) {
	  $id_status{$a_neighbor_id} //= 0;
	  $neighbor_ids->{$a_neighbor_id} = 1 if($id_status{$a_neighbor_id} == 0); # skip any neighbors which have been checked already.
	}
      }
      for my $an_id (keys %$active_ids) { # these have been checked and the set of their neighbors which will need to be checked has been defined.
	$id_status{$an_id} = 2;
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
    # output the best and next-best matches: 
    #  my ($best_id, $best_dist) = $pq->best();
    $search_out_string .= sprintf("%4d    ", $count_d_calcs);
 #   my ($nextbest_id, $nextbest_dist) = $pq->best();
    #  $search_out_string .= sprintf("%5s %7.5f      ",  $nextbest_id, $nextbest_dist);
    #   $initid_bestmatchiddist{$gobj->id()} = $best_id . "_" . sprintf("%f6.4", $best_dist);
    my @bests = $pq->peek_n_best(4);
    for (@bests) {
      $search_out_string .= sprintf("%5s %8.5f   ", @$_);
    }
  }
  #  $search_out_string .= "\n";
  #  print $search_out_string;
  return $search_out_string;	# \%initid_bestmatchiddist
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
   #  my $idA__idB_distance = $self->distances();
   my @ids = sort {$a <=> $b} keys %{$self->distances()}; # idA__idB_distance;
   my $n_nodes = scalar @ids;
   my $d_matrix_string = '# ' . $multiplier ."\n";
   $d_matrix_string .= join(" ", @ids) . "\n";
   while ( my ($i, $id1) = each @ids) {
    #  my $id2_dist = $self->distances()->{$id1}; # a MyPriorityQueue  # $idA__idB_distance->{$id1};
      my @n_bests = $self->distances()->{$id1}->peek_n_best($self->n_keep());
      my %id2_dist = map(@$_ , @n_bests);
#      print STDERR "AAAA: ", join('   ',  map(@$_ , @n_bests)), "\n";

  #    my @id2s = sort { $a <=> $b } keys %id2_dist;
 #     print STDERR "BBBB: ", join(' ', keys %id2_dist), "     ", join(' ', @id2s), "\n";
 #     print STDERR "ASDFASDFASDF: ", join(" ", @id2s), "\n";
      $d_matrix_string .= sprintf("%2d  ", $id1); # , join(" ", map (int($multiplier*$id2_dist{$_} + 0.5), @id2s) ) );
      for (my $j=$i+1; $j < scalar @ids; $j++) {
	my $id2 = $ids[$j];
	my $d = $id2_dist{$id2} // $self->distances()->{$id2}->priority($id1);
         $d_matrix_string .= (defined $d)? sprintf("%d ", int($multiplier*$d + 0.5)) : '- ';
      }
      $d_matrix_string .= "\n";
   }
   return $d_matrix_string;
}

####   ordinary subroutines  ###

sub get_extra_ids{
  my $id2_dist = shift; # may be PQ
  my $nearid_dist = shift; # the ids already found to be the nearest (not symmetrized yet), and dists (hash ref)
    my $n_extras = shift;
  my $n_try = shift // 5;
  
  if (ref $id2_dist eq 'HASH') {
    #  my $near_ids = shift;
  
    for (1..$n_extras) {
      for (1..$n_try) {
	my $extra_id = (keys %$id2_dist)[int(rand keys %$id2_dist)];
	#	(keys %$id2_dist)[0 - $_];
	if (! exists $nearid_dist->{$extra_id}) {
	  #	push @$near_ids, $extra_id;
	  $nearid_dist->{$extra_id} = $id2_dist->{$extra_id};
	  last;
	}
      }
    }
  } elsif (ref $id2_dist eq 'MyPriorityQueue') {
 for (1..$n_extras) {
   for (1..$n_try) {
     my ($extra_id, $extra_d) = $id2_dist->i_th_best(int(rand $id2_dist->get_size()));
 #    print STDERR "xid xd: $extra_id  $extra_d \n";
	#	(keys %$id2_dist)[0 - $_];
	if (! exists $nearid_dist->{$extra_id}) {
	  #	push @$near_ids, $extra_id;
	  $nearid_dist->{$extra_id} = $id2_dist->priority($extra_id); # add to the set of ids to be considered neighbors 
	  last;
	}
      }
    }
}
 #print STDERR "XXX: ", join("  ", map($_ . " " . $nearid_dist->{$_}, keys %$nearid_dist)), "\n";
}

sub quickselect{		# get the nearest $k ids.
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

sub construct_nodes{ # construct the nodes, each with the appropriate neighbors
  my $id_gobj = shift;
  my $idAidBdist = shift;
  my $n_near = min(shift, (scalar keys %$id_gobj) - 1);
  my $n_extras = shift;

  my $id_node = {};
  for my $id1 (keys %$id_gobj) {
    my $id2_dist = $idAidBdist->{$id1};
#    print STDERR "ref id2_dist:  [", ref $id2_dist, "]\n";
    if (ref $id2_dist eq 'HASH') {
      my $count = 0;
      my @neighbor_ids;
      my $m = 'pq'; # controls which method is used to find the neighboring genotypes
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
      get_extra_ids($id2_dist, \%neighborid_dist, $n_extras);

      $id_node->{$id1} = GenotypeGraphNode->new( {
						  id => $id1,
						  genotype => $id_gobj->{$id1}, # genotype object
						  #		    neighbor_ids => \@neighbor_ids,
						  neighbor_id_distance => \%neighborid_dist,
						 } );
    } elsif (ref $id2_dist eq 'MyPriorityQueue') {
      my @n_best = map( ($_->[0], $_->[1]), $id2_dist->peek_n_best($n_near));
  #     print STDERR "#### id1:  $id1 \n";
  #    print STDERR "n best:  ", join('  ', @n_best), "\n";
      my %neighborid_dist = @n_best; # id:dist hash
      #map(($_ => $id2_dist->{$_}), @neighbor_ids); # hash w ids, distances for just nearest $n_near

      get_extra_ids($id2_dist, \%neighborid_dist, $n_extras);
   #   print STDERR "\n";
      $id_node->{$id1} = GenotypeGraphNode->new( {
						  id => $id1,
						  genotype => $id_gobj->{$id1}, # genotype object
						  #		    neighbor_ids => \@neighbor_ids,
						  neighbor_id_distance => \%neighborid_dist,
						 } );
    }

  }				# end node construction loop

  return $id_node;
}

###########################################
__PACKAGE__->meta->make_immutable;

1;
