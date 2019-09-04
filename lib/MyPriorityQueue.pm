package MyPriorityQueue;

our $VERSION = '0.01';

use strict;
use warnings;

sub new {
  my $classname = shift;
 # print $classname, "\n";
  my $size_limit = shift;
  return bless {
		queue => [],
		prios => {},	# by payload
		size_limit => undef,
	       }, $classname;
  
}

sub best {	  # shift the 'best' (highest priority) element queue.
  my ($self) = @_;
  if (@{$self->{queue}} == 0) {
    return undef;
  }
  delete($self->{prios}->{$self->{queue}->[0]});
  return shift(@{$self->{queue}});
}

sub worst { # pop the 'worst' (lowest priority) element off the queue.
  my ($self) = @_;
  if (@{$self->{queue}} == 0) {
    return (undef, undef);
  }
  my $worst_priority = delete($self->{prios}->{$self->{queue}->[-1]});
  return (pop(@{$self->{queue}}), $worst_priority);
}

sub set_size_limit{
  my $self = shift;
  $self->{size_limit} = shift;
}

sub get_size_limit{
  my $self = shift;
  return $self->{size_limit};
}

sub get_size {
  my ($self) = @_;
  return scalar @{$self->{queue}};
}

sub quickselect{ # 
  my $self = shift;
  my $id_list = shift;
  my $k = shift;
 my $id_distance = shift; # hash ref of all N-1 id:distance pairs.

 #   print "id_list: ", join(' ', @$id_list), "  k: $k \n";
my $rand_index =  
# int ( 0.5*@$id_list ) ;
  int rand @$id_list;
#  int (0.5* (rand @{ $id_list } + rand @{$id_list}));
#  int (0.333* (rand @{ $id_list } + rand @{$id_list} + rand @{$id_list})) - 1;
#  int (0.4*@$id_list + 0.2*rand @$id_list) - 1;

  my $pivot_id = $id_list->[$rand_index ]; 
  my $pivot = $id_distance->{$pivot_id};

  my @lefts = (); my @rights = (); my @equals = ();
  if (0) { # use built-in grep, but run through array 2 (or sometimes 3) times; a bit slower.
    @lefts  = grep { $id_distance->{$_} < $pivot } @$id_list;
    @rights = grep { $id_distance->{$_} > $pivot } @$id_list;
    # my @equals = grep { $_ == $pivot } @$list;
    # my @equals = ();
    if (@lefts + @rights + 1 == scalar @$id_list) {
      push @equals, $pivot_id;
    } else {
      @equals = grep { $id_distance->{$_} == $pivot } @$id_list;
    }
  } else { # a bit faster
    for (@$id_list) { # store in separate arrays the ids with distance <, ==, and > the pivot 
      if ($id_distance->{$_} < $pivot) {
	push @lefts, $_;
      } elsif ($id_distance->{$_} > $pivot) {
	push @rights, $_;
      } else {
	push @equals, $_;
      }
    }
 }

  if ($k < @lefts) {  # kth will be in @lefts, but lefts has too many.
    return $self->quickselect(\@lefts, $k, $id_distance);
  } elsif ($k > @lefts + @equals) { # kth will be in @rights
    return ((@lefts, @equals), $self->quickselect(\@rights, $k - @lefts - @equals, $id_distance));
  } elsif ($k == @lefts) {	# done 
    return @lefts;
  } elsif ($k <= @lefts + @equals) { # just @lefts plus 1 or more from @equals
    push @lefts, @equals[0..$k-@lefts-1];
    return @lefts;
  } else {
    die "???\n";
  }
}

sub size_limited_array_insert{
   my $self = shift;
   my $id_distance = shift;
   my $k = shift;
   my @ids = keys %$id_distance;
   my $top_k_ids = $self->quick_select(\@ids, $k, $id_distance);
   for my $topid (@$top_k_ids) {
      $self->insert($topid, $id_distance->{$topid}); # 
   }
}

sub size_limited_insert{
  my ($self, $payload, $priority) = @_;
  my $worst_priority = $self->{prios}->{$self->{queue}->[-1]};
  #  assert ($self->size() <= $self->{size_limit});
  if(defined $self->{size_limit}){
  my $at_limit = scalar @{$self->{queue}} == $self->{size_limit};
  if ($priority < $worst_priority) {
    $self->insert($payload, $priority);
    $self->worst() if($at_limit); # if exceeds size limit pop worst.
  } elsif (!$at_limit) {
    $self->insert($payload, $priority);
  }
}else{ # insert without regard to size of queue
  $self->insert($payload, $priority);
}
}
sub unchecked_insert {
  my ($self, $payload, $priority, $lower, $upper) = @_;
  $lower = 0                             unless defined($lower);
  $upper = scalar(@{$self->{queue}}) - 1 unless defined($upper);

  # first of all, map the payload to the desired priority
  # run an update if the element already exists
  $self->{prios}->{$payload} = $priority;

  # And register the payload in the queue. There are a lot of special
  # cases that can be exploited to save us from doing the relatively
  # expensive binary search.

  # Special case: No items in the queue.  The queue IS the item.
  if (@{$self->{queue}} == 0) {
    push(@{$self->{queue}}, $payload);
    return;
  }

  # Special case: The new item belongs at the end of the queue.
  if ($priority >= $self->{prios}->{$self->{queue}->[-1]}) {
    push(@{$self->{queue}}, $payload);
    return;
  }

  # Special case: The new item belongs at the head of the queue.
  if ($priority < $self->{prios}->{$self->{queue}->[0]}) {
    unshift(@{$self->{queue}}, $payload);
    return;
  }

  # Special case: There are only two items in the queue.  This item
  # naturally belongs between them (as we have excluded the other
  # possible positions before)
  if (@{$self->{queue}} == 2) {
    splice(@{$self->{queue}}, 1, 0, $payload);
    return;
  }

  # And finally we have a nontrivial queue.  Insert the item using a
  # binary seek.
  # Do this until the upper and lower bounds crossed... in which case we
  # will insert at the lower point
  my $midpoint;
  while ($upper >= $lower) {
    $midpoint = ($upper + $lower) >> 1;

    # We're looking for a priority lower than the one at the midpoint.
    # Set the new upper point to just before the midpoint.
    if ($priority < $self->{prios}->{$self->{queue}->[$midpoint]}) {
      $upper = $midpoint - 1;
      next;
    }

    # We're looking for a priority greater or equal to the one at the
    # midpoint.  The new lower bound is just after the midpoint.
    $lower = $midpoint + 1;
  }

  splice(@{$self->{queue}}, $lower, 0, $payload);
}

sub _find_payload_pos {
  my ($self, $payload) = @_;
  my $priority = $self->{prios}->{$payload};
  if (!defined($priority)) {
    return undef;
  }

  # Find the item with binary search.
  # Do this until the bounds are crossed, in which case the lower point
  # is aimed at an element with a higher priority than the target
  my $lower = 0;
  my $upper = @{$self->{queue}} - 1;
  my $midpoint;
  while ($upper >= $lower) {
    $midpoint = ($upper + $lower) >> 1;

    # We're looking for a priority lower than the one at the midpoint.
    # Set the new upper point to just before the midpoint.
    if ($priority < $self->{prios}->{$self->{queue}->[$midpoint]}) {
      $upper = $midpoint - 1;
      next;
    }

    # We're looking for a priority greater or equal to the one at the
    # midpoint.  The new lower bound is just after the midpoint.
    $lower = $midpoint + 1;
  }

  # The lower index is now pointing to an element with a priority higher
  # than our target.  Scan backwards until we find the target.
  while ($lower-- >= 0) {
    return $lower if ($self->{queue}->[$lower] eq $payload);
  }
}

sub delete {
  my ($self, $payload) = @_;
  my $pos = $self->_find_payload_pos($payload);
  if (!defined($pos)) {
    return undef;
  }

  delete($self->{prios}->{$payload});
  splice(@{$self->{queue}}, $pos, 1);

  return $pos;
}

sub unchecked_update {
  my ($self, $payload, $new_prio) = @_;
  my $old_prio = $self->{prios}->{$payload};

  # delete the old item
  my $old_pos = $self->delete($payload);

  # reinsert the item, limiting the range for the binary search (if needed)
  # a bit by checking how the priority changed.
  my ($upper, $lower);
  if ($new_prio - $old_prio > 0) {
    $upper = @{$self->{queue}};
    $lower = $old_pos;
  } else {
    $upper = $old_pos;
    $lower = 0;
  }
  $self->unchecked_insert($payload, $new_prio, $lower, $upper);
}

sub update {
  my ($self, $payload, $prio) = @_;
  if (!defined($self->{prios}->{$payload})) {
    goto &unchecked_insert;
  } else {
    goto &unchecked_update;
  }
}
*insert = \&update;

1;

__END__

=head1 NAME

List::PriorityQueue - high performance priority list (pure perl)

=head1 SYNOPSIS

 my $prio = new List::PriorityQueue;
 $prio->insert("foo", 2);
 $prio->insert("bar", 1);
 $prio->insert("baz", 3);
 my $next = $prio->pop(); # "bar"
 # I decided that "foo" isn't as important anymore
 $prio->update("foo", 99);

=head1 DESCRIPTION

This module implements a high-performance priority list. It's written in pure
Perl.

Available functions are:

=head2 B<new>()

Obvious.

=head2 B<insert>(I<$payload>, I<$priority>)

=head2 B<update>(I<$payload>, I<$new_priority>)

Adds the specified payload (anything fitting into a scalar) to the priority
queue, using the specified priority. Smaller means more important.

If the item already exists in the queue, it is assigned the new priority.
It's optimized to perform better than a delete followed by insert.

These names are actually the same function. The alternative name is provided
so you can make clear which operation you intended to be executed.

=head2 B<pop>()

Removes the most important item (numerically lowest priority) from the queue
and returns it. If no element is there, returns I<undef>.

=head2 B<delete>(I<$payload>)

Deletes an item known by the specified payload from the queue.

=head2 B<unchecked_insert>(I<$payload>, I<$priority>)

=head2 B<unchecked_update>(I<$payload>, I<$new_priority>)

These functions are provided as an alternative to the safe (or "checked")
functions described above. By bypassing some checks, they gain you a speed
advantage, but if you don't know what you're doing, using these might
corrupt/confuse the queue.

You can use these if you I<definitely> know that the element doesn't exist
yet for B<unchecked_insert>, or that the element definitely already exists
for B<unchecked_update>.

=head1 DIFFERENCES TO POE::QUEUE::ARRAY

There are some things I disliked about POE::Queue::Array, which ultimately
led to the creation of this derivative.

First, it stores data in a package global variable. The author brings up a
valid argument why this is not bad in this case. However I still was somehow
not happy with the fact it used a global variable. For example, serializing
a queue would not work as the actual queue reference only stores numerical
IDs into the package variable containing all the data, but that one wouldn't
be saved.

Second, for some operations to be carried out efficiently, you have to carry
the internal IDs around in your program, else you have to do a full search
for your element everytime you want to delete or update it. While carrying it
around is relatively simple, there is no reason why the class itself shouldn't
manage this and relieve the programmer from this work.

A benchmark with POE::Queue::Array and the described payload-to-ID mapping and
this module revealed that they are equally fast.

=head1 BUGS

Maybe.

=head1 SEE ALSO

L<POE::Queue::Array>

=head1 AUTHORS & COPYRIGHTS

Made 2009 by Lars Stoltenow.
List::PriorityQueue is free software; you may redistribute it and/or modify it
under the same terms as Perl itself.

POE::Queue::Array is Copyright 1998-2007 Rocco Caputo. All rights reserved.
Same license.

=cut
