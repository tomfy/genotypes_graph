package ChunkSet;
use Moose;
#use Mouse;
use namespace::autoclean;
use Carp::Assert;
use List::Util qw ( min max sum );
#use Inline 'C';

has n_chunks => (
		 isa => 'Maybe[Int]',
		 is => 'ro',
		 required => 1,
		);

has chunk_size => (
		   isa => 'Int',
		   is => 'ro',
		   required => 1,
		  );

has sequence_length => (
			isa => 'Int',
			is => 'ro',
			required => 1,
		       );

has chunk_specifiers => ( # array ref of strings, e.g. ['1_22_108_23', ...
			 isa => 'ArrayRef',
			 is => 'ro',
			 required => 0,
			);

has chunk_spec_arrays => ( # array ref of array refs, e.g. [ [1,22,108,23], ...
                          isa => 'ArrayRef',
                          is => 'ro',
                          required => 0,
                         );

has seqid_seq => (
		  isa => 'HashRef',
		  is => 'ro',
		  required => 1,
		 );

has chunkspec__seq_ids => (
			   isa => 'HashRef',
			   is => 'ro',
			   default => sub { {} },
			  );


sub BUILD{
  my $self = shift;
  my ($chunk_size, $seq_length) = ($self->chunk_size(), $self->sequence_length());
  my $n_chunks = $self->n_chunks() // int( $seq_length / $chunk_size );
#  print "AAA: $chunk_size $n_chunks $seq_length \n";

  my $n_snps_to_use = min($seq_length, $n_chunks*$chunk_size);
  my @chunk_specs = @{ chunk_index_strings([0 .. $n_snps_to_use], $chunk_size ) };
  while (@chunk_specs < $n_chunks) {
    $n_snps_to_use = min($seq_length, $n_chunks*$chunk_size);
    push @chunk_specs, @{ chunk_index_strings( randomize_array([0 .. $n_snps_to_use - 1]), $chunk_size ) };
  }
 # print join('  ', @chunk_specs), "\n";

  $self->{chunk_specifiers} = \@chunk_specs;
  for (@chunk_specs) {
    $self->chunkspec__seq_ids()->{$_} = {};
  }
  my @chunk_spec_arrays = map( [split('_', $_)], @chunk_specs ); # array of array refs, each of which holds indices of one chunk
  $self->{chunk_spec_arrays} = \@chunk_spec_arrays;

  while ( my($seqid, $seq) = each %{$self->seqid_seq()} ) { # loop over sequences
    my @seq_chars = split('', $seq);
    while ( my ($ich, $ch_indices) = each @chunk_spec_arrays) { # index, and array ref.
      my $ch_spec = $self->{chunk_specifiers}->[$ich]; # as a string
      my $chunk_seq = join('', @seq_chars[@$ch_indices]);
      push @{ $self->{chunkspec__seq_ids}->{$ch_spec}->{$chunk_seq} //= []  }, $seqid;
    }
  }
}


# sub BUILD{
#   my $self = shift;

#   for (@{$self->chunk_specifiers()}) {
#     $self->chunkspec__seq_ids()->{$_} = {};
#   }

#   my @chunk_spec_arrays = map( [split('_', $_)], @{$self->chunk_specifiers()} );
#   $self->{chunk_spec_arrays} = [map( [split('_', $_)], @{$self->chunk_specifiers()} )]; 
#   #  $self->chunk_size( scalar @{ $chunk_spec_arrays[0] } );
#   while ( my($seqid, $seq) = each %{$self->seqid_seq()} ) { # loop over sequences
#     my @seq_chars = split('', $seq);
#     while ( my ($ich, $ch_indices) = each @chunk_spec_arrays) { # index, and array ref.
#       my $ch_spec = $self->{chunk_specifiers}->[$ich]; # as a string
#       my $chunk_seq = join('', @seq_chars[@$ch_indices]);
#       push @{ $self->{chunkspec__seq_ids}->{$ch_spec}->{$chunk_seq} //= []  }, $seqid;
#     }
#   }
# }

sub get_chunk_match_counts{
  my $self = shift;
  my $sequence = shift;		# the sequence to match
  my $id_matchcount = shift // {};

  my $matches_count = 0;
  my @seq_chars = split('', $sequence);
  while ( my ($ich, $ch_indices) = each @{$self->{chunk_spec_arrays}}) { # index, and array ref.
    my $ch_spec =  $self->{chunk_specifiers}->[$ich]; # as a string
    my $chunk_seq = join('', @seq_chars[@$ch_indices]);
    my $id_matches = $self->{chunkspec__seq_ids}->{$ch_spec}->{$chunk_seq} // [];
    for ( @{ $id_matches } ) {
      $id_matchcount->{$_}++;
    }
    $matches_count += scalar @{ $id_matches };

  }
  return $matches_count;
}

# sub n_chunks{
#   my $self = shift;
#   return scalar @{ $self->chunk_specifiers() };
# }


###################

sub chunk_index_strings{	# 
  my $indices = shift;
  my $chunk_size = shift;
  #  print "Indices: ", join(", ", @$indices), "\n";
  my @chunk_specs = ();
  # print join(",", @$indices), "\n";
  while (scalar @$indices >= $chunk_size) {
    my $chunk_indices_str = join("_", splice( @$indices, 0, $chunk_size ));
    push @chunk_specs, $chunk_indices_str;
    # join("_", splice( @$indices, 0, $chunk_size )); # @chunk_indices);
  }
  #  print STDERR "XXX: ", join(";", @chunk_specs), "\n";
  return \@chunk_specs;
}

sub scrambled_chunks{
  my $sequence = shift;
  my $chunk_size = shift;
  my $indices = shift;	# array ref of indices: [13,2,45,5,17,6, ...] 

  my @chunk_seqs = ();
  my @chars = split('', $sequence);
  while (scalar @$indices >= $chunk_size) {
    my @chunk_indices = splice( @$indices, 0, $chunk_size);
    my $chunk_sequence = join('', @chars[@chunk_indices]);
    push @chunk_seqs, $chunk_sequence;
  }
  return \@chunk_seqs;
}

sub randomize_array{
  my $xs = shift;
  my $n_to_do = scalar @$xs;
  my $jswap = -1;     # where randomly chosen elem will be swapped to.
  while ($n_to_do > 1) {
    my $rand_index = int(rand($n_to_do));
    my $tmp = $xs->[$jswap];
    $xs->[$jswap] = $xs->[$rand_index];
    $xs->[$rand_index] = $tmp;
    $n_to_do--;
    $jswap--;
  }
  return $xs;
}

###########################################

__PACKAGE__->meta->make_immutable;

1;

