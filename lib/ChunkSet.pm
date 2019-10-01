package ChunkSet;
#use Moose;
use Mouse;
use namespace::autoclean;
use Carp::Assert;
use List::Util qw ( min max sum );
#use Inline 'C';
#use constant MISSING_DATA => 'X';
#use constant BIG_NUMBER => 1_000_000_000;

# a class for genotypes
# each has an identifier (integer),
# and a genotype represented as a string
# containing characters 0, 1, 2, MISSING_DATA
# 0 = AA, 1 = Aa, 2 = aa.

# has chunk_size => (
# 		   isa => 'Int',
# 		   is => 'rw',
# 		   required => 0,
# 		  );

has chunk_specifiers => ( # array ref of strings, e.g. ['1_22_108_23', ...
			 isa => 'ArrayRef',
			 is => 'rw',
			 required => 1,
			);

has chunk_spec_arrays => ( # array ref of array refs, e.g. [ [1,22,108,23], ...
			 isa => 'ArrayRef',
			 is => 'rw',
			 required => 0,
			);

has seqid_seq => (
		  isa => 'HashRef',
		  is => 'rw',
		  required => 1,
		 );

has chunkspec__seq_ids => (
			   isa => 'HashRef',
			   is => 'rw',
			   default => sub { {} },
			  );

sub BUILD{
  my $self = shift;

  for (@{$self->chunk_specifiers()}) {
    $self->chunkspec__seq_ids()->{$_} = {};
  #  print "chunk specifier: $_ \n";
  }

  my @chunk_spec_arrays = map( [split('_', $_)], @{$self->chunk_specifiers()} );
  $self->chunk_spec_arrays( [map( [split('_', $_)], @{$self->chunk_specifiers()} )] ); 
#  $self->chunk_size( scalar @{ $chunk_spec_arrays[0] } );
  
  while ( my($seqid, $seq) = each %{$self->seqid_seq()} ) { # loop over sequences
    my @seq_chars = split('', $seq);
    while ( my ($ich, $ch_indices) = each @chunk_spec_arrays) { # index, and array ref.
      my $ch_spec = $self->chunk_specifiers()->[$ich]; # as a string
      my $chunk_seq = @seq_chars[@$ch_indices];
      $self->chunkspec__seq_ids()->{$ch_spec}->{$chunk_seq} //= [];
      push @{ $self->chunkspec__seq_ids()->{$ch_spec}->{$chunk_seq} }, $seqid;
    }

  }
}

sub get_chunk_match_counts{
  my $self = shift;
  my $sequence = shift; # the sequence to match
  my @seq_chars = split('', $sequence);
  my %ids_matchcounts = ();
    while ( my ($ich, $ch_indices) = each @{$self->chunk_spec_arrays()}) { # index, and array ref.
      my $ch_spec = $self->chunk_specifiers()->[$ich]; # as a string
    #  print STDERR "ch spec: $ch_spec ", join(',', @$ch_indices), "\n";
      my $chunk_seq = @seq_chars[@$ch_indices];
      my $matching_ids = $self->chunkspec__seq_ids()->{$ch_spec}->{$chunk_seq} // []; # array ref of ids matching sequence in this chunk.
      for(@$matching_ids){
	$ids_matchcounts{$_} += 1;
      }

#      $self->chunkspec__seq_ids()->{$ch_spec}->{$chunk_seq} //= [];
#      push @{ $self->chunkspec__seq_ids()->{$ch_spec}->{$chunk_seq} }, $seqid;
    }
  return \%ids_matchcounts;
}

  ###################

sub scrambled_chunks{
      my $sequence = shift;
      my $chunk_size = shift;
      my $indices = shift; # array ref of indices: [13,2,45,5,17,6, ...] 

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
      my $jswap = -1; # where randomly chosen elem will be swapped to.
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

