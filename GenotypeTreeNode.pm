package GenotypeTreeNode;
use Moose;
use namespace::autoclean;
use Carp;
use List::Util qw ( min max sum );


# a class for nodes of a tree.
has ids => (
            isa => 'ArrayRef',
            is => 'rw',
            default => sub { [] },
            );

has count => (
              traits => ['Counter'],
              isa => 'Int',
              is => 'rw',
              default => 0,
              handles => {
                          inc_counter => 'inc',
                         },
             );

has depth => (               # 0 for root, 1 for root's children, etc.
              isa => 'Int',
              is => 'ro',
              required => 1,
             );

has genotype => ( # the single-snp genotype associated with this node. 'R' for the root, otherwise 0, 1, 2, or '-'.
                 isa => 'Str',
                 is => 'ro',
                 required => 1,
                );

has parent => (
               isa => 'Maybe[Object]',
               is => 'ro',
               default => undef,
              );

has children => (
                 isa => 'HashRef', # keys are 0, 1, 2, or '-', values are GenotypeTreeNode objects
                 is => 'rw',
                 default => sub { {} },
                );

sub add_child{            # 
   my $self = shift;
   my $id = shift;
   my $g = shift;               # 0, 1, 2, or '-'
   my $d = $self->depth();      # depth of this node
   my $child_node;

   if (!defined $self->children()->{$g}) {
      $child_node = GenotypeTreeNode->new( {parent => $self, depth => $d+1, genotype => $g, ids => [$id] } );
      $self->children()->{$g} = $child_node;
   } else {
      $child_node = $self->children()->{$g};
      my @ids_array = push @{$child_node->ids()}, $id;
      $child_node->ids( \@ids_array );
   }
   $child_node->inc_counter();
   return $child_node;
}

sub as_string{
   my $self = shift;
   my $string = '';

   my @child_gs = keys  %{$self->children()};
   $string .= 'depth: ' . $self->depth() . '  ';
   $string .= 'count: ' . $self->count() . '  ';
   $string .= 'ids: ' . join(',', @{$self->ids()}) . '  ';
   $string .= 'genotype: ' . $self->genotype() . '  ';
   $string .= 'children: ' . join(' ', @child_gs) . "\n";
   return $string;
}

sub as_string_recursive{
   my $self = shift;
   my $string = $self->as_string();
   my @gs = sort keys %{$self->children()};
   print "child genotypes: ", join('  ', @gs), "\n"; 
   for my $g (@gs) {
      my $child_node = $self->children()->{$g};
      print "g: $g ", $child_node->genotype(), "  depths: ", $self->depth(), "  ", $child_node->depth(), "\n";
      $string .= $child_node->as_string_recursive();
   }
   return $string;
}

############################################

__PACKAGE__->meta->make_immutable;

1;
