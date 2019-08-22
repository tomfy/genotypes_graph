package GenotypeTree;
use Moose;
use namespace::autoclean;
use Carp;
use List::Util qw ( min max sum );
#use Readonly;

use constant BIG_NUMBER => 1_000_000_000;
#use constant MULTIPLIER => 1000;

has root => (
             isa => 'Object', # should be a GenotypeTreeNode Object
             is => 'ro',
             default => sub { GenotypeTreeNode->new( { depth => 0, genotype => 'R' } ); },
            );

has depth => ( # equal to length of sequences
           isa => 'Int',
           is => 'ro',
           required => 1,
          );


# class to represent a tree of genotypes

sub add_genotype{
   my $self = shift;
   my $gobj = shift ;
   my $genotype = $gobj->sequence(); # e.g. (0, 0, 1, 2, 1, 2, 2, 0, 0, 0, 1, 0, '-', 1, 0);
   my $id = $gobj->id();
   my $root = $self->root();
   $root->inc_counter();
   my @ids_array = push @{$root->ids()}, $id;
   $root->ids( \@ids_array );

   my $current_node = $self->root();
   for my $g (@$genotype){
     my $next_node =  $current_node->add_child($id, $g);
       print STDERR $current_node->as_string(), "   ", $next_node->as_string(), "\n";
     $current_node = $next_node;
   }
   return $current_node->depth() . " " . join(',', $current_node->ids());
}

sub as_string{
   my $self = shift;
   my $string = $self->root()->as_string_recursive();
   return $string;
}



###########################################
__PACKAGE__->meta->make_immutable;

1;


