package GenotypeTree;
use Moose;
use namespace::autoclean;
use Carp;
use List::Util qw ( min max sum );
#use Readonly;

use constant BIG_NUMBER => 1_000_000_000;
#use constant MULTIPLIER => 1000;

has root => (
             isa => 'Object',   # should be a GenotypeTreeNode Object
             is => 'ro',
             default => sub { GenotypeTreeNode->new( { depth => 0, genotype => 'R' } ); },
            );

has depth => (                  # equal to length of sequences
              isa => 'Int',
              is => 'ro',
              required => 1,
             );


# class to represent a tree of genotypes

sub add_genotype{
   my $self = shift;
   my $gobj = shift ;
   my $genotype = $gobj->sequence(); # ArrayRef, e.g. [0, 0, 1, 2, 1, 2, 2, 0, 0, 0, 1, 0, '-', 1, 0];
   my $id = $gobj->id();
   my $root = $self->root();
   $root->inc_counter();
   my @ids_array = push @{$root->ids()}, $id;
   $root->ids( \@ids_array );

   my $current_node = $self->root();
   for my $g (@$genotype) {
      my $next_node =  $current_node->add_child($id, $g);
      #  print STDERR $current_node->as_string(), "   ", $next_node->as_string(), "\n";
      $current_node = $next_node;
   }
   return $current_node->depth() . " " . join(',', $current_node->ids());
}

sub add_genotype_compact{
   my $self = shift;
   my $gobj = shift;
   my $use_s = shift;
   my $genotype_string = $gobj->sequence(); # entire (all snps) genotype as string.
   my $id = $gobj->id();
   my $root = $self->root();

   my $ghead = substr($genotype_string, 0, 1);
   if (exists $root->children()->{$ghead}) {
      my $child = $root->children()->{$ghead};
      #  print STDERR "$genotype_string \n";
      if ($use_s) {
         $child->add_genotype_compact_s($id, $genotype_string);
      } else {
         $child->add_genotype_compact($id, $genotype_string);
      }
   } else {
      my $new_node = GenotypeTreeNode->new( {parent => $root, depth => 1, 
                                             genotype => $genotype_string, 
                                             #  genotype_ar => @$genotype,
                                             ids => [$id] } );
      $root->children()->{$ghead} = $new_node;
   }
}


sub as_newick{
   my $self = shift;
   return $self->root()->newick_recursive();
}

sub as_string{
   my $self = shift;
   my $leaves_only = shift // 0;
   my $string = ($leaves_only)? 
     $self->root()->leaves_as_string_recursive($leaves_only) : 
       $self->root()->as_string_recursive($leaves_only);
   return $string;
}



###########################################
__PACKAGE__->meta->make_immutable;

1;


