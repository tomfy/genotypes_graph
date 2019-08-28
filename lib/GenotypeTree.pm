package GenotypeTree;
use Moose;
use namespace::autoclean;
use Carp;
use List::Util qw ( min max sum );
#use constant MISSING_DATA => '3';
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
print "$id  ", join(',', @ids_array), "\n";
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
   my $genotype_string = $gobj->sequence(); # entire (all snps) genotype as string.
   my $id = $gobj->id();
   my $root = $self->root();
   $root->add_id($id);
   my $ghead = substr($genotype_string, 0, 1);
print "Adding genotype to tree: $genotype_string \n";
   if (exists $root->children()->{$ghead}) {
      my $child = $root->children()->{$ghead};
         $child->add_genotype_compact($id, $genotype_string);
   } else {
      my $new_node = GenotypeTreeNode->new( {parent => $root, depth => 1, 
                                             genotype => $genotype_string, 
                                             ids => [$id] } );
      $root->children()->{$ghead} = $new_node;
   }
}

sub search{
   my $self = shift;
   my $gobj = shift;
   my $root = $self->root();
   print "root id: ", join(',', @{$root->ids()}), "\n";
   print "search id, genotype: ", $gobj->id(), "  ", $gobj->sequence(), "\n";
   my $matching_ids = $root->search_recursive($gobj->id(), $gobj->sequence());
   print "id: ", $gobj->id(), "  matches ids: $matching_ids \n";
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


