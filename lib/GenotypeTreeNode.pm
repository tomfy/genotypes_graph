package GenotypeTreeNode;
use Moose;
use namespace::autoclean;
#use Carp::Assert;
use List::Util qw ( min max sum );
use constant MISSING_DATA => 'X';

no warnings 'recursion';

#use constant DEBUG => 0;


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

has depth => ( # 0 for root, counts not number of nodes below root, but number of snps, so all leaves have depth == lenght of genotype sequences.
              isa => 'Int',
              is => 'ro',
              required => 1,
             );

has genotype => ( # the one or more snp genotype associated with this node. Characters are 'R' for the root, otherwise 0, 1, 2, or '-'.
                 isa => 'Str',
                 is => 'rw',
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

sub add_child{                  # 
   my $self = shift;
   my $id = shift;
   my $g = shift;               # string of 0, 1, 2, or '-'
   my $d = $self->depth();      # depth of this node
   my $child_node;

   if (!defined $self->children()->{$g}) {
      $child_node = GenotypeTreeNode->new( {parent => $self, depth => $d+1, genotype => $g, ids => [$id] } );
      $self->children()->{$g} = $child_node;
   } else {
      $child_node = $self->children()->{$g};
      push @{$child_node->ids()}, $id;
   }
   $child_node->inc_counter();
   return $child_node;
}

sub add_genotype_compact{
   my $self = shift;
   my $id2 = shift;
   my $gt2s = shift;
   my $gt1s = $self->genotype();
   my $d = $self->depth();

   my $n_equal_characters = 0;
   for my $i (0 .. (length $gt1s) - 1) {
      last if(substr($gt1s, $i, 1) ne substr($gt2s, $i, 1));
      $n_equal_characters++;
   }
#   assert ($n_equal_characters > 0) if DEBUG;

   if ($n_equal_characters < length $gt1s) { # $gt1s and beginning of $gt2s are not the same.
      my $g_common = substr($gt1s, 0, $n_equal_characters, '');
      my $g2_common = substr( $gt2s, 0, $n_equal_characters, ''); 
#     assert ($g2_common eq $g_common) if DEBUG;
      my @new_child1_ids = @{ $self->ids() };
      my $new_child1 = GenotypeTreeNode->new( { parent => $self, 
                                                depth => $d+$n_equal_characters, 
                                                genotype => $gt1s, 
                                                ids => \@new_child1_ids,
                                                children => $self->children() } );
      my $new_child2 = GenotypeTreeNode->new( { parent => $self, 
                                                depth => $d+$n_equal_characters, 
                                                genotype => $gt2s, 
                                                ids => [$id2] } );
      $self->genotype($g_common);
      push @{$self->ids()}, $id2;
      $self->children( { substr($gt1s, 0, 1) => $new_child1, substr($gt2s, 0, 1) => $new_child2 } );
   } else {              # $gt1s and beginning of $gt2 are the same, so
      substr($gt2s, 0, $n_equal_characters, '');
      if ((length $gt2s) > 0) { # theres a bit of $gt2s left 
         my $gt2head = substr($gt2s, 0, 1);
         if (exists $self->children()->{$gt2head}) {
            $self->children()->{$gt2head}->add_genotype_compact($id2, $gt2s);
         } else {
            my $new_child = 
              GenotypeTreeNode->new( { parent => $self, depth => $d+$n_equal_characters, 
                                       genotype => $gt2s,
                                       ids => [$id2] } );
            $self->children()->{$gt2head} = $new_child;
         }
      } else {                  # this must be a leaf node
         push @{$self->ids()}, $id2;
      }
   }
}

sub search_recursive{
   my $self = shift;
   my $id2 = shift;
   my $gt2s = shift;            # searching for this one.
   my $gt1s = $self->genotype();
   my ($L1, $L2) = (length $gt1s, length $gt2s);

   my $id1s = join(",", @{$self->ids()});
#   assert ($L2 >= $L1) if DEBUG;
   my $matching_id_string = '';
   my $n_equal_characters = 0;
   for my $i (0 .. $L1 - 1) {
      my ($g1, $g2) = (substr($gt1s, $i, 1), substr($gt2s, $i, 1));
      last if($g1 ne MISSING_DATA  and  $g2 ne MISSING_DATA  and  $g2 ne $g1 );
      $n_equal_characters++;
   }
   if ($L2 == $L1) {     # moment of truth - this must be a leaf node.
      if ($n_equal_characters == $L1) { # these are identical genotypes!!
         $matching_id_string .= join(',', @{$self->ids()}) . ',';
      } else { # the genotype searched for is not present in this branch.
      }
   } elsif ($L2 > $L1) {        # Ok so far, but must look further ...
      substr($gt2s, 0, $n_equal_characters, '');
      my $g2head = substr($gt2s, 0, 1);
      if ($g2head eq MISSING_DATA) { # must check all subtrees
         while ( my ($gh, $child) = each %{ $self->children() }) {
            $matching_id_string .= $child->search_recursive($id2, $gt2s);
         }
      } else {
         for my $gh ($g2head,MISSING_DATA) {
            my $child = $self->children()->{$gh} // undef;
            if (defined $child) {
               if ($gh eq $g2head  or  $gh eq MISSING_DATA) {
                  $matching_id_string .= $child->search_recursive($id2, $gt2s);
               }
            }
         }
      }
   }
   return $matching_id_string;
}

sub newick_recursive{
   my $self = shift;
   my $newick_str = '';

   if (keys %{ $self->children() } ) { # not a leaf node.
      $newick_str .= '(';
      for my $g (0,1,2,MISSING_DATA) {
         $newick_str .= ( exists $self->children()->{$g} )? $self->children()->{$g}->newick_recursive() . ',' : '';
      }
      $newick_str =~ s/,$//;
      $newick_str .= ')' . $self->genotype() . ':1';
   } else {                     # leaf
      my $gstring = $self->genotype();
      substr($gstring, 4, -4, '...') if(length $gstring > 8);
      $newick_str .= join('-', @{$self->ids()}) . '_' . $gstring . ':1';
   }
   return $newick_str;
}

sub as_string{
   my $self = shift;
   my $leaves_only = shift;
   my $string = '';
   my @child_gs = keys  %{$self->children()};
   if ($leaves_only) {
      if (scalar @child_gs == 0) {
         $string .= 'd: ' . $self->depth() . '  ';
         $string .= 'g: ' . $self->genotype() . '  ';
         $string .=  $self->count() . ' ';
         $string .=  join(',', @{$self->ids()}) . "\n";
      }
   } else {
      $string .= 'depth: ' . $self->depth() . '  ';
      $string .= 'count: ' . $self->count() . '  ';
      $string .= 'ids: ' . join(',', @{$self->ids()}) . '  ';
      $string .= 'genotype: ' . $self->genotype() . '  ';
      $string .= 'children: ' . join(' ', @child_gs) . "\n";
   }
   return $string;
}

sub as_string_recursive{
   my $self = shift;
   my $leaves_only = shift;
   my $string = $self->as_string($leaves_only);
   my @gs = sort keys %{$self->children()};
   #  print "child genotypes: ", join('  ', @gs), "\n"; 
   for my $g (@gs) {
      my $child_node = $self->children()->{$g};
      #   print "g: $g ", $child_node->genotype(), "  depths: ", $self->depth(), "  ", $child_node->depth(), "\n";
      $string .= $child_node->as_string_recursive($leaves_only);
   }
   return $string;
}

sub leaves_as_string{
   my $self = shift;
   my $string = '';
   my @child_gs = keys  %{$self->children()};
   if (scalar @child_gs == 0) {
      $string .= 'd: ' . $self->depth() . '  ';
      $string .= 'g: ' . $self->genotype() . '  ';
      $string .=  $self->count() . ' ';
      $string .=  join(',', @{$self->ids()}) . "\n";
   }
   return $string;
}

sub leaves_as_string_recursive{
   my $self = shift;
   my $string = $self->leaves_as_string();
   my @gs = sort keys %{$self->children()};
   #  print "child genotypes: ", join('  ', @gs), "\n"; 
   for my $g (@gs) {
      my $child_node = $self->children()->{$g};
      #   print "g: $g ", $child_node->genotype(), "  depths: ", $self->depth(), "  ", $child_node->depth(), "\n";
      $string .= $child_node-> leaves_as_string_recursive();
   }
   return $string;
}

sub add_id{
my $self = shift;
my $id_to_add = shift;

my @ids = @{$self->ids()};
push @ids, $id_to_add;
$self->ids(\@ids);
}

sub id_as_string{
my $self = shift;
return join(',', @{$self->ids()});
}

############################################

__PACKAGE__->meta->make_immutable;

1;
