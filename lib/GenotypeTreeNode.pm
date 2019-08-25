package GenotypeTreeNode;
use Moose;
use namespace::autoclean;
use Carp::Assert;
use List::Util qw ( min max sum );
no warnings 'recursion';

# use constant DEBUG => 0;


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

has genotype => ( # the one or more snp genotype associated with this node. Characters are 'R' for the root, otherwise 0, 1, 2, or '-'.
                 isa => 'Str',
                 is => 'rw',
                 required => 1,
                );

# has genotype_ar => ( # the one or more snp genotype associated with this node. Characters are 'R' for the root, otherwise 0, 1, 2, or '-'.
#                  isa => 'Array',
#                  is => 'rw',
#                  required => 0,
#                 );

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

# sub BUILD{
#    my $self = shift;
#    my @xs = split('', $self->genotype());
#    print STDERR "[" , ref @xs, "]   ", scalar @xs, "  ", join("; ", @xs), "  ", $xs[0], "\n";
#    $self->genotype_ar( split('', $self->genotype()) );
#    return $self;
# }

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
   my $gt2add = shift;
   my @gt2 = split('', $gt2add);      # new multi-snp genotype ArrayRef
#   my @gt1 = @{ $self->genotype_ar() }; 
   my @gt1 = split('', $self->genotype()); # just the snps associated with this node.
   my $d = $self->depth();
   my @g_common = undef;        # the initial part
#   my $i_first_not_equal = undef;
   my $n_equal_characters = 0;
   while (my ($i, $g1) = each @gt1) { # find first non-equal character. Maybe use Inline::C ?
      last if ($gt2[$i] ne $g1);
      $n_equal_characters++;
   }
 #  die "In add_genotype_compact. first chars are not equal.\n" if($n_equal_characters == 0);
   assert ($n_equal_characters > 0) if DEBUG;
   if ($n_equal_characters < scalar @gt1) { # @gt1 and beginning of @gt2 are not the same.
      @g_common = splice( @gt1, 0, $n_equal_characters); # common to start of both @gt1 and @gt2
      splice( @gt2, 0, $n_equal_characters);             # 
      my @new_child1_ids = @{ $self->ids() };
      my $new_child1 = GenotypeTreeNode->new( { parent => $self, 
                                                depth => $d+$n_equal_characters, 
                                                genotype => join('', @gt1), 
                                              #  genotype_ar => @gt1,
                                                ids => \@new_child1_ids, 
                                                children => $self->children() } );
      my $new_child2 = GenotypeTreeNode->new( { parent => $self, 
                                                depth => $d+$n_equal_characters, 
                                                genotype => join('', @gt2), 
                                           #     genotype_ar => @gt2,
                                                ids => [$id2] } );
      $self->genotype(join('', @g_common));
      push @{$self->ids()}, $id2;
      $self->children( { $gt1[0] => $new_child1, $gt2[0] => $new_child2 } ); 
   } else {              # @gt1 and beginning of @gt2 are the same, so
   #   print STDERR "YYY: ", join('', @gt1), " ", join('', @gt2), "\n";
      splice (@gt2, 0, $n_equal_characters); # delete $n_equal_characters characters from beginning of @gt2
      if (@gt2) {                            # theres a bit of @gt2 
       #  print STDERR "XXX: ", join("; ", @gt2), "\n";
         if (exists $self->children()->{$gt2[0]}) {
            $self->children()->{$gt2[0]}->add_genotype_compact($id2, join('', @gt2));
         } else {
            my $new_child = 
              GenotypeTreeNode->new( { parent => $self, depth => $d+$n_equal_characters, 
                                       genotype => join('', @gt2), 
                                   #    genotype_ar => @gt2,
                                       ids => [$id2] } );
            $self->children()->{$gt2[0]} = $new_child;
         }
      } else {                  # this must be a leaf node
         push @{$self->ids()}, $id2;
      }
   }
   
}

sub add_genotype_compact_s{
   my $self = shift;
   my $id2 = shift;
   my $gt2s = shift;
 #  my @gt2 = @{ $gt2add };      # new multi-snp genotype ArrayRef
 #  my $gt2s = join('', @$gt2add);
#   my @gt1 = @{ $self->genotype_ar() }; 
#   my @gt1 = split('', $self->genotype()); # just the snps associated with this node.
   my $gt1s = $self->genotype();
   my $d = $self->depth();

 #  my $i_first_not_equal = undef;
   my $n_equal_characters = 0;
   for my $i (0 .. (length $gt1s) - 1) {
      last if(substr($gt1s, $i, 1) ne substr($gt2s, $i, 1));
      $n_equal_characters++;
   }
   assert ($n_equal_characters > 0) if DEBUG;

   my $g_common;
   if ($n_equal_characters < length $gt1s) { # @gt1 and beginning of @gt2 are not the same.
     #  = splice( @gt1, 0, $n_equal_characters); # common to start of both @gt1 and @gt2
    #  print "$gt1s \n", "$gt2s \n";
       $g_common = substr($gt1s, 0, $n_equal_characters, '');
      my $g2_common = substr( $gt2s, 0, $n_equal_characters, '');             # 
  #    print "$g_common $gt1s \n";
   #   print "$g2_common $gt2s \n\n";
       assert ($g2_common eq $g_common) if DEBUG;
      my @new_child1_ids = @{ $self->ids() };
      my $new_child1 = GenotypeTreeNode->new( { parent => $self, 
                                                depth => $d+$n_equal_characters, 
                                                genotype => $gt1s, 
                                             #   genotype_ar => \@gt1,
                                                ids => \@new_child1_ids,
                                                children => $self->children() } );
      my $new_child2 = GenotypeTreeNode->new( { parent => $self, 
                                                depth => $d+$n_equal_characters, 
                                                genotype => $gt2s, 
                                             #   genotype_ar => \@gt2,
                                                ids => [$id2] } );
      $self->genotype($g_common);
      push @{$self->ids()}, $id2;
      $self->children( { substr($gt1s, 0, 1) => $new_child1, substr($gt2s, 0, 1) => $new_child2 } );
   } else {              # @gt1 and beginning of @gt2 are the same, so
   #   print STDERR "YYY: ", join('', @gt1), " ", join('', @gt2), "\n";
      #splice (@gt2, 0, $n_equal_characters); # delete $n_equal_characters characters from beginning of @gt2
      substr($gt2s, 0, $n_equal_characters, '');
      if ((length $gt2s) > 0) {                            # theres a bit of @gt2 
       #  print STDERR "XXX: ", join("; ", @gt2), "\n";
         my $gt2head = substr($gt2s, 0, 1);
         if (exists $self->children()->{$gt2head}) {
            $self->children()->{$gt2head}->add_genotype_compact_s($id2, $gt2s);
         } else {
            my $new_child = 
              GenotypeTreeNode->new( { parent => $self, depth => $d+$n_equal_characters, 
                                       genotype => $gt2s, 
                                     #  genotype_ar => \@gt2,
                                       ids => [$id2] } );
            $self->children()->{$gt2head} = $new_child;
         }
      } else {                  # this must be a leaf node
         push @{$self->ids()}, $id2;
      }
   }
   
}


sub newick_recursive{
   my $self = shift;

   my $newick_str = '';

   if (keys %{ $self->children() } ) { # not a leaf node.
      $newick_str .= '(';
      for my $g (0,1,2,'-') {
         $newick_str .= ( exists $self->children()->{$g} )? $self->children()->{$g}->newick_recursive() . ',' : '';
      }
      $newick_str =~ s/,$//;
      $newick_str .= 
     #   join('-', @{$self->ids()}) . '_' .
          #  $self->genotype() . 
              ')' . $self->genotype() . ':1';
   } else {                     # leaf
      $newick_str .= join('-', @{$self->ids()}) . '_' . $self->genotype() . ':1';
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

############################################

__PACKAGE__->meta->make_immutable;

1;
