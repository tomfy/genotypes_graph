package Genotype;
#use Moose;
use Mouse;
use namespace::autoclean;
use Carp::Assert;
use List::Util qw ( min max sum );
use Inline 'C';
use constant MISSING_DATA => 'X';
use constant BIG_NUMBER => 1_000_000_000;

# a class for genotypes
# each has an identifier (integer),
# and a genotype represented as a string
# containing characters 0, 1, 2, MISSING_DATA
# 0 = AA, 1 = Aa, 2 = aa.

has id => (
           isa => 'Int',
           is => 'ro',
           required => 1,
          );

has generation => (
                   isa => 'Int',
                   is => 'ro',
                   default => sub { undef },
                  );

has pedigree => (
                 isa => 'Str',
                 is => 'ro',
                 default => sub { undef },
                );

has sequence => (
		 isa => 'Str',
		 is => 'ro',
		 required => 1,
                );

around BUILDARGS => sub {
   my $orig = shift;
   my $class = shift;

   if (@_ == 1) {
      if (!ref $_[0]) { # argument should be string with >id generation pedigree \n sequence.
         my $arg = $_[0];
         my ($idline, $sequence_string) = split(/\n/, $arg);
         my ($id, $generation, $pedigree) = (undef, undef, undef);
         if ($idline =~ /^>(\S+)\s+(\d+)\s+(\S+)\s*$/x) {
            ($id, $generation, $pedigree) = ($1, $2, $3);
         } elsif ($idline =~ /^ > (\S+) /x) {
            $id = $1;
         } else {
            die "in Genotype BUILDARGS argument $arg does not have expected >id format.\n";
         }
         $sequence_string =~ s/\s+//gx;
         #   my @sequence_array = split(//, $sequence_string);
         return {id => $id, sequence => $sequence_string, generation => $generation, pedigree => $pedigree};
      } else {
         return $_[0];		# should be hash ref of arguments
      }
   }
   # otherwise, don't modify arguments, should be array ref with id, sequence, and optionally generation and pedigree
};

sub clone{
   my ($self, %params) = @_;
   return $self->meta->clone_object($self, %params);
}


sub distance{ # calculate distance between this genotype obj. and another
   my $self = shift;
   my $other_genotype = shift;
   my $this_gt = $self->sequence();	      # string
   my $other_gt = $other_genotype->sequence(); # string

   my $Dcalc = 1;

   if ($Dcalc == 1) {		       # use inline::C function: (much faster)
      my ($d, $agmr, $hgmr, $ogmr);
      distances_C(0, BIG_NUMBER, $this_gt, $other_gt, $d, $agmr, $hgmr, $ogmr);
#      printf("%10.6f %10.6f, %10.6f, %10.6f \n", $d, $agmr, $hgmr, $ogmr);
      return ($d, $agmr, $hgmr, $ogmr);
    }elsif($Dcalc == 2){
      my ($d, $agmr, $hgmr, $ogmr) = (distance_C($this_gt, $other_gt), 5, 5, 5);
      return ($d, $agmr, $hgmr, $ogmr);
   } else {			# pure perl. much slower.
      my $distance = 0;
      my $count_both = 0; # count of snps with data present in both sequences
      my $count_missing = 0; # count of snps with data absent in one or both sequences
      if (length $this_gt == length $other_gt) {
         for my $i (0 .. (length $this_gt) - 1) {
            my $c1 = substr($this_gt, $i, 1);
            my $c2 = substr($other_gt, $i, 1); # $other_gt->[$i];
            if ( ($c1 eq MISSING_DATA)  or ($c2 eq MISSING_DATA) ) {
               $count_missing++;
            } else {
               $count_both++;
               if ($c2 != $c1) {
                  if ($c1 == 0) {
                     $distance += $c2 - $c1;
                  } elsif ($c1 == 1) {
                     $distance += 1;
                  } elsif ($c1 == 2) {
                     $distance += $c1 - $c2;
                  } else {
                     die "c1 has unhandled value: $c1 \n";
                  }
               }
            }
         }
      }
      return ($count_both > 0)? $distance/$count_both : 1000;
   }
}

sub count_mismatches_up_to_limit{
   my $self = shift;
   my $other_genotype = shift;
   my $max_mismatches = shift;
   my ($n_ok, $n_mm);
   count_mismatches_up_to_limit_C($self->sequence(), $other_genotype->sequence(), $max_mismatches, $n_ok, $n_mm);
   # print "aaa: $n_ok  $n_mm \n";
   return ($n_ok, $n_mm);
}

sub mean{
   my $self = shift;
   my $other_genotype = shift;
   my $this_gt = $self->sequence();	      # array ref of 0,1,2,3
   my $other_gt = $other_genotype->sequence(); # array ref of 0,1,2,3
   my @mean_gt = ();
   if ( length $this_gt == length $other_gt) { # check that lengths are equal
      for my $i (0 .. (length $this_gt) - 1) {
         my $c1 = substr($this_gt, $i, 1);
         my $c2 = substr($other_gt, $i, 1);
         if ($c1 eq MISSING_DATA) {
            push @mean_gt, $c2;
         } elsif ($c2 eq MISSING_DATA) {
            push @mean_gt, $c1;
         } else {
            push @mean_gt, 0.5*($c1 + $c2);
         }
      }
   }
   return \@mean_gt;
}

sub add_noise{
   my $self = shift;
   my $error_prob = shift;
   my $sequence = $self->sequence();
   if (1) {
      #  print "$sequence \n";
      for my $i (0 .. (length $sequence)-1) {
         #   print "$i  ", rand(1), "\n";
         if (rand(1) < $error_prob) { # error
            #    print "   ***\n";
            my $c = substr($sequence, $i, 1);
            if ($c == 0) {
               if (rand(1) < 0.5) {
                  $c = 1;
               } else {
                  $c = 2;
               }
            } elsif ($c == 1) {
               if (rand(1) < 0.5) {
                  $c = 0;
               } else {
                  $c = 2;
               }
            } elsif ($c == 2) {
               if (rand(1) < 0.5) {
                  $c = 0;
               } else {
                  $c = 1;
               }
            }
            substr($sequence, $i, 1, $c);
         }                      # end if error branch
      }
   } else {			# error model 2
      my $trans_prob1 = $error_prob*(1.0 - $error_prob);
      my $trans_prob2 = $error_prob*$error_prob;
      for my $i (0 .. (length $sequence)-1) {
         my $c = substr($sequence, $i, 1);
         my $rand_number = rand();
         if ($c == 1) {
            if ($rand_number < $trans_prob1) {
               $c = 0;
            } elsif ($rand_number < 2*$trans_prob1) {
               $c = 2;
            }
         } elsif ($c == 0) {
            if ($rand_number < 2*$trans_prob1) {
               $c = 1;
            } elsif ($rand_number < 2*$trans_prob1 + $trans_prob2) {
               $c = 2;
            }
         } elsif ($c == 2) {
            if ($rand_number < 2*$trans_prob1) {
               $c = 1;
            } elsif ($rand_number < 2*$trans_prob1 + $trans_prob2) {
               $c = 0;
            }
         }
      }
   }
   $self->{sequence} = $sequence;
}
  


############################################

__PACKAGE__->meta->make_immutable;

1;

__DATA__
############################################
########### inline C stuff #################
__C__

  /* distance between two genotype strings  */
/* 0-1, 1-2 -> d=1; 0-2 -> d=2  */
double distance_C(char* str1, char* str2) {
  int i = 0;
  char c1;
  char c2;
  int dist = 0;
  int count = 0;
  while(c1 = str1[i]) {
    c2 = str2[i];

    if (c1 == '0') {
      if(c2 == '0'){
        count++;
      }else if(c2 == '1'){
        dist++;
        count++;
      }else if(c2 == '2'){
        dist += 2;
        count++;
      } // else c2 is missing data
} else if(c1 == '1'){

  if (c2 == '0') {
  dist++;
count++;
} else if(c2 == '1'){
          count++;
       } else if(c2 == '2'){
                 dist += 1;
                 count++;
              }

} else if(c1 == '2'){

   if (c2 == '0') {
       dist += 2;
       count++;
    } else if(c2 == '1'){
              dist++;
              count++;
           } else if(c2 == '2'){
                     count++;
                  }

}

i++;
}
//  printf("%d  %d \n", dist, count);
if (count > 0) {
   //  printf("%g\n", 1.0*dist/count);
   return 1.0*dist/count;
} else
  return 1000;
}



double distances_C(int start, int n_homozyg, char* str1, char* str2, SV* dist, SV* agmr, SV* hgmr, SV* ogmr) {
     // calculates 3 kinds of distance: dist, agmr and hgmr
     // start at snp with index start, go until end or until n_homozyg pairs with both homozygous.
       int i = start;
     char c1;
     char c2;
    
       int count_02 = 0; // 0<->2 (both directions)
         int count_0022 = 0; // 0<->0 or 2<->2
           int count_11 = 0; // 1<->1
             int count_1x = 0; // 1<->0 or 1<->2
       int hcount = 0; // count only snps homozygous in both
         while (c1 = str1[i]) {
            c2 = str2[i];

            if (c1 == '0') {

               if (c2 == '0') {
                   count_0022++;
                } else if(c2 == '1'){
                          count_1x++;
                       } else if(c2 == '2'){
                                 count_02++;
                              }
               // else c2 is missing data
            } else if(c1 == '1'){

               if (c2 == '0') {
                  count_1x++;
                } else if(c2 == '1'){
                          count_11++;
                       } else if(c2 == '2'){
                                count_1x++;
                              }

            } else if(c1 == '2'){

               if (c2 == '0') {
                  count_02++;
                } else if(c2 == '1'){
                          count_1x++;
                       } else if(c2 == '2'){
                                 count_0022++;
                              }
            }
            i++;
//	     printf("%d %d %d \n", start, n_homozyg, count_0022 + count_02);
            if((count_0022 + count_02) >= n_homozyg){ break; }  // stop after n_homozyg homozygous-only pairs
         } // end of while

         double homozyg_denom = count_0022 + count_02; // both genotypes of pair homozygous
         double other_denom = count_11 + count_1x; // 0 or 1 genotypes of pair are homozygous ( no missing data)
         double all_denom = homozyg_denom + other_denom; // all except if there is missing data

//	 printf("%8d %8d %8d %8d   %8.1f\n", count_0022, count_02, count_11, count_1x, all_denom);
     double homozyg_dist = (homozyg_denom > 0)? count_02/homozyg_denom : 1000;
     double other_dist = (other_denom > 0)? count_1x/other_denom : 1000;
     double all_dist = (all_denom > 0)? (count_1x + count_02)/all_denom : 1000;
     double d_dist = (all_denom > 0)? (count_1x + 2*count_02)/all_denom : 1000;

     sv_setnv(hgmr, homozyg_dist);
     sv_setnv(ogmr, other_dist);
     sv_setnv(agmr, all_dist);
     sv_setnv(dist, d_dist);

  }


/* number of mismatches between to strings up to some max  */
  /* returns 0,1,2,...,max_mismatches or -1 if there are more mismatches */
void count_mismatches_up_to_limit_C(char* str1, char* str2, int max_mismatches, SV* n_ok, SV* n_mm) {
  int i = 0;
  char c1;
  char c2;
  int mismatch_count = 0;
  while(c1 = str1[i]) {
     c2 = str2[i];
     if (c1 != c2) {
       mismatch_count++;
     }
     if(mismatch_count > max_mismatches){
       mismatch_count = -1;
       break;
     }
     i++;
/* i is now the number of characters which have been compared without exceeding the mismatch limit. */
  }
 /* Inline_Stack_Vars;
Inline_Stack_Reset;
Inline_Stack_Push(sv_2mortal(newSVuv(i)));
Inline_Stack_Push(sv_2mortal(newSVuv(mismatch_count)));
Inline_Stack_Done; */

  sv_setiv(n_ok, i);
  sv_setiv(n_mm, mismatch_count);

}


