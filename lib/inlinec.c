#define MISSING_DATA 'X'

/* number of mismatches between two strings up to some max  */
/* returns 0,1,2,...,max_mismatches */
void count_mismatches_C(char* str1, char* str2,
					    int max_mismatches, int mismatches_so_far,
					    SV* n_before_limit, SV* n_usable_to_limit, SV* n_mm){
  int i = 0;
  int n_usable_pair = 0; // count the number of good pairs (i.e. with neither being MISSING_DATA
  char c1;
  char c2;
  int mismatch_count = mismatches_so_far;
  while (c1 = str1[i]) {
    if (c1 != MISSING_DATA) {
      c2 = str2[i];
      if (c2 != MISSING_DATA) {
	n_usable_pair++;
	if (c1 != c2) {
	  mismatch_count++;
	}
	if (mismatch_count >= max_mismatches) {
	  break; // i not incremented
	}
      }
    }
    i++; //  the number of characters compared before hitting the mismatch limit. 
  } // end of while
  sv_setiv(n_before_limit, i);
  sv_setiv(n_usable_to_limit, n_usable_pair);
  sv_setiv(n_mm, mismatch_count);
}

void count_homozygous_mismatches_C(char* str1, char* str2,
						       int max_mismatches, int mismatches_so_far,
						       SV* n_before_limit, SV* n_usable_to_limit, SV* n_mm){
  int i = 0;
  int n_usable_pair = 0; // count the number of usable pairs (i.e. with neither being heterozygous or MISSING_DATA
  char c1;
  char c2;
  int mismatch_count = mismatches_so_far;
  while (c1 = str1[i]) {
    if ((c1 != '1') && (c1 != MISSING_DATA)) {
      c2 = str2[i];
      if ((c2 != '1') && (c2 != MISSING_DATA)) {
        n_usable_pair++;
	if (c1 != c2) {
	  mismatch_count++;
	}
	if (mismatch_count >= max_mismatches) {
	  break; // i not incremented
	}
      }
    }
    i++; // the number of characters compared before hitting the mismatch limit. 
  } // end of while
  sv_setiv(n_before_limit, i);
  sv_setiv(n_usable_to_limit, n_usable_pair);
  sv_setiv(n_mm, mismatch_count);
}


void agmr_C(char* str1, char* str2, SV* numer, SV* denom){
  char c1;
  char c2;
  int i=0;
  int n=0;
  int d=0;
  while(c1 = str1[i]){
    if(c1 != MISSING_DATA){
      c2 = str2[i];
      if(c2 != MISSING_DATA){
	d++;
	if(c1 != c2){
	  n++;
	}
      }
    }
    i++;
  }
  sv_setiv(numer, n);
  sv_setiv(denom, d);
}

void hgmr_Cx(char* str1, AV* qindices, char* str2, SV* numer, SV* denom){ // slower!
  /* qindices is array of the indices for which str1 (query) is homozygous */
  char c1;
  char c2;
  int n=0;
  int d=0;
  int L = av_len(qindices);
  for(int i=0; i<L; i++){
    int idx = SvNV( * av_fetch( qindices, i, 0) );
    c1 = str1[idx];
    c2 = str2[idx];
    if((c2 != '1') && (c2 != MISSING_DATA)){
      d++;
      if(c1 != c2){
	n++;
      }
    }
  }
  sv_setiv(numer, n);
  sv_setiv(denom, d);
}

void hgmr_C(char* str1, char* str2,
	    SV* numer, SV* denom){
  char c1;
  char c2;
  int i=0;
  int n=0;
  int d=0;
  while(c1 = str1[i]){
    if((c1 != '1') && (c1 != MISSING_DATA)){
      c2 = str2[i];
      if((c2 != '1') && (c2 != MISSING_DATA)){
	d++;
	if(c1 != c2){
	  n++;
	}
      }
    }
    i++;
  }
  sv_setiv(numer, n);
  sv_setiv(denom, d);
}

void agmr_hgmr_C(char* str1, char* str2, SV* a_numer, SV* a_denom, SV* h_numer, SV* h_denom){
  char c1;
  char c2;
  int i=0;
  int n_a=0;
  int r_a=0;
  int n_h=0;
  int r_h=0;
  while(c1 = str1[i]){
    if(c1 != MISSING_DATA){
      c2 = str2[i];
      if(c2 != MISSING_DATA){
	if(c1 != c2){
	  n_a++;
	  if((c1 != '1')  &&  (c2 != '1')){
	    n_h++;
	  }
	}else{
	  r_a++;
	  if((c1 != '1')  &&  (c2 != '1')){
	    r_h++;
	  }
	}

      }
    }
    i++;
  }
  sv_setiv(a_numer, n_a);
  sv_setiv(a_denom, r_a+n_a);
  sv_setiv(h_numer, n_h);
  sv_setiv(h_denom, r_h+n_h);
}




