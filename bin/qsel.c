#include <stdio.h>
#include <string.h>

int qselect(int* ix, int* v, int len, int k)
{
#	define SWAP(a, b) { tmp = v[a]; v[a] = v[b]; v[b] = tmp; }
  int i, st, tmp;
 
  for (st = i = 0; i < len - 1; i++) {
    if (v[i] > v[len-1]) continue;
    //    SWAP(i, st);
    tmp = v[i]; v[i] = v[st]; v[st] = tmp;
    tmp = ix[i]; ix[i] = ix[st]; ix[st] = tmp;
    st++;
  }
  tmp = v[len-1]; v[len-1] = v[st]; v[st] = tmp;
  tmp = ix[len-1]; ix[len-1] = ix[st]; ix[st] = tmp;
 
  return k == st	?v[st]
    :st > k	? qselect(ix, v, st, k)
    : qselect(ix + st, v + st, len - st, k - st);
}


int main(){

#	define N (sizeof(x)/sizeof(x[0]))
  int z[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 
      10, 11, 12, 13, 14, 15, 16, 17, 18, 19};
  int x[] = {9, 8, 7, 6, 5,
             10, 19, 13, 14, 18, 
             17, 11, 15, 16, 12, 
             0, 1, 2, 3, 4,};

 for(int i=0; i<N; i++){
      printf("%d  %d \n", z[i], x[i]);
 }printf("\n");

  int k = 13;
 
  int i;
  //  for (i = 0; i < 10; i++) {
  //   memcpy(y, x, sizeof(x)); // qselect modifies array
    printf("%d: %d\n", k, qselect(z, x, N, k));
    for(int i=0; i<N; i++){
      printf("%d  %d \n", z[i], x[i]);
    }
    // }
 
  return 0;
}

