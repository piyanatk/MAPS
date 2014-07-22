/*
 * A function that generates random numbers between 0 and 1.  The 
 * generator is seeded by passing a negative number and can be
 * reseeded at any point by passing another negative number.  This
 * was copied from Numerical Recipes in C.
 *
 * Sean Ting 09 July 05
 */


#include <stdlib.h>
#include <math.h>

long int random(void);
void srandom(unsigned int seed);


float ran0(idum)
     int *idum;

{
  static float y, maxran, v[7027];
  float dum;
  static int iff = 0;
  int j;
  unsigned i, k;

  if(*idum < 0 || iff == 0) {
    iff=1;
    i = 2;
    do {
      k=i;
      i<<=1;
    } while(i);
    maxran=k;
    srandom(*idum);
    *idum=1;
    for(j = 1; j <= 7027; j++) dum=random();
    for(j = 1; j <= 7027; j++) v[j]=random();
    y = random();
  }
  j = 1+7027.0*y/maxran;
  y=v[j];
  v[j]=random();
  return y/maxran;
}
