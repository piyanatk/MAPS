/*    Basic functions for complex variables   */
/*                  6/19/91 cmn     	     	  */

#include <math.h>
#include <stdio.h>
#include "comp_func.h"

#ifdef _COMPLEX_H

/* complex arithmetic already supported by compiler, or */

#else

double c_phase(complex z) 			/* Phase of complex number */
{
    if (z.re == 0.0)
       {
       if (z.im > 0.0)
	  return (1.570796327);
       else
	  return (-1.570796327);
       }
    else 
	return (atan2(z.im,z.re));
}

complex rect(double mag,double phase)		/* Reconstruct complex number from phase and magnitude */
{
    complex z;
    z.re = mag * cos(phase);
    z.im = mag * sin(phase);
    return(z);
}

void c_print(complex z) {
    if ((z.re==0.) && (z.im==0.)) printf("ZERO\n");
    else printf("%lf + %lf i ; %lf , %lf rad\n",z.re,z.im,c_mag(z),c_phase(z));
    fflush(stdout);
}

#endif
