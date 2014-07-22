/*    Basic type definition for complex variables   */
/*                  6/19/91 cmn     	     	  */

#ifndef type_comp_done
#define type_comp_done

/* use built-in C99 complex type if available */
#include <complex.h>

/* Check for __COMPLEX_H__ as this is what is defined on some platforms */
/* such as on Mac OS X 10.7 */
#ifdef __COMPLEX_H__
#define _COMPLEX_H
#endif

#ifdef _COMPLEX_H
/* nothing to do here since complex is already defined */
#else

/* fallback: define our own complex type */
typedef struct complex_tag
{
  double re;
  double im;
} complex;

#endif


#endif
