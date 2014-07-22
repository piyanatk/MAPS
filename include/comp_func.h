/* public function prototypes for functions in "comp_func.c"
created Jan 2007. Randall Wayth. */


#ifndef COMP_FUNC_H
#define COMP_FUNC_H

#include "type_comp.h"
#include <math.h>

/* if C99 complex is already included, then just create inline wrapper 
 * functions using standard C99 complex */
#ifdef _COMPLEX_H

static inline double c_mag (complex x) {
    return cabs(x);
}
static inline double c_phase(complex x) {
    return carg(x);
}
static inline complex c_add(complex a, complex b) {
    return (a+b);
}
static inline complex c_sub(complex a, complex b) {
    return (a-b);
}
static inline complex c_mult(complex a, complex b) {
    return a*b;
}
static inline complex rect(double mag, double ph) {
    return cos(ph)*mag + I*mag*sin(ph);
}
static inline complex c_conj(complex z) {
    return conj(z);
}
static inline complex c_exp(double theta) {
    return cos(theta) + I*sin(theta);
}
static inline complex c_zero() {
    return 0.0;
}
static inline complex s_mult(complex z, double i) {return z*i; }
static inline double  c_real(complex z) { return creal(z); }
static inline double  c_imag(complex z) { return cimag(z); }
static inline complex c_make(double re, double im) {return (re + I*im);}

#else
/* use exising (deprecated) types. Changed to inline for maximum efficiency. */
static inline double c_mag(complex z) {    return (sqrt(z.re * z.re + z.im * z.im));}
double c_phase(complex z);
static inline complex c_add(complex a, complex b) { complex z; z.re = a.re + b.re; z.im = a.im + b.im; return(z); }
static inline complex c_sub(complex a, complex b) { complex z; z.re = a.re - b.re; z.im = a.im - b.im; return(z); }
static inline complex c_mult(complex a, complex b){ complex z; z.re = a.re*b.re - a.im*b.im; z.im = a.im*b.re + a.re*b.im; return(z); }
complex rect(double mag, double ph);
static inline complex c_conj(complex z) { complex c; c.re = z.re; c.im = -(z.im); return(c); }
static inline complex c_exp(double theta) { complex z; z.re = cos(theta); z.im = sin(theta); return(z); }
static inline complex c_zero(void) { complex z; z.re=0.; z.im=0.; return(z); }
static inline complex s_mult(complex z, double i) { complex c; c.re=z.re * i; c.im=z.im * i; return(c); }
void    c_print(complex z);
static inline double c_real(complex z) { return z.re; }
static inline double c_imag(complex z) { return z.im; }
static inline complex c_make(double re, double im) {complex z; z.re = re; z.im = im; return z;}
#endif /* _COMPLEX_H */

#endif
