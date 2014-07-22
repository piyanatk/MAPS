/************************************************************************
* Functions to support polarised response in the array
*
* Created July 2007 by Randall Wayth
*
************************************************************************/
#ifndef POL_RESPONSE_H
#define POL_RESPONSE_H


/* private function prototypes */
void init_stokes_2_xy();
void init_stokes_I_only();
int pol_make_mueller(int n_prod, complex *j1, complex *j2, complex *m);
int pol_make_MS(const int n_prod,const complex *m,const complex *s, complex *result);
int pol_make_MI(const int n_prod,const complex *m, float *flux, complex *result);
void print_complex_matrix(int num_pol_prod, complex *m);

extern complex stokes_2_receptor[16];
extern int     pol_init;

#endif
