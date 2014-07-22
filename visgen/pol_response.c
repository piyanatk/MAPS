/************************************************************************
* Functions to support polarised response in the array
*
* Created July 2007 by Randall Wayth
*
************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include "utility.h"
#include "comp_func.h"
#include "pol_response.h"


/* public global variables */
complex stokes_2_receptor[16]; /* [pol1*I,pol1*Q,pol1*U,pol1*V,pol2*u etc ] */
int     pol_init=0;

/* private global variabels */
static const int m_conv[16] = {0,1,4,5,2,3,6,7,8,9,12,13,10,11,14,15};


/******************************
 create a Mueller matrix from two Jones matrices. this is the Kronecker (outer) product of the two
 matrices with the conj of the second. The Jones matrices are 4 element, 1D arrays with index px,py,qx,qy. This code
 has been adapted from the code in compute_beamconvl(). They MUST be kept in sync. 
 ******************************/
int pol_make_mueller(int n_prod, complex *j1, complex *j2, complex *m) {
  int i1,i2,m_index;

  for (i1=0; i1 < n_prod; i1++) {
    for (i2=0; i2 < n_prod; i2++) {
      m_index = m_conv[i1*4 + i2];
      /* multiply terms including complex conj of second term */
      m[m_index] = c_mult(j1[i1],c_conj(j2[i2]));
    }
  }
  return 0;
}


/******************************
 multiply the Mueller and Stokes converter matrices into a single matrix
 ****************************** */
int pol_make_MS(const int n_prod,const complex *m,const complex *s, complex *result) {
  int j,r,c,m_ind,s_ind;
  complex temp;

  for (r=0; r<n_prod; r++){
    for (c=0; c<n_prod; c++) {
      temp = c_zero();
      for (j=0; j<n_prod; j++) {
        m_ind = j+n_prod*r;
        s_ind = c+n_prod*j;
        temp = c_add(temp,c_mult(m[m_ind],s[s_ind]));
      }
      result[c+n_prod*r] = temp;
    }
  }
  return 0;
}


/*****************************
  multiply a mueller type matrix (4x4) with the input stokes vector to form MI. (Stokes vector on the right).
******************************/
int pol_make_MI(const int n_prod,const complex *m, float *flux, complex *result) {
  int i,j;
  complex temp;

  for (j=0; j<n_prod; j++) {
    temp = c_zero();
    for (i=0; i<n_prod; i++) {
      temp = c_add(temp,s_mult(m[i+n_prod*j],flux[i]));
    }
    result[j] = temp;
  }
  return 0;
}


/* place complex elements in stokes to XY converter. overall matrix is:
S = | 1 1 0  0| |I|
    | 0 0 1  i| |Q|
    | 0 0 1 -i| |U|
    | 1 -1 0 0| |V|
*/
void init_stokes_2_xy() {
  pol_init = 1;

  /* first row of matrix */
  stokes_2_receptor[0] = c_make(1.0,0.0);
  stokes_2_receptor[1] = c_make(1.0,0.0);
  stokes_2_receptor[2] = c_make(0.0,0.0);
  stokes_2_receptor[3] = c_make(0.0,0.0);

  /* second row */
  stokes_2_receptor[4] = c_make(0.0,0.0);
  stokes_2_receptor[5] = c_make(0.0,0.0);
  stokes_2_receptor[6] = c_make(1.0,0.0);
  stokes_2_receptor[7] = c_make(0.0,1.0);

  /* third row */
  stokes_2_receptor[8] = c_make(0.0,0.0);
  stokes_2_receptor[9] = c_make(0.0,0.0);
  stokes_2_receptor[10] = c_make(1.0,0.0);
  stokes_2_receptor[11] = c_make(0.0,1.0);

  /* 4th row of matrix */
  stokes_2_receptor[12] = c_make(1.0,0.0);
  stokes_2_receptor[13] = c_make(-1.0,0.0);
  stokes_2_receptor[14] = c_make(0.0,0.0);
  stokes_2_receptor[15] = c_make(0.0,0.0);

}


/* place complex elements in stokes to RL converter. overall matrix is:
S = | 1 0 0  1| |I|
    | 0 -j 1 0| |Q|
    | 0 -j -1 0| |U|
    | 1 0 0 -1| |V|
*/
void init_stokes_2_rl() {
  pol_init = 1;
  msg("init_stokes_2_rl: not implemented yet",3);
}


void init_stokes_I_only() {
  pol_init = 1;

  /* first row of matrix */
  stokes_2_receptor[0] = c_make(1.0,0.0);
  stokes_2_receptor[1] = c_make(0.0,0.0);
  stokes_2_receptor[2] = c_make(0.0,0.0);
  stokes_2_receptor[3] = c_make(0.0,0.0);

  /* second row */
  stokes_2_receptor[4] = c_make(0.0,0.0);
  stokes_2_receptor[5] = c_make(0.0,0.0);
  stokes_2_receptor[6] = c_make(0.0,0.0);
  stokes_2_receptor[7] = c_make(0.0,0.0);

  /* third row */
  stokes_2_receptor[8] = c_make(0.0,0.0);
  stokes_2_receptor[9] = c_make(0.0,0.0);
  stokes_2_receptor[10] = c_make(0.0,0.0);
  stokes_2_receptor[11] = c_make(0.0,0.0);

  /* 4th row of matrix */
  stokes_2_receptor[12] = c_make(0.0,0.0);
  stokes_2_receptor[13] = c_make(0.0,0.0);
  stokes_2_receptor[14] = c_make(0.0,0.0);
  stokes_2_receptor[15] = c_make(0.0,0.0);
}
 
/******************************
*******************************/
void print_complex_matrix(int num_pol_prod, complex *m) {
  int i,j;
  for (j=0; j<num_pol_prod; j++) {
    for (i=0; i<num_pol_prod; i++) {
      fprintf(stdout,"(%g,%g) ", c_real(m[i+4*j]),c_imag(m[i+4*j]) );
    }
    fprintf(stdout,"\n"); 
  }
}

