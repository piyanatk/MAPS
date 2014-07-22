/************************************************************************/
/*                                                                      */
/* Interpolates gridded uv data to a specific point in the uv plane.    */
/*                                                                      */
/*      Inputs:     patch       grid of uv data (complex type)          */
/*                  u,v         coordinates (wavelengths) of point      */
/*                              in uv plane at which to interpolate     */
/*                                                                      */
/*      Output:     interp_vis() returns complex visibility             */
/*                                                                      */
/* Created 2002.2.20 by rjc                                             */
/*         2003.9.30 - modified 16 point algorithm to use 25 points ??? */
/*                                                                      */
/* Performance: this routine uses 90 floating point multiplies,         */
/*              12 divides, and 62 adds, for a total of 164 flops.      */
/************************************************************************/
#include <stdlib.h>
#include <math.h>
#include "uvgrid.h"
#include "type_comp.h"
#include "comp_func.h"
#include "utility.h"
#include "cspline.h"
                                    // macro to fetch ij-th element of the
                                    // desired patch cell.
/* the macro includes the current pol product number */
#define gridre(i,j) c_real((patch->patch[k])[(i)*patch->vsize+j])
#define gridim(i,j) c_imag((patch->patch[k])[(i)*patch->vsize+j])
#define grid(i,j)   (patch->patch[k])[(i)*patch->vsize+j]

#define SIZ_SPLINE 4

/* private function prototypes */
static double sinc(double x);
static double qsinc(double x);
static int debug=0;



/***********************************************************/
/* Cublic spline interp, ACTUALLY called from integrate.c  */
/***********************************************************/
int interp_vis(struct uvpatch *patch, double u, double v, complex result[4]) {
  int m,n,i,j,k=0;
  double p,q;
  float sum_re=0,sum_im=0;
  static float xv[SIZ_SPLINE] = {-1., 0., 1., 2.}; /* Vertical ruler,   V */
  static float xu[SIZ_SPLINE] = {-1., 0., 1., 2.}; /* Horizontal ruler, U */
  static float *z[SIZ_SPLINE], *p2z[SIZ_SPLINE];
  static int init=0;

  // calculate distances to nearest grid points. p & q are in the
  // range 0. to 1.0, and specify where in the central cell the 
  // point (u, v) is located
  p = (u - patch->umin) / patch->cellsize;
  q = (v - patch->vmin) / patch->cellsize;

  /* nearest integer cell numbers such that offset is positive */
  m = (int)p;   /* in U direction */
  n = (int)q;   /* in V direction */
  p -= m;	 /* p,q become offsets to the nearest cell */
  q -= n;        /* in the patch coordinates */

  /* printf("interp_vis: p,q,m,n: %g,%g,%d,%d \n",p,q,m,n); */

  if (debug) {
  printf("interp_vis: p,q,m,n: %g,%g,%d,%d \n",p,q,m,n);
  for(j=0; j<SIZ_SPLINE; j++) {
    for(i=0; i<SIZ_SPLINE; i++) {
      printf("(%g,%g) ", gridre(m+i-1,n+j-1), gridim(m+i-1,n+j-1));
    }
    printf("\n");
  }
  }

  /* initialze the x and y arrays for cubic spline interp */
  if (!init) {
    for (i=0; i<SIZ_SPLINE; i++){
      z[i]  = calloc(SIZ_SPLINE,sizeof(float));
      p2z[i] = calloc(SIZ_SPLINE,sizeof(float));
    }
    init=1;
  }

  /* Removed. Instead, the rulers xu[4] and xv[4] are declared as static */
  /* L. Benkevitch */ 
  /* for (i=0; i<SIZ_SPLINE; i++){ */
  /*   xv[i]=i-1.0; /\* coords of columns (V) from UV plane *\/ */
  /*   xu[i]=i-1.0; /\* coords of rows (U) from UV plane *\/ */
  /* } */

  /* loop over pols */
  for (k = 0; k < 4; k++) {
    if (patch->patch[k] == NULL) continue;

    sum_re = sum_im = 0.0;

    /* set y values for spline (REAL part) */
    for (j=0; j<SIZ_SPLINE; j++){         /* loop over rows */
      for (i=0; i<SIZ_SPLINE; i++){       /* loop over cols */
        z[j][i] = gridre(m+i-1,n+j-1);
      }
    }

    /* initialise spline coeffs; returns p2z[4,4] array of second derivatives */
    splie2(xv,xu,z,SIZ_SPLINE,SIZ_SPLINE,p2z);
    /* now do the interp (REAL part) */
    splin2(xv,xu,z,p2z,SIZ_SPLINE,SIZ_SPLINE,q,p,&sum_re);

    /* set z values for spline (IMAG part) */
    for (j=0; j<SIZ_SPLINE; j++){
      for (i=0; i<SIZ_SPLINE; i++){
        z[j][i] = gridim(m+i-1,n+j-1);
      }
    }

    /* initialise spline coeffs; returns p2z[4,4] array of second derivatives */
    splie2(xv,xu,z,SIZ_SPLINE,SIZ_SPLINE,p2z);
    /* now do the interp (IMAG part) */
    splin2(xv,xu,z,p2z,SIZ_SPLINE,SIZ_SPLINE,q,p,&sum_im);
    result[k] = c_make(sum_re,sum_im);
  }
  return 0;
}




/************************
*************************/
int interp_vis_sinc(struct uvpatch *patch, double u, double v, 
		    complex result[4]) {
  int m,n,siz_window = 0,i,j,k = 0;
  double p,q,weight_u,weight_v,offset;
  complex sum = c_zero();

  // calculate distances to nearest grid points. p & q are in the
  // range 0. to 1.0, and specify where in the central cell the 
  // point (u,v) is located
  p = (u - patch->umin) / patch->cellsize;
  q = (v - patch->vmin) / patch->cellsize;

  m = p + 0.5;         /* nearest integer cell number m for u*/
  n = q + 0.5;         // n for v.
  p -= m;	       /* p,q become offsets to the nearest cell */
  q -= n;

  /* create a symmetric window 2*siz_window+1 cells wide which is the 
   * largest set of samples we can use that are symmetric around the 
   * central point */
  siz_window=(3);

  if (debug) {
    printf("interp_vis: p,q,m,n: %g,%g,%d,%d \n",p,q,m,n);
    for(j=-siz_window; j<=siz_window; j++) {
      for(i=-siz_window; i<=siz_window; i++) {
    	printf("(%g,%g) ",gridre(m+i,n+j),gridim(m+i,n+j));
      }
      printf("\n");
    }
    printf("interp:\n");
  }

  /* loop over pols */
  for (k=0; k< 4; k++) {
    complex interp_v[40];   // temporary space for result of interp over v
    if (patch->patch[k] == NULL) continue;

    /* interp over v. */
    for(i=-siz_window; i<= siz_window; i++) {  // loop over columns of v
      sum=c_zero();
      for(j=-siz_window; j<= siz_window; j++) { 
	// create interp point for a single column of v (fixed u)
        offset = q - j;
        weight_v = sinc(offset);
        sum = c_add(sum,s_mult(grid(m+i,n+j),weight_v));
      }
      interp_v[i+siz_window] = sum;
      if (debug) {
        printf("(%g,%g) ",c_real(interp_v[i+siz_window]),
	       c_imag(interp_v[i+siz_window]));
      }
    }
    if (debug) printf("\n");

    // we now have an interp over v, do a final, single interp over u
    sum=c_zero();
    for(i=-siz_window; i<= siz_window; i++) {
      offset = p-i;
      weight_u = sinc(offset);
      sum = c_add(sum,s_mult(interp_v[i+siz_window],weight_u));
    }
    result[k] = sum;
  }
  return 0;
}


/******************************
 ******************************/
/* old, 5 point Lagrange interp which fails badly for data that is not significantly oversampled */
int interp_vis_old(struct uvpatch *patch,
		double u,
		double v,
		complex result[4])
{
  int m,n,k;
        
    double p, pm2, pm1, pp1, pp2,
           q, qm2, qm1, qp1, qp2,
           gm2_p, gm1_p, g_p, gp1_p, gp2_p,
           gm2_q, gm1_q, g_q, gp1_q, gp2_q;


                                    // calculate distances to nearest
                                    // grid points. p & q are in the
                                    // range -0.5 to 0.5, and specify where in
                                    // the central cell the point (u, v)
                                    // is located
    p = (u - patch->umin) / patch->cellsize;
    q = (v - patch->vmin) / patch->cellsize;

    m = p + 0.5; /* nearest integer cell number */
    n = q + 0.5;
    p -= m;	 /* p,q become offsets to the nearest cell */
    q -= n;

    pm2 = p - 2.;
    pm1 = p - 1.;
    pp1 = p + 1.;
    pp2 = p + 2.;
    
    qm2 = q - 2.;
    qm1 = q - 1.;
    qp1 = q + 1.;
    qp2 = q + 2.;

    gm2_p =  pm2 * pm1 * p * pp1       / 24.;
    gm1_p = -pm2 * pm1 * p *       pp2 / 6.;
    g_p   =  pm2 * pm1     * pp1 * pp2 / 4.;
    gp1_p = -pm2       * p * pp1 * pp2 / 6.;
    gp2_p =        pm1 * p * pp1 * pp2 / 24.;

    gm2_q =  qm2 * qm1 * q * qp1       / 24.;
    gm1_q = -qm2 * qm1 * q *       qp2 / 6.;
    g_q   =  qm2 * qm1     * qp1 * qp2 / 4.;
    gp1_q = -qm2       * q * qp1 * qp2 / 6.;
    gp2_q =        qm1 * q * qp1 * qp2 / 24.;

                                    // tombstone on (u,v) too close to edge
                                    // after program is debugged, this test
                                    // could be taken out for efficiency's sake
    if (m<2 || n<2 || m>patch->usize-3 ||  n>patch->vsize-3){
      msg ("patch boundary overstepped in interp_vis, quitting. m %d n %d", 2, m, n);
      msg ("(u, v) coordinates (%lg, %lg)", 2, u, v);
      exit (1);
    }

    result[0] = result[1] = result[2] = result[3] = c_zero();

    for (k=0; k< 4; k++) {
      if (patch->patch[k] == NULL) continue;

      /* test whether interp actually does more harm than good. */
      /*      result[k].re = gridre(m, n);
      result[k].im = gridim(m, n);
      continue;*/
        
      result[k] = c_make( (gm2_p * (gm2_q * gridre (m-2, n-2) + gm1_q * gridre (m-2, n-1) + g_q   * gridre (m-2, n)    
				+ gp1_q * gridre (m-2, n+1) + gp2_q * gridre (m-2, n+2))

		       + gm1_p * (gm2_q * gridre (m-1, n-2) + gm1_q * gridre (m-1, n-1) + g_q   * gridre (m-1, n)    
				  + gp1_q * gridre (m-1, n+1) + gp2_q * gridre (m-1, n+2))
                    
		       + g_p   * (gm2_q * gridre (m, n-2) + gm1_q * gridre (m, n-1) + g_q   * gridre (m, n)    
				  + gp1_q * gridre (m, n+1) + gp2_q * gridre (m, n+2))

		       + gp1_p * (gm2_q * gridre (m+1, n-2) + gm1_q * gridre (m+1, n-1) + g_q   * gridre (m+1, n)   
				  + gp1_q * gridre (m+1, n+1) + gp2_q * gridre (m+1, n+2))
		       
		       + gp2_p * (gm2_q * gridre (m+2, n-2) + gm1_q * gridre (m+2, n-1) + g_q   * gridre (m+2, n)    
				  + gp1_q * gridre (m+2, n+1) + gp2_q * gridre (m+2, n+2))
		       )
		       ,( gm2_p * (gm2_q * gridim (m-2, n-2) + gm1_q * gridim (m-2, n-1) + g_q   * gridim (m-2, n)    
				+ gp1_q * gridim (m-2, n+1) + gp2_q * gridim (m-2, n+2))

		       + gm1_p * (gm2_q * gridim (m-1, n-2) + gm1_q * gridim (m-1, n-1) + g_q   * gridim (m-1, n)    
				  + gp1_q * gridim (m-1, n+1) + gp2_q * gridim (m-1, n+2))
	
		       + g_p   * (gm2_q * gridim (m, n-2) + gm1_q * gridim (m, n-1) + g_q   * gridim (m, n)    
				  + gp1_q * gridim (m, n+1) + gp2_q * gridim (m, n+2))

		       + gp1_p * (gm2_q * gridim (m+1, n-2) + gm1_q * gridim (m+1, n-1) + g_q   * gridim (m+1, n)   
				  + gp1_q * gridim (m+1, n+1) + gp2_q * gridim (m+1, n+2))
	
		       + gp2_p * (gm2_q * gridim (m+2, n-2) + gm1_q * gridim (m+2, n-1) + g_q   * gridim (m+2, n)    
				  + gp1_q * gridim (m+2, n+1) + gp2_q * gridim (m+2, n+2))
		       ));
      }

                                    // return the interpolated complex value
    return 0;
}


/***********************
************************/
static double sinc(double x) {
  if (x==0) return 1.0;
  return sin(M_PI*x)/(M_PI*x);
}


/**********************
quick sinc based on lookup table
***********************/
static double qsinc(double x) {
  static int init=0,siz_tab=1000;
  static float  max_x=10, *table=NULL;
  int i;

  if (!init) {
    table = calloc(siz_tab,sizeof(float));
    for(i=0; i<siz_tab; i++) {
      table[i] = sinc((i*max_x)/(double)siz_tab);
    }
    init=1;
  }
  /* table lookup */
  i=(fabs(x)/max_x)*siz_tab;
  if(i >= siz_tab) return sinc(x);
  else return table[i];
}


/* These comments have been removed from the header because the info 
 * is inconsistent with the real 4x4 spline inerpolation */
/* To do so, it utilizes a Langrange 25 point interpolation scheme      */
/* on a 5x5 array of complex visibility data.                           */
