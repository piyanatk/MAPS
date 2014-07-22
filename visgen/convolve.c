/************************************************************************/
/*                                                                      */
/* This routine performs a convolution, with the target array and       */
/* the convolving function defined on the same grid size.  It is the    */
/* responsibility of the routine which creates the convolving function  */
/* to ensure that it is appropriately normalized - this routine does    */
/* not perform normalization.                                           */
/*                                                                      */
/*      Inputs:     target      function to be convolved                */
/*                  ntx, nty    Size of target array in cells           */
/*                  cfunc       Convolving function                     */
/*                                                                      */
/*      Output:     target      Convolved result                        */
/*                  return value   0=OK, else bad                       */
/*                                                                      */
/* Created 18 July 2002 by CJL                                          */
/* Modified 3 Dec 2002 by SSD to read in patch by columns, not rows.    */
/* Modified 18 Jul 06 by SSD: make 'result' static external var         */
/*                                                                      */
/************************************************************************/
#include <string.h>
#include <stdlib.h>
#include "utility.h"
#include "comp_func.h"
#include "uvgrid.h"
#include "convolve.h"
#include "type_comp.h"
#include "pol_response.h"


/*********************************
**********************************/
int
do_convl (// patch, beamconvl, ionoconvl)
	  struct uvpatch *input,
	  struct uvpatch *result,
	  struct conv_fn *beamconvl,
	  int n_pols) {
  int ret;
  
  //    msg ("before beam_convl\n",-2);
  //   msg ("patch->usize=%d   patch->vsize=%d\n", -2, 
  //        patch->usize, patch->vsize);
  //   msg ("beamconvl_xsize = %d   beamconvl_ysize = %d\n",-2, 
  //        beamconvl->xsize,beamconvl->ysize);
  ret = convolve (input, result, beamconvl,n_pols);
  if (ret != 0) {
    msg ("Problem convolving with station beams", 2);
  }
  //    msg ("after beam_convl\n",-2);
  return (ret); 
}


/************************************
 ************************************/
int
convolve(
	 struct uvpatch *input,  /* Patch from get_patch(). contains an array *
				  * of 4 pointers to patch cutouts from the   *
				  * Stokes input UV grids. The pointers will  *
				  * be NULL if there is no input for the      *
				  * corresponding Stokes parameter.           */

	 struct uvpatch *result, /* result of convolution goes here. There    *
				  * may be more pol products than inputs      */
	 struct conv_fn *cfunc,  /* baseline beam from compute_beamconvol().  *
				  * Contains Mueller matrix elements          */
	 int n_pols)             /* number of output pol products (1 or 4)    *
				  * e.g. XX,XY,YX,YY=4 etc                    */
{
  int i,j,k,l, nx, ny, tx, ty, cx, cy, x, y, xcoord, ycoord, rindex;
  int ncx, ncy, xref, yref,ntx,nty,s_index=0,p_index=0;
  complex prod;

  /* init the Stokes to XY matrix */
  if (!pol_init) {
    msg("polarisation config not initialised",3);
    return 1;
  }

                                 // Copy to local variables for neatness
    ncx = cfunc->xsize;
    ncy = cfunc->ysize;
    xref = cfunc->xref_pixel;
    yref = cfunc->yref_pixel;
    ntx = input->usize;
    nty = input->vsize;

                                 // Check for valid inputs
    if ((ntx > MAXPATCH) || (nty > MAXPATCH)) {
      msg ("Maximum array size exceeded in convolve()", 2);
      return (1);
    }
    if ((ncx > ntx) || (ncy > nty)) {
      msg ("Array sizes incompatible in convolve()", 2);
      return (1);
    }

    /* make space for the result if it hasn't been made already */
    for (i=0; i<n_pols; i++) {
      if (result->patch[i]==NULL) {
	// arrays are zeroed below
        result->patch[i] = malloc(MAXPATCH*MAXPATCH*sizeof(complex)); 
        if (result->patch[i]==NULL) {
            msg("no malloc in convolve()",3);
            return -1;
        }
      }
    }

    /* the following nested loops perform the equivalent of the M*S*I *
     * matrix mult where I is a 4-element input Stokes vector, S is   *
the 4x4 array that converts Stokes to crossed linear pol and M is the Mueller matrix
       response of a baseline to the crossed linear pols. The resulting pol products will be XX,XY,YX,YY (or equiv) */

    /* pol product loop */
    for (k=0; k<n_pols; k++) {

      /* clear any previous results */
      memset(result->patch[k],0,sizeof(complex)*MAXPATCH*MAXPATCH);

      /* Stokes to receptor converter loop */
      for (l=0; l<n_pols; l++){

        p_index = k*4+l; /* index to final pol product Mueller matrix */

        /* input Stokes loop */
        for (j=0; j<n_pols; j++) {

          s_index = l*4+j; /* index to Stokes-to-receptor converter matrix */

          /* skip to the next input pol if there is no data for this one
           * or if the Stokes input does not contribute to this product */
          if (input->patch[j]==NULL || (c_real(stokes_2_receptor[s_index])==0.0 && c_imag(stokes_2_receptor[s_index])==0.0)) continue;

          /* Target array double loop - number of subarrays to sum over */
          nx = ntx - ncx + 1;
          ny = nty - ncy + 1;
          for (tx=0; tx<nx; tx++) {
            for (ty=0; ty<ny; ty++) {
              /* Target array coordinate for result value */
              xcoord = tx + xref;
              ycoord = ty + yref;
              rindex = xcoord * nty + ycoord;
              /* Convolution loop */
              prod = c_zero();
              for (cx=0; cx<ncx; cx++) {
                register complex z1,z2;
                for (cy=0; cy<ncy; cy++)  {
                  x = xcoord - xref + cx;
                  y = ycoord - yref + cy;
                  /* printf ("tindex, value = %d, %g\n", 
		   *         tindex, target[tindex].re); */
                  /* printf ("cindex, value = %d, %g\n", 
		   * cindex, cfunc->convl[cindex].re); */
                  /* this code is executed zillions of times, 
		   * so explicity write out the
                   * complex math to avoid function calls */
                  z1 = input->patch[j][x * nty + y];
                  z2 = cfunc->convl[p_index][cy * ncx + cx];
                  /* prod = c_mult (target[tindex], cfunc->convl[cindex]); */
                    /* does not include conj
                  prod.re += z1.re*z2.re - z2.im*z1.im;
                  prod.im += z2.re*z1.im + z1.re*z2.im;
                  */
                  /* includes conj */
#ifdef _COMPLEX_H
                    /* use in-built complex type */
                  prod += z1*conj(z2);
#else
                    /* mult and add explicitly written out 
		     * with old complex definition */
                  prod.re += z1.re*z2.re + z2.im*z1.im;
                  prod.im += z2.re*z1.im - z1.re*z2.im;
#endif
                  //printf("cx,cy: %d,%d, z1: %g,%g z2: %g,%g, prod: %g,%g\n",
		  //       cx,cy,c_real(z1),c_imag(z1), c_real(z2),
		  //       c_imag(z2),c_real(prod),c_imag(prod));
                }
              }
              prod = c_mult(prod,stokes_2_receptor[s_index]);
              result->patch[k][rindex] = c_add(prod,result->patch[k][rindex]);

            }
          }
        }
      }
    }
    return (0);
}


