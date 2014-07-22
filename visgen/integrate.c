/************************************************************************/
/*                                                                      */
/* Given a convolved uv patch and the necessary geometrical information */
/* this routine will numerically integrate the data to generate a       */
/* complex visibility value                                             */
/*                                                                      */
/*      Inputs:     patch       grid of data to be integrated           */
/*                  stn1, stn2  Station info structs for this baseline  */
/*                  time        secs since Jan 1.0 2000                 */
/*                  freq        MHz                                     */
/*                  pobs        ptr to observation struct               */
/*                                                                      */
/*      Output:     visibility struct                                   */
/*                  return value    0=OK, else bad                      */
/*                                                                      */
/* Created 2002.2.15 by rjc                                             */
/* RBW: This function appears to integrate via the extended Simpson's rule,
    which breaks the integral up into an even number of chunks.
    The code has been updated to work with no integration at all.
    Checks for even number of cells are made in read_obs_spec. */
/************************************************************************/
#define TCELL_MAX 65
#define FCELL_MAX 65

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <errno.h>
#include <math.h>
#include "visibility.h"
#include "uvgrid.h"
#include "array_specification.h"
#include "observing_param.h"
#include "comp_func.h"
#include "utility.h"
#include "interp_vis.h"

static int debug=0;

int integrate (struct uvpatch *patch,
               struct stationspec *stn1,
               struct stationspec *stn2,
               double time,
               double freq,
    	       int    scan_index,
    	       int    integration_index,
               struct observing_param *pobs,
               complex pvis[NPOL],
               double *umean,
               double *vmean)
{
  int m, n;

    static int first = 1;
    static double weight[TCELL_MAX][FCELL_MAX];

    double f=0, f0, f1,
           t=0, t0, t1,
           tcell_size, fcell_size,
           gha,
           u0, v0, w0,
           u=0, v=0, sumwt=0,normaliser=1.0;
    complex z[4], sum[4];

    // precalculate weights first time through. These weights are for the 
    // extended Simpson's rule, which breaks the integral up into an even 
    // number of chunks (N),  hence requires the N+1 abscissa values and 
    // weights. N can be zero here, which is just no integration at all.
    if (first) {
                                    // tombstone if array not large enough
      if (pobs->nt_cells > (TCELL_MAX -1) || pobs->nf_cells > (FCELL_MAX -1)) {
        msg ("weight array bound exceeded in integrate(); increase it", 2);
        exit (1);
      }

      if(debug) printf("integrate: setting weights. nt_cells=%d, nf_cells=%d\n",
		       pobs -> nt_cells,pobs -> nf_cells);	
      for (m=0; m <= pobs -> nt_cells; m++)
        for (n=0; n <= pobs -> nf_cells; n++)
          weight[m][n] = ((m==0 || m==pobs->nt_cells) ? 1 : 2 * (m % 2 + 1))
            * ((n==0 || n==pobs->nf_cells) ? 1 : 2 * (n % 2 + 1));
      first = 0;
    }

                                    // calculate useful quantities
                                    // time limits
    t0 = time - pobs->integ_time/2.;
    t1 = t0 +  pobs->integ_time;
    tcell_size =  pobs->integ_time / (pobs->nt_cells==0?1:pobs->nt_cells);

                                    // frequency limits pcal_prob.log
    f0 = freq - pobs -> corr_chan_bw / 2.;
    f1 = f0 +  pobs -> corr_chan_bw;
    fcell_size =  pobs -> corr_chan_bw / 
      (pobs -> nf_cells==0? 1:pobs -> nf_cells) ;


    sum[0] = sum[1] = sum[2] = sum[3] = c_zero();
    z[0] = z[1] = z[2] = z[3] = c_zero();
                                    // cell loop over time
    *umean = *vmean = 0.0;
    sumwt = 0.0;
    for (m=0; m < pobs->nt_cells + 1; m++){

      // calculate time at this edge of cell
      t = t0 +  tcell_size * m;
                                  
      // get Greenwich HA for this time
      if (pobs->scan[scan_index].gha_used) {
        gha = pobs->scan[scan_index].gha_start + 
	  (integration_index+m/(pobs->nt_cells == 0? 
				1.0:(float)pobs->nt_cells))*
	      pobs->integ_time/3600.0*M_PI/12.0*(24.0/23.9344696);
      } else {
        compute_gha(pobs, t, &gha);
      }

      // calculate u(t) and v(t)
      compute_uvw (stn1, stn2, gha, pobs, &u0, &v0, &w0);
      if (debug) printf("integrate: gha = %g, uvw = %g,%g,%g\n",gha,u0,v0,w0);
                                    // cell loop over frequency
      for (n=0; n < pobs->nf_cells + 1; n++) {
        register double wt;
                      // calculate frequency at this edge of cell
        f = f0 + fcell_size * n;
            // change u and v units (meters -> wavelengths)
        u = u0 * f * 1.0e6 / VLIGHT;
        v = v0 * f * 1.0e6 / VLIGHT;
            
            // interpolate complex visibility z at (u, v)
	/* printf("&&&&&&& integrate: u = %g, v = %g\n"); */
        interp_vis(patch, u, v, z);

        // add weighted contribution into integral
        wt = weight[m][n];
#ifdef _COMPLEX_H
        sum[0] += z[0]*wt;
        sum[1] += z[1]*wt;
        sum[2] += z[2]*wt;
        sum[3] += z[3]*wt;
#else
    /* old definition of complex */
        sum[0].re += z[0].re * wt;
        sum[0].im += z[0].im * wt;
        sum[1].re += z[1].re * wt;
        sum[1].im += z[1].im * wt;
        sum[2].re += z[2].re * wt;
        sum[2].im += z[2].im * wt;
        sum[3].re += z[3].re * wt;
        sum[3].im += z[3].im * wt;
#endif

        // Compute mean uv coordinate
        *umean += u * wt;
        *vmean += v * wt;
        sumwt += wt;

        if (debug) {
          printf("integrate: wt: %g, vis. (%g,%g),(%g,%g),(%g,%g),(%g,%g)\n",
		 wt,
                 c_real(z[0]), c_imag(z[0]),c_real(z[1]), c_imag(z[1]),
                 c_real(z[2]), c_imag(z[2]),c_real(z[3]), c_imag(z[3]));
        }
      }
    }

                                    // normalize to find average per cell
    if (pobs->nt_cells > 0) normaliser *= 1./(3.0*pobs -> nt_cells);
    if (pobs->nf_cells > 0) normaliser *= 1./(3.0*pobs -> nf_cells);
    // this does not work for N=0:
    //normaliser = 1. / (9. * pobs -> nt_cells *  pobs -> nf_cells);  
    sum[0] = s_mult (sum[0], normaliser);
    sum[1] = s_mult (sum[1], normaliser);
    sum[2] = s_mult (sum[2], normaliser);
    sum[3] = s_mult (sum[3], normaliser);

                                    // return the integrated visibility
    pvis[0] = sum[0];
    pvis[1] = sum[1];
    pvis[2] = sum[2];
    pvis[3] = sum[3];
    /*    msg("f %g t %7.4f u %g v %g integ. vis. (%g,%g),(%g,%g),(%g,%g)," \
          "(%g,%g)", 
	  0, f, t, u, v,
          sum[0].re, sum[0].im, sum[1].re, sum[1].im, sum[2].re, sum[2].im, 
	  sum[3].re, sum[3].im);*/
    msg("f %g t %7.4f u %g v %g ", -1, f, t, u, v);

    *umean /= sumwt;
    *vmean /= sumwt;

    return (0);
}
