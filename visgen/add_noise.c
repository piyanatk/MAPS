/************************************************************************/
/*                                                                      */
/* The job of this routine is to figure out the noise level             */
/* corresponding to the current baseline, integration time, bandwidth   */
/* and frequency, generate gaussian noise samples for real and          */
/* imaginary components, and add them to the simulated visibility       */
/* measurement.  System noise is a uv plane effect, so can be added     */
/* after uv plane convolution and integration.  The routine takes the   */
/* elevation dependence of the antenna response into account by         */
/* interpolating a lookup table stored in an external file.  It also    */
/* includes the effects of weighting in the beamforming at each station */
/* (via a weight sum computed in makebeam.c)                            */
/*                                                                      */
/* Created 27 June 2002 by CJL                                          */
/*                                                                      */
/************************************************************************/

/* From synthesis imaging workshop book, chapter 9 (sensitity) 
 * equation 9-10, the general form of additive noise (in Janksy) 
 * for a visibility from baseline i,j is:
 *
 *   delta_S = 1/[eta*sqrt(2*delta_nu*tau)]*sqrt[S_c^2 + S_t^2 + 
 *             S_t*(SEFD_i + SEFD_j) + SEFD_i*SEFD_j]
 *
 * where:
 * eta = overall system efficiency (between 0 and 1)
 * delta_nu = channel bandwidth (Hz)
 * tau = integration time (sec)
 * S_c = correlated power (Jy)
 * S_t = total power received by a station (Jy)
 * SEFD= system equivalent flux desnity of a station (Jy)
 *
 * In the above expression, it is assumed that each antenna sees the 
 * same S_t, which may not be true if different types of antennas are 
 * used in the array. (e.g. for VLBI). It is also assumed that the
 * receiver temp is not the same as "antenna temp" which is often 
 * defined to include the signal from astronomical sources. Receiver 
 * temp is just the noise intrinsic to the receiver itself.
 *
 * In this function, the visibility is assumed to already contain the 
 * correlated power, in Jy. This is used directly as S_c.
 *
 * How to find the total power into the antenna (S_t)? 
 * Answer: Use Parseval's theorem.
 *
 * The total power into the station is:
 *
 *   \int_{sky} I(l,m) B(l,m) dldm
 * 
 * where I() is the sky brightness in Jy/sr and B() is the primary 
 * beam power pattern on the sky.
 * Using Parsevals' theorem, this is equal to:
 *
 *   \int_{u,v plane} I'(u,v) B'(u,v) dudv
 *
 * where I' and B' are the fourier transforms of the sky brightness 
 * and primary beam power patterns, respectively. For our pixellised
 * sky and sampled primary beams, we have almost everything we
 * need to do this integral. We just have to FFT the primary beam 
 * power, plonk this down in the middle of the u,v plane and sum 
 * over the pixels.
 *
 */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "array_specification.h"
#include "observing_param.h"
#include "type_comp.h"
#include "comp_func.h"
#include "uvgrid.h"
#include "station_beam.h"
#include "get_patch.h"
#include "add_noise.h"

/* private function prototypes */
void gauss_noise (double *value1, double *value2);


/* private global vars */
int debug_addnoise=0;


void addnoise_set_debug(int level) {
  debug_addnoise = level;
}


int seed_noise (int seed, int myid)
{
  long seedval;
  // initialization function 
  seedval = seed * (myid+1); 
  srand48(seedval);
  if(debug_addnoise) 
    printf("seed: %d, myid: %d, seed value is %ld \n ", seed, myid, seedval);
  return 0;
}


int add_noise (// st1spec, st2spec, freq, obsparam, visibility)
	       const struct stationspec *st1spec,
	       const struct stationspec *st2spec,
	       const struct stnbeam     *st1beam,
	       const struct stnbeam     *st2beam,
	       const float   bandwidth_Hz,         // in Hz
	       const float   integration_time,     // in seconds.
	       struct uvgrid_parms *uvgrid,
	       const float   global_sefd,          // in Jy
	       const float   system_efficiency,    // dimensionless, [0..1]
	       complex *visibility,
	       int n_pols,
	       float *vnoise)
{
  float st1_sefd,st2_sefd,st1_totalpower=0,st2_totalpower=0;
  float corr_power=0,sigma=0,denom=0,inside_part=0;
  double rnoise,inoise;
  //    struct stnbeam *beam1,*beam2;
  //    struct conv_fn *beamconvl;
  int     pol;
  /*    static struct uvpatch *uvcenter=NULL;  // Not used? -- benkev */
    
  /* sanity checks */
  if (bandwidth_Hz<=0.0 || integration_time<=0.0 || 
      system_efficiency<=0.0 || global_sefd<=0.0) {
    msg("add_noise: invalid input params. bw: %g, " \
	"tau: %g, eta: %g, SEFD: %g\n", 3,
	bandwidth_Hz, integration_time, system_efficiency, global_sefd);
    return 1;
  }

  /* set default values */
  st1_sefd = (st1spec->sefd > 0.0 ? st1spec->sefd : global_sefd);
  st2_sefd = (st2spec->sefd > 0.0 ? st2spec->sefd : global_sefd);
  if (debug_addnoise) {
    fprintf(stdout,"bw: %g, tau: %g, eta: %g\n",
	    bandwidth_Hz, integration_time, system_efficiency);
    fprintf(stdout,"st1sefd: %g, st2sefd: %g\n", st1_sefd, st2_sefd);
  }

  /* if we have access to the pixellised u,v plane and station beams, 
   * then use them to calculate the total received power in each station */
  /******* NOT finished yet
   * if(uvgrid !=NULL && st1beam!=NULL && st2beam!=NULL && 
   *    st1beam->valid_flag && st2beam->valid_flag) {
   * // clip out the central patch of the u,v plane. Only need to do this once
   *   if (uvcenter ==NULL) {
   *     uvcenter = calloc(1,sizeof(struct uvpatch));
   *     if (uvcenter==NULL) {
   *	   msg("add_noise: no malloc\n", 3);
   *	   return -1;
   *     }
   *   get_patch(uvgrid, 0.0, 0.0, 0.0, st1beam->size, 
   *             st1beam->size, uvcenter);
   *   }
   *   // create FFT of station power pattern     
   * }
   *************/

  /* pre-calculate the constant parts of sigma */
  denom = 1.0/(system_efficiency*sqrt(2.0*bandwidth_Hz*integration_time));
  inside_part = st1_totalpower*st2_totalpower + st1_totalpower*st1_sefd
    + st2_totalpower*st2_sefd + st1_sefd*st2_sefd;

  /* now calculate sigma for each pol product and add a random noise term */
  for (pol=0; pol < n_pols; pol++) {
    corr_power = c_mag(visibility[pol]);
        
    /* caculation of error magnitude from equation above. 
     * This has been generalised for
     *  stations that don't have the same total power */
    sigma = sqrt(corr_power*corr_power + inside_part)*denom;
        
    // Gaussian distribution of unit variance
    gauss_noise (&rnoise, &inoise);
    if (debug_addnoise) {
      fprintf(stdout,"pol %d, sigma: %g, noiseval: (%g,%g)\n",
	      pol, sigma, sigma*rnoise, sigma*inoise);
    }

    // Add it to the simulated visibility
    visibility[pol] = c_add(visibility[pol], 
			    c_make(sigma*rnoise,sigma*inoise));
  }
  // return sigma;
  *vnoise = sigma;
  return 0;
}



void
gauss_noise (/* value1, value2) */
	     double *value1,
	     double *value2)
{
  double fac, rsq, v1, v2;

  while (1)
    {
      // Random location in unit square
      v1 = 2.0 * drand48() - 1.0;
      v2 = 2.0 * drand48() - 1.0;
      rsq = v1*v1 + v2*v2;
      // Take only those inside unit circle
      if (rsq < 1.0) break;
    }
  // Following the recipe ...
  fac = sqrt (-2.0 * log (rsq) / rsq);
  *value1 = v1 * fac;
  *value2 = v2 * fac;
}
