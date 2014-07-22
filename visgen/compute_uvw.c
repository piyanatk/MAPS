/************************************************************************/
/*                                                                      */
/* compute_uvw.c                                                        */
/*                                                                      */
/* Given information on two stations, the time, and the observation     */
/* parameters, this routine computes u, v, w in meters for the baseline */
/*                                                                      */
/*      Inputs:     st1, st2    Filled station information structs      */
/*                  gha         Radians                                 */
/*                  obsparam    Observing parameter struct              */
/*                                                                      */
/*      Output:     u_len, v_len, w_len    In meters                    */
/*                  return value    0=OK, >0 = el. limit, <0 = bad      */
/*                                                                      */
/* Created 5 Feb 2002 by CJL                                            */
/* changed sign of longitude to make it East in el calc.  rjc 2003.7.21 */
/************************************************************************/

/* A NOTE on the definition of a baseline:
 * The de-facto standard for UVFITS files follows the AIPS convention:
 * - a baseline is defined as location(ant1)-location(ant2) where the 
 *     antenna indices are such that ant2 > ant1
 * - using this definition of a baseline, a visibility (for a point 
 *     source) is V(u,v) = I*exp(-2*pi*(ul+vm))
 * - this means that if a baseline is defined as (ant2-ant1),
 *     then the exponent in the exponential must be +2pi, not -2pi
 * - The AIPS standard is consistent with the definition of a visibility 
 *     from a coherence average point of view (e.g section 14.1 of Thompson, 
 *     Moran and Swenson or Chapter 1 of "Synthesis imaging in radio 
 *     astronomy II".)
 */

#include <stdio.h>
#include <math.h>
#include "array_specification.h"
#include "observing_param.h"
#include "utility.h"

int
compute_uvw (/* st1, st2, gha, obsparam, u_len, v_len, w_len) */
	     struct stationspec *st1,
	     struct stationspec *st2,
	     double gha,
	     struct observing_param *obsparam,
	     double *u_len,
	     double *v_len,
	     double *w_len)
{
  double xb, yb, zb, cgha, sgha, cdec, sdec;
  double long1, long2, lat1, lat2, sinel1, sinel2,el1=0,el2=0;
  double radeg = 180./3.14159265;  /* radians-to-degrees conversion factor */

  // Baseline vector in meters (see note above)
  xb = st1->x_coord - st2->x_coord;
  yb = st1->y_coord - st2->y_coord;
  zb = st1->z_coord - st2->z_coord;
  msg ("xyz %lf %lf %lf. GHA: %g", -2, xb, yb, zb, gha);

  // Just to be tidy/efficient
  cgha = cos(gha); sgha = sin(gha);
  cdec = cos(obsparam->phase_cent_DEC); sdec = sin(obsparam->phase_cent_DEC);

  // Here basic coordinate transformation
  // bracketed for clarity
  *u_len = (xb * sgha) + (yb * cgha);
  *v_len = (xb * -sdec * cgha) + (yb * sdec * sgha) + (zb * cdec);
  *w_len = (xb * cdec * cgha) - (yb * cdec * sgha) + (zb * sdec);
  msg ("(u, v, w) = (%g, %g, %g) meters", -2, *u_len, *v_len, *w_len);
    
  // Compute and check elevations
  long1 = st1->longitude; long2 = st2->longitude;
  lat1 = st1->latitude; lat2 = st2->latitude;
  sinel1 = (sin(lat1) * sdec) + (cos(lat1) * cdec * cos(gha + long1));
  /*  sinel1 = (sin(lat1) * sdec) + (cos(lat1) * cdec * cos(gha - long1)); */
  el1 = asin(sinel1);
  sinel2 = (sin(lat2) * sdec) + (cos(lat2) * cdec * cos(gha + long2));
  /*  sinel2 = (sin(lat2) * sdec) + (cos(lat2) * cdec * cos(gha - long2)); */
  el2 = asin(sinel2);


  msg ("Elevations %g, %g (radians), or %g, %g (degrees)",-2, 
       el1, el2, radeg*el1, radeg*el2);
  /* Return 1 if the elevation angles is less than the station limits */
  if (el1 <st1->low_el_limit || el1 > st1->high_el_limit
      || el2 <st2->low_el_limit || el2 > st2->high_el_limit ) return (1);

  return (0);
}
