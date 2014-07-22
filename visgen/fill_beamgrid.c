/******************************************************************************/
/* fill_beamgrid routine:						      */
/*			                                    		      */
/* Fills in grid of (ra,dec) points at which to sample station beams.	      */
/* Each point is a radec structure of two doubles.  A 2-D rectangular grid    */
/* is mapped onto a sphere and rotated to the ra,dec of the center of the     */
/* field of view.  Station beams are evaluated at each of these grid points.  */
/*									      */
/* Inputs : beam_grid structure pointer, cellsize			      */
/*                                                                            */
/* Output : filled in beam structure                                          */
/* Created Oct 14, 2002 SSD                                                   */
/******************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "utility.h"
#include "station_beam.h"
#include "observing_param.h"

#define NUM_SAMPS 32 /* size along an edge of the sky grid.      *
		      * This should probably be a power of 2     *
                      * or it might break the fourc FFT function */

int
fill_beamgrid(
	      struct beamgrid *bgrid,
	      double cellsize,
	      struct observing_param *obsparam)
{ 
  int dx, dy, index;
  double delx,dely;
  double theta, dist, samp_spacing;
  double x,y,z,xp,yp,zp,sr,sd,cr,cd;

  /* set up 2-D grid of points over which beam will be sampled for plotting  */
  /* The idea is to map the 2-D grid onto the surface of a sphere using      */
  /* a simple alg. that takes each grid point and lays it down on the sphere */
  /* preserving its distance from the grid center and its angle from North   */
  /* this is all just to get a grid of ra and dec for sampling the beam      */
  /* Define spacing of samples on sky - must sample PADDED sky array.        */
  samp_spacing = 1.0/((NUM_SAMPS)*cellsize);

  /* allocate space for the 2D array */
  bgrid->grid = calloc(NUM_SAMPS*NUM_SAMPS,sizeof(radec));
  if (bgrid->grid ==NULL) {
    msg("ERROR: no malloc in fill_beamgrid",3);
    return -1;
  }

  /* pre-calculate sines and cosines */
  sr = sin(obsparam->phase_cent_RA);
  cr = cos(obsparam->phase_cent_RA);
  sd = sin(obsparam->phase_cent_DEC);
  cd = cos(obsparam->phase_cent_DEC);

  for (dy=0;dy<NUM_SAMPS;dy++){
    for (dx=0;dx<NUM_SAMPS;dx++){

      /* Get offsets from beam center defined as NUM_SAMPS/2 */
      /* Note that RA increases to the right (East). This is because     *
       * the grid which these coords are being used for is ultimately    *
       * Fourier transformed and applied to the UV plane. The UV plane   *
       * is a right-handed coord system, so this grid must also be one.  *
       * The easiest way to visualise is to think of looking down on     *
       * the array from afar, through the beamgrid. In this case u,v is  *
       * in the x,y direction of a normal cartesian plane. u points in   *
       * the direction of +RA (east), so RA increases to the right.      *
       * This is the opposite of the usual "view", where we are looking  *
       * out onto the sky and RA increases to the left.                  */
      delx = (dx-NUM_SAMPS/2.0)*samp_spacing;
      dely = (dy-NUM_SAMPS/2.0)*samp_spacing;
	  
      /* Get index of this beam pixel in 1-D complex array.   */
      /* Note that lower left element is (0,0)		*/
      index = dy*NUM_SAMPS + dx;

      theta = atan2(delx,dely); 
      dist = sqrt(delx*delx+dely*dely);

      /* We need to include a SINE projection in this grid. The grid runs * 
       * over conceptual FOV with equally spaced points on an angular     *
       * grid. In reality it needs to be equally spaced points on the l,m *
       * grid. We therefore need to introduce a SINE projection. However  *
       * the grid does run from -FOV/2 to FOV/2 in l,m (although it is    *
       * set up as angle) so we need to shrink the beam accordingly.      *
       * Since l = sin(theta) and the points on the grid are angle units, *
       * we have to take the asin of the grid unit so that the non        *
       * equally spaced angles on the grid make equally spaced units in   *
       * l,m. For large FOVs, this will map some of the grid points more  *
       * than pi/2 radian away, so these must be discarded. */

      /* scale the distance according to the desired SINE projection */
      if (fabs(dist) <= 1.0) {
	dist = asin(dist);

	/* Get x,y,z of grid point assuming great circle distance */
	/* from grid center(tangent point of 2-D plane)	*/
	x = cos(dist);
	y = sin(dist)*sin(theta);
	z = sin(dist)*cos(theta);
	  
	/* Now rotate x,y,z to desired ra,dec 	*/  
	xp = -y*sr + cr*(x*cd - z*sd);
	yp = sr*(x*cd - z*sd) + y*cr;
	zp = x*sd + z*cd;

	/* Fill appropriate grid point with (ra,dec) 	*/
	bgrid->grid[index].ra = atan2(yp,xp);
	bgrid->grid[index].dec = atan2(zp,sqrt(xp*xp+yp*yp));
      } else {
	/* grid point is more than pi/2 radian from the field centre. */
	bgrid->grid[index].ra = BAD_POINT;
	bgrid->grid[index].dec = BAD_POINT;
      }

      msg("fill_beamgrid dx,dy: %d,%d (ra,dec): (%.6f  %.6f), " \ 
	  "dist: %g, angle: %g",-1,
	  dx, dy, bgrid->grid[index].ra, bgrid->grid[index].dec, dist, theta);
    }
  }

  bgrid->size = NUM_SAMPS;
  bgrid->xref_pixel = NUM_SAMPS/2;
  bgrid->yref_pixel = NUM_SAMPS/2;
  return 0;
}
