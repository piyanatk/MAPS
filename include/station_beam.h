/************************************************************************/
/*                                                                      */
/* Defines structures needed to store station beam parameters           */
/*                                                                      */
/* Created 12 Feb 2002 by SSD                                           */
/* Modified 24 July 2002 SSD to complex array                           */
/* Modified 14 Oct 2002 SSD to add stngrid structure			*/
/************************************************************************/
#ifndef STATION_BEAM_H
#define STATION_BEAM_H 1

#include "type_comp.h"
#include "type_radec.h"

#define BAD_POINT -99.0

struct stnbeam {
  complex  *beam[4]; /* complex array, one for each element of 
		      * 2x2 Jones matrix */
                     /* index of beam is px,py,qx,qy */
  int       station; /* Station number */
  int 	    size;    /* Length of 1 dimension of array (i.e. total is size^2) */
  int	    xref_pixel;	/* Beam Center */
  int	    yref_pixel; 
  int	    valid_flag;	/* flag=0 if stn beam not computed	*/
                        /* flag=1 if stn beam computed	    */
};

struct beamgrid
    {
	radec  *grid;	    /* array of (ra,dec) pairs	*/
	int	size;	    /* Dimension of grid (one edge of square) */
	int	xref_pixel; /* Center of beam sampling grid	*/
	int	yref_pixel;
    };

#endif

