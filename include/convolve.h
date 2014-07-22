/************************************************************************/
/*                                                                      */
/* Defines structures needed to store convolving function information   */
/*                                                                      */
/* Created 18 june 2002  by SSD                                         */
/* Changed to complex numbers, 19 July 2002 by CJL                      */
/* Deal with special needs of ionosphere, 22 July 2002 by CJL           */
/*                                                                      */
/************************************************************************/
#ifndef CONVOLVE_H
#define CONVOLVE_H 1

#include "type_comp.h"
#include "uvgrid.h"

struct conv_fn
{
    complex         *convl[16];     /* Complex convolution function, one for 
				     * each Mueller matrix element */
                                    /* Stored by row, lower left first */
    int             xsize;          /* Size of function */
    int             ysize;          /* xsize*ysize  */
    int             xref_pixel;     /* Center of convolving fn. */
    int             yref_pixel;
};


/* public function prototype */
int convolve (struct uvpatch *input, struct uvpatch *result, 
	      struct conv_fn *cfunc, int n_pols);
int do_convl (struct uvpatch *input, struct uvpatch *result,
	      struct conv_fn *beamconvl, int n_pol);
#endif
