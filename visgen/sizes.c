/************************************************************************/
/*                                                                      */
/* Utility routines to compute the sizes of various structures which    */
/* have variable length arrays in them.  Sizes are used for memory      */
/* allocation and pointer arithmetic in the output routines             */
/*                                                                      */
/* Created 18 Jan 2002 by CJL                                           */
/*                                                                      */
/************************************************************************/
#include "visibility.h"
#include "sizes.h"
#include "type_comp.h"

int
size_of_tbh (/* nchan) */
int nchan)
{
    return sizeof (struct time_block_header) + ((nchan - 1) * sizeof (struct freq_span));
}


int
size_of_vg (int n_pol)
{
    return sizeof (struct visgroup);
}

int
size_of_vf (/* nchan, stokes) */
int nchan, int n_pol)
{
    int size;

    size = sizeof (struct visfreq) + nchan*size_of_vg(n_pol) - size_of_vg(1);

    return (size);
}


int
size_of_vb (/* nfreq, maxchan, stokes) */
int nfreq,
int maxchan,
int npol)
{
                                        /* Visblock struct has enough space for */
                                        /* nfreq visfreq structs, each one having */
                                        /* maxchan correlator channels in it. */
                                        /* Some may have less than that, in which */
                                        /* case some of the space goes unused */
    return sizeof (struct visblock) + nfreq*size_of_vf(maxchan,npol) - size_of_vf(1,1);
}
