/************************************************************************/
/*                                                                      */
/* Defines structures needed to support handling of uv patches          */
/*                                                                      */
/* Created Feb 11 2002 by CJL                                           */
/* Modified Feb 13 2002 by Peter Sherwood                               */
/*                                                                      */
/************************************************************************/
#ifndef UVGRID
#define UVGRID 1

#define MAXPATCH 50

#include <stdio.h>
#include "type_comp.h"
//#include "LargeFiles.h"

struct uvgrid_header
{                                   // Minimum essential information
    int         skygridsize;        // Cells, need not be power of 2
    int         uvgridsize;         // >= skygridsize, must be power of 2
/*     char        modelfile[256];  // Full pathname */
/*     char        imagefile[256];  // Full pathname */
};

struct uvgrid_parms
{
  int         ncells;     // Total size in cells (normally power of 2)
  double      padfactor;  // Ratio of sky array size to sky model size
  double      cellsize;   // In wavelengths
  FILE        *fh[4];     // Used for IO. Have up to 4 files, one for each pol.
  int         n_pol_inputs; // number of pol inputs. should be between 1 and 4.
  int         dummy;      // Place additional information here defining
                          // origin of this uvgrid dataset
};

struct uvpatch
{
  double      frequency;       // MHz
  double      umin,vmin;       // u,v of lower left corner (wavelengths)
  int         usize,vsize;     // u,v dimension in cells
  double      cellsize;        // In wavelengths
  complex     *patch[4];       // pointers to 2-D complex array of visibilities,
                               // up to 4 pol products
};

#endif
