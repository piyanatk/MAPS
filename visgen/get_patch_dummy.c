// get_patch.c: Copies out a small patch of main uv grid suitable for 
// convolution and integration.

// Author:      Peter Sherwood  sherwood@computer.org   (617) 244-0836

// 2/12/02      create from VisibilityView
// 3/3/05	changed PatchColumn routine to use standard large 
//              file routines. SSD
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include <errno.h>
#include "uvgrid.h" // includes LargeFiles.h
#include "utility.h"
#include "get_patch.h"


// Local functions
long ArrayPosition (long arraySize, double value, double spacing);
void PatchColumn_dummy (FILE *fh, long i, long j0, long j1, long nc, 
			complex **pp);
void ExactUV (long i, long j, long nc, double dc, double *pu, double *pv);

/***********************************************************************/
/*                                                                     */
/* Given a u,v coordinate and frequency, this routine goes to the main */
/* uv grid and copies out a small patch suitable for convolution and   */
/* integration.  The size of the patch depends on the sizes of both    */
/* the integration region and of the convolution functions             */
/*                                                                     */
/* Inputs:  uvgrid->.  grid parameters, including fh, from             */
/*                       open_gridfile()                               */
/*                  u, v        center of patch, in wavelengths        */
/*                  freq        MHz                                    */
/*                  udim, vdim  patch size in wavelengths              */
/*                                                                     */
/*      Output:     patch       section of main uv grid, with header   */
/*                  return      0=OK                                   */
/*                              1=udim or vdim is not positive         */
/*                              2=u out of range                       */
/*                              3=v out of range                       */
/*                              4=excessive patch size requested       */
/*                              5=bad read                             */
/*                                                                     */
/***********************************************************************/

int 
get_patch_dummy (/* uvgrid, u, v, freq, udim, vdim, patch) */
	   struct uvgrid_parms *uvgrid, 
	   double u, 
	   double v, 
	   double freq, 
	   double udim, 
	   double vdim, 
	   struct uvpatch *patch) {
  int k;
  long i,  i0, i1, j0, j1, nu=0, nv=0, nc, npu, npv;
  double dc, u0, u1, v0, v1, umin, umax, vmin, vmax;
  complex *ptchk;

  nc = uvgrid->ncells;           // For brevity (commonly used)
  dc = uvgrid->cellsize;

  npu = ceil (udim  / dc); /* Patch U size in cells */
  npv = ceil (vdim  / dc); /* Patch V size in cells */

  if (udim <= 0.0 || vdim <= 0.0) {
        patch->usize = npu; /* cells, added for better diagnostics by benkev */
	patch->vsize = npv; /* cells, added for better diagnostics by benkev */
	return 1; // sizes must be positive:  1=udim or vdim is not positive 
  }

  /* make space for patch if not done already */
  for (i=0; i<uvgrid->n_pol_inputs; i++) {
    if (patch->patch[i] == NULL) {
      patch->patch[i]  = calloc(MAXPATCH*MAXPATCH,sizeof(complex));
      if (patch->patch[i] == NULL) {
	msg("get_patch: no malloc",2);
      }
    }
  }

  // Check patch size
  if (npu*npv > MAXPATCH*MAXPATCH) {
    patch->usize = npu; /* in cells, for better diagnostics by benkev */
    patch->vsize = npv; /* in cells, for better diagnostics by benkev */
    // msg ("Requested patch size %dx%d cells too big", 2, npu, npv);
    return 4; /* 4 = excessive patch size requested       */
  }
  // uv grids are constrained to be square 
  umax = vmax = (nc / 2) * dc;
  umin = vmin = -(nc / 2) * dc;
  // Get boundary points of patch
  u0 = u - udim/2; 
  u1 = u0 + udim; 
  v0 = v - vdim/2; 
  v1 = v0 + vdim; 

  if (u0 < umin || u1 > umax) {
    msg ("u out of range, u0, u1, umin, umax = %g %g %g %g", 
	 -1, u0, u1, umin, umax);
    return 2; /*  2 = u out of range */
  }
  if (v0 < vmin || v1 > vmax) {
    msg ("v out of range, v0, v1, vmin, vmax = %g %g %g %g", 
	 -1, v0, v1, vmin, vmax);
    return 3;  /* 3 = v out of range */
  }
            
                                   // Calculate the positions in the array
                                   // For u0 and u1 with the same sign, 
                                   // aray positions will be consecutive. For 
                                   // different signs, array position
                                   // will wrap at visibilityzN (z=X or Y); 
                                   // e g, N-2, N-1, 0, 1, ...
  i0 = ArrayPosition (uvgrid->ncells, u0, uvgrid->cellsize); 
  i1 = ArrayPosition (uvgrid->ncells, u1, uvgrid->cellsize);
  j0 = ArrayPosition (uvgrid->ncells, v0, uvgrid->cellsize); 
  j1 = ArrayPosition (uvgrid->ncells, v1, uvgrid->cellsize);
  msg ("u0,u1,v0,v1 = %g,%g,%g,%g",-2,u0,u1,v0,v1);
  msg ("i0,i1,j0,j1 = %d,%d,%d,%d",-2,i0,i1,j0,j1);

  for (k=0; k<uvgrid->n_pol_inputs; k++) {
    // Populate the patch
    ptchk = patch->patch[k]; // ->output
    if (i1 > i0) {
      for (i=i0; i<=i1; i++)
	PatchColumn_dummy (uvgrid->fh[k], i, j0, j1, nc, &ptchk);
      nu = i1 - i0 + 1;
    }
                   // Wrap around edge of grid, which is stored
                   // with zero at the edge
    else {
      for (i=i0; i<uvgrid->ncells; i++)  
	PatchColumn_dummy (uvgrid->fh[k], i, j0, j1, nc, &ptchk);
      for (i=0; i<=i1; i++)              
	PatchColumn_dummy (uvgrid->fh[k], i, j0, j1, nc, &ptchk);
      nu = (uvgrid->ncells - i0) + (i1 - 0 + 1);
    }
                  //    msg ("past populating the patch",-2);
  }
                  // Fill in the patch information
  patch->frequency = freq;
  ExactUV(i0, j0, nc, dc, &patch->umin, &patch->vmin);
                 // nu=# of points in the u direction
                 // Check nv from pointer change
  patch->usize = nu;
  if (j1 > j0) nv = j1 - j0 + 1;
  else nv = (nc - j0) + (j1 + 1);
    
   
                     /* this only tests the last patch now */
                     /* Assumption is that all pol fields are the same size */
                     /* A condition asserted in open_gridfile */
                     /* I had to do this as the previous logic was flawed */

  if (nv != (ptchk - patch->patch[uvgrid->n_pol_inputs-1]) / nu) {
    msg ("Incorrect number of bytes read in get_patch(): " \
	 "nu = %ld, nv = %ld, " \
	 "ptchk - patch->patch[uvgrid->n_pol_inputs-1])/nu = %ld", 2, nu, nv, 
	 (long)(ptchk - patch->patch[uvgrid->n_pol_inputs-1])/nu);
    msg (" ",2);
    return 5; /* 5 = bad read; not needed here in get_patch_dummy() */
  }
 

  patch->vsize = nv;
  patch->cellsize = dc;             // i.e. copied from uvgrid
  /* printf("nu = %ld, npu = %ld, nv = %ld, npv = %ld\n", nu, npu, nv, npv); */
  return 0;
}


//
// Imitate getting a single column (v varies) of the patch
//
void PatchColumn_dummy (FILE *fh, long i, long j0, long j1, long nc,    
		  complex **pp) 
{
  off_t filePosition;
  int nr, hdrsize;

  hdrsize = sizeof (struct uvgrid_header);
  filePosition = ((off_t)(i * nc + j0)) * sizeof(complex) + hdrsize;
  /* fseeko (fh, filePosition, SEEK_SET);     // Remove disc access  */
  if (j1 > j0) {
    /* nr = fread (*pp, sizeof(complex), j1 - j0 + 1, fh); // No file ops! */ 
    nr = j1 - j0 + 1;
    *pp += nr;                // update the pointer
  }
  // read from two different parts of
  // the array, wrapping at n
  else {
                            // Read from j0 to end of column
    /* nr = fread (*pp, sizeof(complex), nc - j0, fh); // No file ops! */  
    nr = nc - j0;  /* Added by benkev */
    *pp += nr;              // update the pointer
                            // Wrap back to start of column & read to j1
    filePosition = ((off_t)i*nc+0) * sizeof(complex) + hdrsize;
    /* fseeko (fh, filePosition, SEEK_SET); // No file ops! */ 
    /* nr = fread (*pp, sizeof (complex), j1 - 0 + 1, fh); // No file ops! */ 
    nr = j1 - 0 + 1;  /* Added by benkev */
    *pp += nr;                   // update the pointer
  }
  return;
}
        
/* long  */
/* ArrayPosition (long arraySize, double value, double spacing) */
/*     { */
/*     long i; */

/*     if (value < 0) i = arraySize - ((int)(-value/spacing + 0.5)); */
/*     else i = (int)(value/spacing + 0.5); */
/*                               // Trap case of small neg values */
/*                               // which should wrap back to zero in array */
/*     if (i == arraySize) i = 0; */
/*     return i; */
/*     } */

/* // Calculate the (u,v) corresponding to (i,j) in the visibility array */
/* void ExactUV (long i, long j, long nc, double dc, double *pu, double *pv) */
/* { */
/*                // make the ambiguous center freq positive */
/*     if (i <= nc/2) *pu = i * dc; */
/*     else *pu = (i - nc) * dc; */
/*     if (j <= nc/2) *pv = j * dc; */
/*     else *pv = (j - nc) * dc; */
/*     return; */
/* } */
