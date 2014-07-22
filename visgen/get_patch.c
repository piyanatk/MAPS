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
// void PatchColumn (FileHandle fh, long i, long j0, long j1, long nc, 
//                   complex **pp);
void PatchColumn (FILE *fh, long i, long j0, long j1, long nc, complex **pp);
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
get_patch (/* uvgrid, u, v, freq, udim, vdim, patch) */
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
  complex *p;

  if (udim <= 0.0 || vdim <= 0.0) return 1; // sizes must be positive

  /* make space for patch if not done already */
  for (i=0; i<uvgrid->n_pol_inputs; i++) {
    if (patch->patch[i] == NULL) {
      patch->patch[i]  = calloc(MAXPATCH*MAXPATCH,sizeof(complex));
      if (patch->patch[i] == NULL) {
	msg("get_patch: no malloc",2);
      }
    }
  }

  nc = uvgrid->ncells;           // For brevity (commonly used)
  dc = uvgrid->cellsize;

  // Check patch size
  npu = ceil (udim  / dc);
  npv = ceil (vdim  / dc);
  //printf("&&&&&&& get_patch: npu = %d, npv = %d\n", (int)npu, (int)npv);
  if (npu*npv > MAXPATCH*MAXPATCH) {
    patch->usize = npu;   /* Added for better diagnostics by benkev */
    patch->vsize = npv;   /* Added for better diagnostics by benkev */
    msg ("Requested patch size %dx%d cells too big", 2, (int)npu, (int)npv);
    return (4);
  }
  // uv grids are constrained to be square 
  umax = vmax = (nc / 2) * dc;
  umin = vmin = -(nc / 2) * dc;
  // Get boundary points of patch
  u0 = u - udim/2; 
  u1 = u0 + udim; 

  if (u0 < umin || u1 > umax) {
    msg ("u out of range, u0, u1, umin, umax = %g %g %g %g", 
	 2, u0, u1, umin, umax);
    return 2; 
  }
  v0 = v - vdim/2; 
  v1 = v0 + vdim; 
  if (v0 < vmin || v1 > vmax) {
    msg ("v out of range, v0, v1, vmin, vmax = %g %g %g %g", 
	 2, v0, v1, vmin, vmax);
    return 3; 
  }
            
                     // Calculate the positions in the array
                     // For u0 and u1 with the same sign, 
                     // array positions will be consecutive. For 
                     // different signs, array position
                     // will wrap at visibilityzN (z=X or Y); 
                     // e g, N-2, N-1, 0, 1, ...
  i0 = ArrayPosition (nc, u0, dc); 
  i1 = ArrayPosition (nc, u1, dc);
  j0 = ArrayPosition (nc, v0, dc); 
  j1 = ArrayPosition (nc, v1, dc);


  /* printf("u0,u1,v0,v1 = %g,%g,%g,%g\n",u0,u1,v0,v1); */
  /* printf("i0,i1,j0,j1 = %d,%d,%d,%d\n",i0,i1,j0,j1); */

  msg ("u0,u1,v0,v1 = %g,%g,%g,%g",-2,u0,u1,v0,v1);
  msg ("i0,i1,j0,j1 = %d,%d,%d,%d",-2,i0,i1,j0,j1);

  for (k=0; k<uvgrid->n_pol_inputs; k++) {
    // Populate the patch
    p = patch->patch[k]; // ->output
    if (i1 > i0) {
      for (i=i0; i<=i1; i++)
	PatchColumn (uvgrid->fh[k], i, j0, j1, nc, &p);
      nu = i1 - i0 + 1;
    }
                   // Wrap around edge of grid, which is stored
                   // with zero at the edge
    else {
      for (i=i0; i<nc; i++)  
	PatchColumn (uvgrid->fh[k], i, j0, j1, nc, &p);
      for (i=0; i<=i1; i++)              
	PatchColumn (uvgrid->fh[k], i, j0, j1, nc, &p);
      nu = (nc - i0) + (i1 - 0 + 1);
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

  if (nv != (p - patch->patch[uvgrid->n_pol_inputs-1]) / nu) {
    msg ("Incorrect number of bytes read in get_patch()", 2);
    return (5);
  }
 

  patch->vsize = nv;
  patch->cellsize = dc;             // i.e. copied from uvgrid

  //    for (i=0;i<patch->vsize;i++){
  //	    for (j=0;j<patch->usize;j++){
  //		    printf("%.2g %.2g  ",patch->patch[j*patch->usize + i].re,
  //				    patch->patch[j*patch->usize + i].im);
  //	    }
  //	    printf("\n");
  //    }
  //    msg ("at end of get_patch\n",-2);

  return (0);
}


                            // Get a single column (v varies) of the patch
void PatchColumn (FILE *fh, long i, long j0, long j1, long nc, complex **pp) {
  off_t filePosition; 
  int nr, hdrsize;

  hdrsize = sizeof (struct uvgrid_header);
  filePosition = ((off_t)(i * nc + j0)) * sizeof(complex) + hdrsize;
  // printf("in patchcolumn: filePosition is %lld\n",filePosition);
  fseeko (fh, filePosition, SEEK_SET);
  if (j1 > j0) {
    nr = fread (*pp, sizeof(complex), j1 - j0 + 1, fh);
    *pp += nr;                // update the pointer
    // msg ("fileposition = %lld, i,nc,j0=%d,%d,%d\n", 
    //      -3, filePosition,i,nc,j0);
  }
  // read from two different parts of 
  // the array, wrapping at n
  else {
                            // Read from j0 to end of column
    nr = fread (*pp, sizeof(complex), nc - j0, fh);
    *pp += nr;              // update the pointer
                            // Wrap back to start of column & read to j1
    filePosition = ((off_t)i*nc+0) * sizeof(complex) + hdrsize;
    fseeko (fh, filePosition, SEEK_SET);
    nr = fread (*pp, sizeof (complex), j1 - 0 + 1, fh);
    *pp += nr;                   // update the pointer
  }
  return;
}         
        
long ArrayPosition (long arraySize, double value, double spacing) {
  long i;

  if (value < 0) 
    i = arraySize - ((int)(-value/spacing + 0.5));
  else 
    i = (int)(value/spacing + 0.5);
  // Trap case of small neg values
  // which should wrap back to zero in array
  if (i == arraySize) i = 0;
  return i;
}

// Calculate the (u,v) corresponding to (i,j) in the visibility array
void ExactUV (long i, long j, long nc, double dc, double *pu, double *pv)
{
               // make the ambiguous center freq positive
    if (i <= nc/2) *pu = i * dc;
    else *pu = (i - nc) * dc;
    if (j <= nc/2) *pv = j * dc;
    else *pv = (j - nc) * dc;
    return;
}
