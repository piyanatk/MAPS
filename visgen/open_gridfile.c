/************************************************************************/
/*                                                                      */
/* Opens up the main uv grid file, representing in-beam sources, for    */
/* reading of selected small patches using get_patch()                  */
/*                                                                      */
/*      Inputs:     filename        Name of grid file                   */
/*                  fov             field of view size in radians       */
/*                                                                      */
/*      Output:     uvgrid          Filled struct with essential info   */
/*                  return value    0=OK, else bad                      */
/*                                                                      */
/* Created 2/12/02 by Peter Sherwood, based on VisibilityView           */
/* Modified to use file header, and simplified, Feb 17 2002 by CJL      */
/* Removed obsolete large file handling, Mar 4, 2005 by SSD             */
/************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "uvgrid.h"
#include <sys/stat.h>
#include <errno.h>
#include "utility.h"

int 
open_gridfile (// filename, fov, uvgrid)
	       char *filename[4], 
	       double fov, 
	       struct uvgrid_parms *uvgrid)
{
  int er, hdrsize,i; 
  off_t ncomplex, nc2, data_size=0;
  struct uvgrid_header uh;
  struct stat fileInfo;

  hdrsize = sizeof (struct uvgrid_header);

  for (i = 0; i < uvgrid->n_pol_inputs; i++) {

    if (filename[i] == NULL) continue;

    //    printf ("open_gridfile is reading filename: %s\n",filename[i]);
    if (stat (filename[i], &fileInfo)) {
      er = errno;
      if (er == ENOENT) msg ("File '%s' does not exist", 2, filename[i]);
      else if (er == EACCES) msg ("No permission to access file '%s'", 2, 
				  filename[i]);
      else msg ("Unable to get information for file '%s'", 2, filename[i]);
      return (1);
    }
    // Get the size of the matrix represented 
    // by the file, and do sanity check
    msg("size of uv gridfile %s is %lld bytes",1,filename[i],
	(long long) fileInfo.st_size);
    ncomplex = (fileInfo.st_size - hdrsize) / sizeof (complex);
    if (ncomplex * sizeof (complex) != (fileInfo.st_size - hdrsize)) {
      msg ("File '%s' size not a multiple of size of complex (%d)", 2, 
	   filename[i], sizeof(complex));
      return (2);
    }
   
    // Open the visibility file, and read header
    uvgrid->fh[i] = fopen(filename[i], "r");
    if (uvgrid->fh[i] == NULL) {
      msg("Failed to open file %s",3,filename[i]);
      return 1;
    }

    if (fread (&uh, hdrsize, 1, uvgrid->fh[i]) != 1) {
      msg ("Failed to read uvgrid header from '%s'", 2);
      return (3);
    }
    msg("uv grid size is %d",0,uh.uvgridsize);
    // Another file size check
    nc2 = uh.uvgridsize * uh.uvgridsize;
    if (nc2 != ncomplex) {
      msg ("Bad file size in open_gridfile (grid=%dx%d, ncomplex=%lld)", 2,
	   uh.uvgridsize, uh.uvgridsize, ncomplex);
      return (4);
    }
    // paranoia check all other polarisation 
    // files are the same size
    if ( i == 0) {		    
      data_size = nc2;
    }
    else if (data_size != nc2) {
      msg("Error: File '%s' data size (%lld complex floats) may not be " \
	  "compatible with data size of '%s' (%lld complex floats)", 2, 
	  filename[i], nc2,filename[0],data_size);
      msg("Error: bad size",2);
      return (4);  
    
    }


    
  }
  // Fill in the remaining uvgrid parameters
  uvgrid->ncells = uh.uvgridsize; 
  uvgrid->padfactor = (double)uh.uvgridsize / (double)uh.skygridsize;
  uvgrid->cellsize = 1.0 / (fov * uvgrid->padfactor);
  msg ("Uvgrid file '%s', (ncells, pad, spacing) = (%d, %g, %g).", 2,
       filename[0], uvgrid->ncells, uvgrid->padfactor, uvgrid->cellsize);
  msg ("Max wavenumber supported by this image: %g (wavelengths)", 2, 
       uvgrid->ncells*uvgrid->cellsize/(2*uvgrid->padfactor));
            
  return 0;
}

