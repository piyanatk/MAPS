/* public function prototype */

#ifndef GET_PATCH_H
#define GET_PATCH_H

int 
get_patch (/* uvgrid, u, v, freq, udim, vdim, patch) */
	   struct uvgrid_parms *uvgrid, 
	   double u, 
	   double v, 
	   double freq, 
	   double udim, 
	   double vdim, 
	   struct uvpatch *patch);

/* Added by L. Benkevitch, 01-Dec-2010 */
int 
get_patch_dummy (/* uvgrid, u, v, freq, udim, vdim, patch) */
		 struct uvgrid_parms *uvgrid, 
		 double u, 
		 double v, 
		 double freq, 
		 double udim, 
		 double vdim, 
		 struct uvpatch *patch);

#endif
