/************************************************************************/
/*                                                                      */
/* patch_size.c                                                         */
/*                                                                      */
/* Figures out what patch size in uv cells is needed to support         */
/* convolution and integration for this point                           */
/*                                                                      */
/*      Inputs:     st1		Station information structures          */	
/*		    st2							*/
/*		    start       Start time of integration       	*/
/*		    freq        Frequency at center of channel		*/
/*                  obsparam    Observing parameters, e.g. intg time    */
/*                  scan_index  number of scan in obsparam.nscan        */
/*                  integration_number index of the time slice          */
/*		    cellsize    (u,v) grid cellsize in wavelengths  	*/
/*                  beamconvl   Station beam convolving function        */
/*                                                                      */
/*      Output:     udim, vdim  Patch sizes in wavelengths              */
/*                                                                      */
/* Created 16 Sept 02 by SSD                                            */
/* Corrected 22 Nov 2012 by L. Benkevitch:        			*/
/*    The description of parameters brought to consistency with the     */
/*    actual parameters                                                 */
/*                                                                      */
/************************************************************************/
#include <stdlib.h>
#include <math.h>
#include "utility.h"
#include "observing_param.h"
#include "convolve.h"
#include "array_specification.h"

#define max(a,b) ((a)>(b) ? (a) : (b))
#define min(a,b) ((a)<(b) ? (a) : (b))

int
patch_size (/* st1, st2, start, freq, obsparam, scan_index, integration_number,
	     * cellsize, *beamconvl, *udim, *vdim */
	    struct stationspec *st1,
	    struct stationspec *st2,
	    double start,
	    double freq,
	    struct observing_param *obsparam,
	    int scan_index,
	    int integration_number,
	    double cellsize,
	    struct conv_fn *beamconvl,
	    double *udim,
	    double *vdim)
{
    int  ret1=0, ret2=0;
    double start_freq, stop_freq; 
    double u1,v1,u2,v2,u3,v3,u4,v4;
    double u_max, v_max, u_min, v_min;
    double u_len1, u_len2, v_len1, v_len2, w_len1, w_len2;
    double gha1, gha2;
    double convlx_max, convly_max;


	    /* Get start and stop frequencies	*/
    start_freq = freq - (obsparam->corr_chan_bw)/2.0;
    stop_freq = freq + (obsparam->corr_chan_bw)/2.0;
 
    /* printf("patch_size: freq=%.7e, start_freq=%.7e, stop_freq=%.7e\n", */
    /*        freq, start_freq, stop_freq);                               */

    	    /* Get start and stop gha		*/

/*     time = (start + stop)/2.0; */
/*     ret = compute_gha (obsparam, time, &gha); */
    if (obsparam->scan[scan_index].gha_used) {
      gha1 = obsparam->scan[scan_index].gha_start + 
	integration_number*obsparam->
	    integ_time/3600.0*M_PI/12.0*(24.0/23.9344696);
      gha2 = obsparam->scan[scan_index].gha_start + 
	(integration_number+1)*obsparam->
	    integ_time/3600.0*M_PI/12.0*(24.0/23.9344696);
    } else {
      ret1 = compute_gha (obsparam, start, &gha1);
      ret2 = compute_gha (obsparam, start+obsparam->integ_time, &gha2);
      if ((ret1!=0)||(ret2!=0))  {
        msg ("Failure in compute_gha in patch_size", 2);
        exit(1);
      }
    }
    
    /* Get u,v,w in meters for gha1 and gha2	*/
/*     ret = compute_uvw(st1, st2, gha, obsparam, &u_len, &v_len, &w_len); */
    ret1 = compute_uvw(st1, st2, gha1, obsparam, &u_len1, &v_len1, &w_len1);
    ret2 = compute_uvw(st1, st2, gha2, obsparam, &u_len2, &v_len2, &w_len2);
    /* we have already checked the elevation limits, so it doesn't really  *
     * matter if we cross it during an integration                         */
//    if ((ret1!=0)||(ret2!=0)) {
//        msg ("Failure in compute_uvw in patch_size", 2);
//        exit(1);
//   }

/* printf("patch_size: u_len = %.7e, v_len = %.7e, w_len = %.7e\n",
        u_len, v_len, w_len);                                    
    printf("patch_size: u_len1 = %.7e, v_len1 = %.7e, w_len1 = %.7e\n", 
	   u_len1, v_len1, w_len1);
    printf("patch_size: u_len2 = %.7e, v_len2 = %.7e, w_len2 = %.7e\n", 
	   u_len2, v_len2, w_len2); */

	    /* Get u,v pairs for each of the four bounding points	*/
	    /* Units are wavelengths					*/
/*    u = u_len * freq * 1.0e6/VLIGHT;
    v = v_len * freq * 1.0e6/VLIGHT; 
    printf("patch_size: u=%.7e, v=%.7e\n",u,v); */
    u1 = u_len1 * start_freq * 1.0e6 / VLIGHT;
    v1 = v_len1 * start_freq * 1.0e6 / VLIGHT;
    u2 = u_len1 * stop_freq * 1.0e6 / VLIGHT;
    v2 = v_len1 * stop_freq * 1.0e6 / VLIGHT;
    u3 = u_len2 * start_freq * 1.0e6 / VLIGHT;
    v3 = v_len2 * start_freq * 1.0e6 / VLIGHT;
    u4 = u_len2 * stop_freq * 1.0e6 / VLIGHT;
    v4 = v_len2 * stop_freq * 1.0e6 / VLIGHT;

    	    /* Find limits of enclosing rectangle	*/
    u_max = max(max(u1,u2),max(u3,u4));
    v_max = max(max(v1,v2),max(v3,v4));
    u_min = min(min(u1,u2),min(u3,u4));
    v_min = min(min(v1,v2),min(v3,v4));

    	    /* Get patch size, padding by 4 cells for integrate *
	     * interpolation	                                */
    *udim = u_max - u_min + 4.0*cellsize;
    *vdim = v_max - v_min + 4.0*cellsize;

    /* printf("u_max = %g, u_min = %g, v_max = %g, v_min = %g\n",  */
    /* 	   u_max, u_min, v_max, v_min); */
    /* printf("u_max - u_min = %g;   v_max - v_min = %g\n",  */
    /* 	   u_max-u_min, v_max-v_min); */
    
    /* printf("u1, v1, u2, v2, u3, v3, u4, v4 = %g %g %g %g %g %g %g %g\n",  */
    /* 	   u1, v1, u2, v2, u3, v3, u4, v4); */
    //printf("*udim = %12.1f, *vdim = %12.1f\n", *udim, *vdim);


    /* Now pad by size of largest convolution dimensions        *
     * (in wavelengths). Note that largest dimensions can come  *
     * from different convolving fn's                      	*/

    convlx_max = beamconvl->xsize * cellsize;
    convly_max = beamconvl->ysize * cellsize;

    /* printf("beamconvl->xsize = %d, cellsize = %g, convlx_max = %g\n", */
    /* 	   beamconvl->xsize, cellsize, convlx_max); */
    /* printf("beamconvl->ysize = %d, cellsize = %g, convly_max = %g\n", */
    /* 	   beamconvl->ysize, cellsize, convly_max); */

    *udim += convlx_max;
    *vdim += convly_max;

    

    /* printf("patch_size: udim = %.3f\tvdim = %.3f\n", *udim, *vdim); 
    printf("u_max=%.5e, u_min=%.5e, v_max=%.5e, v_min=%.5e\n",
           u_max, u_min, v_max, v_min);
    printf("midpoint = %.5f, %.5f\n", (u_max+u_min)/2.0, (v_max+v_min)/2.0); */
    
    return (0);
    }
