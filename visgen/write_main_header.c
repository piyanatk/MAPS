/************************************************************************/
/*                                                                      */
/* Given information in the obsparam struct, this routine writes the    */
/* main header into the visibility output file, using the format        */
/* laid out in visibility.h                                             */
/*                                                                      */
/*      Inputs:     obsparam    Filled observing_param struct           */
/*                  vgparam     General program information             */
/*                  outfp       Open file pointer for output file       */
/*                                                                      */
/*      Output:     Data written into output file                       */
/*                  return value    0=OK, else bad                      */
/*                                                                      */
/* Created 2 Feb 2002 by CJL                                            */
/*                                                                      */
/************************************************************************/
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include "utility.h"
#include "observing_param.h"
#include "visgen.h"
#include "visibility.h"
#include "array_specification.h"

int
write_main_header (/* obsparam, vgparam, outfp) */
struct observing_param *obsparam,
struct vg_param *vgparam,
struct array_specification *arrayspec,
FILE *outfp)
{
    struct main_header mh;
    int i, n, size;
    double scantime, earliest,ra,dec;

    strncpy (mh.simname, vgparam->simname,SIZ_SIMNAME);

    ra  = obsparam->phase_cent_RA*(12.0/M_PI);   /* ra in hours */
    dec = obsparam->phase_cent_DEC*(180.0/M_PI); /* dec in degs */

    mh.field_center.ra_hrs   = floor(ra);
    mh.field_center.ra_mins  = floor((ra-mh.field_center.ra_hrs)*60.0);
    mh.field_center.ra_secs  = (ra-mh.field_center.ra_hrs-mh.field_center.ra_mins/60.0)*3600.0;
    mh.field_center.dec_degs = floor(fabs(dec))*(dec < 0? -1 : 1);
    mh.field_center.dec_mins = floor(fabs(dec-mh.field_center.dec_degs)*60.0)*(dec < 0? -1 : 1);
    mh.field_center.dec_secs = (dec-mh.field_center.dec_degs-mh.field_center.dec_mins/60.0)*3600.0;

    for (i=0; i<8; i++) mh.stokes[i] = '\0';
    for(i=0; i<vgparam->n_pol_products; i++) {
      mh.stokes[2*i  ] = arrayspec->station[0].layout->ants[0].pol_type[i/2];
      mh.stokes[2*i+1] = arrayspec->station[0].layout->ants[0].pol_type[i%2];
    }

    strncpy (mh.field_name,"Undefined" ,SIZ_FIELDNAME);
    mh.ntime_block = obsparam->nscan;
                                        // Find reftime - use start of
                                        // earliest scan, so all offsets
                                        // are positive
    earliest = FLT_MAX;
    for (i=0; i<obsparam->nscan; i++) {
      scantime = time_to_double (obsparam->scan[i].start);
      if (scantime < earliest) earliest = scantime;
    }
    double_to_time(earliest, &(mh.reftime));
                                        // Record in program struct
    vgparam->reftime = earliest;

    size = sizeof (struct main_header);
    msg ("Writing main header, size = %d bytes", -1, size);
    n = fwrite (&mh, size, 1, outfp); 
    if (n != 1) {
      msg ("Error writing main header of output file", 3);
      return (1);
    }
    return (0);
}
