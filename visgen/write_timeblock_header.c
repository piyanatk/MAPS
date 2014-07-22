/************************************************************************/
/*                                                                      */
/* For each scan, the output visibility file needs a new time block     */
/* header, which this routine appends to the file.                      */
/*                                                                      */
/*      Inputs:     obsparam    Filled observing parameter struct       */
/*                  scan_no     Scan number in obsparam                 */
/*                  vgparam     General program information             */
/*                  outfp       File pointer for output file            */
/*                                                                      */
/*      Output:     Appended header record                              */
/*                  return value    0=OK, else bad                      */
/*                                                                      */
/* Created 2 Feb 2002 by CJL                                            */
/* Fixed bug in computation of tbh->start_time, CJL, 17 April 2003      */
/* Added an explicit declaration of time_to_double, CJL, 18 April 2003  */
/* moved decl of time_to_double to utility.h May 2007 RBW.              */
/************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include "utility.h"
#include "observing_param.h"
#include "visgen.h"
#include "visibility.h"
#include "sizes.h"

int
write_timeblock_header (/* obsparam, scan_no, vgparam, outfp) */
struct observing_param *obsparam,
int scan_no,
struct vg_param *vgparam,
FILE *outfp)
{
    int i, n, nchan, size;
    double scantime;
    void *ptr;
    struct scaninfo *scan;
    struct time_block_header *tbh;

    msg ("Writing timeblock header for scan %d", -1 , scan_no);
    scan = obsparam->scan + scan_no;
                                        // Allocate memory for time block header
    nchan = scan->nfreq;
    size = size_of_tbh (nchan);
    msg ("Size of timeblock header = %d bytes", -1, size);
    if ((ptr = (void *)malloc (size)) == NULL)
        {
        msg ("Error allocating memory for time block header", 2);
        return (1);
        }
    tbh = (struct time_block_header *)ptr;
                                        // Start filling in data
    scantime = time_to_double(scan->start);
    msg ("scan start = %g, reftime=%g", -2, scantime, vgparam->reftime);
    tbh->start_time = scantime - vgparam->reftime;
    tbh->duration = scan->duration;
    tbh->intg_time = obsparam->integ_time;
    tbh->num_IFs = nchan;
    msg ("time, duration, intg_time and nchan = %g %g %g %d", -2,
                    tbh->start_time, tbh->duration, tbh->intg_time, tbh->num_IFs);
                                        // Variable length frequency table
    for (i=0; i<nchan; i++) {
      tbh->ftable[i].frequency = scan->freq[i].frequency;
      tbh->ftable[i].bandwidth = scan->freq[i].bandwidth;
      msg ("channel %d, freq, bw = %g %g", -2, i, tbh->ftable[i].frequency,
                                                    tbh->ftable[i].bandwidth);
    }
                                        // Do the IO
    n = fwrite (tbh, size, 1, outfp);
    if (n != 1)
        {
        msg ("Error writing time block header for scan %d", 2, scan_no);
        free (tbh);
        return (1);
        }

    free (tbh);
    return (0);
}
