/************************************************************************/
/*                                                                      */
/* Given a source RA and Dec, plus a date, this routine computes the    */
/* Greenwich hour angle, a quantity used many times in the innermost    */
/* loop of visgen.                                                      */
/*                                                                      */
/*      Inputs:     obsparam    Contains sky coordinates                */
/*                  time        Standard double format                  */
/*                                                                      */
/*      Output:     gha         In radians                              */
/*                                                                      */
/* Created 5 Feb 2002 by CJL                                            */
/* Ideas and code borrowed from PS, SSD and other CJL programs          */
/*                                                                      */
/************************************************************************/
#include <math.h>
#include "utility.h"
#include "observing_param.h"
#include "novas.h"

static double mod_pi(double x);

int
compute_gha (/* obsparam, time, gha) */
struct observing_param *obsparam,
double time,
double *gha)
    {
    double hour, jd, jd_int, jd_frac, ee, gst;
    static int nday[] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
    short month, day_of_month=0, day_of_year;
    struct date indate;

    /* special case: if we are specifying observation time in GHA,
       rather than absolute time, just return the GHA */

    // Decode time
    double_to_time(time, &indate);
                                        // Leap year
    if ((indate.year % 4) == 0) nday[1] = 29;
    else nday[1] = 28;
                                        // Convert doy to month/dom
    day_of_year = indate.day;
    for (month = 0; month < 12; month++)
        {
        if (day_of_year <= nday[month])
            {
            day_of_month = day_of_year;
            month += 1;                 // Convert to 1-relative
            break;
            }
        day_of_year -= nday[month];
        }
                                        // Compute julian date
    hour = (double)indate.hour + (double)indate.minute / 60.0 
                                + (double)indate.second / 3600.0;
    jd = julian_date (indate.year, month, day_of_month, hour);
                                        // Split jd for input to sidereal_time()
    jd_frac = modf (jd, &jd_int);
    if (isnan (jd))
        {
        msg ("Bad JD from julian_date()", 2);
        return (1);
        }
        
    ee = 0.0;
                                        // Compute Greenwich sidereal time in hours
    sidereal_time (jd_int, jd_frac, ee, &gst);

    if (isnan (gst))  {
	  msg ("Bad GST from sidereal_time()", 2);
	  return (1);
	}
                                        // Convert to radians
    gst /= 3.819718634205;
                                        // Subtract RA to get GHA

    *gha = mod_pi(gst - obsparam->phase_cent_RA);

    msg ("gha, gst = %g %g",-2,*gha,gst);

    return (0);
    }

/* make sure a number is between -pi and pi */
static double mod_pi(double x) {

  while (x > M_PI) {
    x -= 2.0*M_PI;
  }
  while (x <= -M_PI) {
    x += 2.0*M_PI;
  }
  return x;
}
