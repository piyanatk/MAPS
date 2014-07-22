/************************************************************************/
/*                                                                      */
/* This routine returns the number of seconds since midnight on         */
/* Jan 1st 2000 represented by the date struct contents.  This makes    */
/* for fast, easy checks on timeranges etc. using a compact, single-    */
/* variable time representation in seconds.  This routine handles only  */
/* dates in the 21st century.                                           */
/*                                                                      */
/*      Inputs:         sim_date        date struct                     */
/*                                                                      */
/*      Output:         return value    Double precision                */
/*                                                                      */
/* Modelled after Mk4 routine of same name, 7 Jan 2002 by CJL           */
/*                                                                      */
/************************************************************************/
#include <stdio.h>
#include "visibility.h"
#include "utility.h"

double
time_to_double (/* sim_date) */
struct date sim_date)
{
    int year, nleaps; 
    double secs_since_00;

    msg ("time_to_double - yr,day,hr,min,sec = %d %d %d %d %f", -2, sim_date.year,
                    sim_date.day,
                    sim_date.hour,
                    sim_date.minute,
                    sim_date.second);

    year = sim_date.year % 100;
    nleaps = (year + 3) /4;   /* Only count if past the leap year */

    secs_since_00 = year * 31536000.0 
                    + (sim_date.day+nleaps-1) * 86400.0
                    + sim_date.hour * 3600.0 
                    + sim_date.minute * 60.0 
                    + sim_date.second;

    return (secs_since_00);
}
