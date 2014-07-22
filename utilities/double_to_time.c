/************************************************************************/
/*                                                                      */
/* This routine receives the number of seconds since midnight on        */
/* Jan 1st 2000 and converts to a date struct.  It is the inverse of    */
/* the time_to_double() routine.                                        */
/*                                                                      */
/*      Inputs:         time             Seconds since Jan 1 2000       */
/*                                                                      */
/*      Output:         sim_date         Filled in date struct          */
/*                                                                      */
/* Created 7 Jan 2002 by CJL                                            */
/*                                                                      */
/************************************************************************/
#include <math.h>
#include "utility.h"

void
double_to_time (/* tim, sim_date) */
double time,
struct date *sim_date)
{
    int i, itime;
    double time_min;

    sim_date->second = fmod (time, 60.0);   /* Floating seconds */
    time_min = (time - sim_date->second) / 60.0;  /* Time in minutes */
    itime = rint (time_min);                /* Convert to integer for the rest */
    sim_date->minute = itime%60;
    itime /= 60;                            /* Time now in hours */
    sim_date->hour = itime%24;
    itime /= 24;                            /* Time now in days */
    for (i=0; i<99; i++) 
        {                                   /* Clumsy but effective */
        if (i%4 == 0 && itime <= 366) break;
        else if (itime <= 365) break;
        if (i%4 == 0) itime -= 366;
        else itime -= 365;
        }
    sim_date->year = i + 2000;
    sim_date->day = itime + 1;              /* Days are 1-relative */
}
