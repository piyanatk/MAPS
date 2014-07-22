/************************************************************************/
/*                                                                      */
/* This routine takes a string identifying an array site, and looks     */
/* for that site in $TEXTDIR/sites.txt, parsing latitude and longitude  */
/* from the corresponding entry in that file.  For array specifications */
/* based on relative coordinates for stations, the site location        */
/* is used in computation of earth-centered coordinates for antennas    */
/* and stations.  Site locations obtained through this routine          */
/* override site locations obtained from array specification files.     */
/*                                                                      */
/*      Inputs:     sitename        keyword in sites.txt file           */
/*                                                                      */
/*      Output:     arrayspec       array center coordinates updated    */
/*                  return value    0=OK, else bad                      */
/*                                                                      */
/* Created Jan 8th 2002 by CJL                                          */
/*                                                                      */
/************************************************************************/
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "array_specification.h"

#define FALSE 0
#define TRUE 1


int
get_site (/* sitename, arrayspec) */
	  char *sitename,
	  struct array_specification *arrayspec) {
  char sitefile[FILENAME_MAX], keyword[32], buf[BUFSIZ];
  char latit[64], longit[64];
  char *mtextdir=NULL;
  int found, min, deg, hour;
  double sec, angle;
  FILE *fp;

	/* get the environment variable "TEXTDIR" */
	mtextdir = getenv("TEXTDIR");
	if(mtextdir==NULL) mtextdir=".";
    sprintf (sitefile, "%s/sites.txt", mtextdir);
    if ((fp = fopen (sitefile, "r")) == NULL) {
        msg ("Could not open site file '%s'", 2, sitefile);
        return (1);
    }
    msg ("Opened the site file '%s'", 2, sitefile);

                                        /* Find line in text file */
    found = FALSE;
    while (fgets (buf, BUFSIZ, fp) != NULL) {
        if (buf[0] == '*') continue;
        if (strlen (buf) == 0) continue;
                                        /* Get first field */
        sscanf (buf, "%s %*s", keyword);
        if (strcmp (keyword, sitename) == 0){
            found = TRUE;
            break;
        }
    }
	fclose(fp);
    if (! found) {
        msg ("Could not find site '%s' in file '%s'", 2, sitename, sitefile);
        return (1);
    }
                                /* Parse lat/long and update struct */
    sscanf (buf, "%*s %s %s", longit, latit);
    sscanf (longit, "%d:%d:%lf", &deg, &min, &sec);
    if (longit[0] == '-') {min = -min; sec = -sec;}
    angle = deg + min/60.0 + sec/3600.0;
    hour = angle / 15.0;
    angle -= hour * 15.0;
    min = angle * 4;
    sec = (angle - min/4.0) * 240.0;

    arrayspec->location.long_hrs = hour;
    arrayspec->location.long_mins = min;
    arrayspec->location.long_secs = sec;
    sscanf (latit, "%d:%d:%lf", &deg, &min, &sec);
    if (latit[0] == '-') {min = -min; sec = -sec;}
    arrayspec->location.lat_degs = deg;
    arrayspec->location.lat_mins = min;
    arrayspec->location.lat_secs = sec;


                                // Convert array center lat/log to radians
                                // The location struct should have already 
                                // been filled in by get_site()
    arrayspec->location.lon_radian = ((double)arrayspec->location.long_hrs
              + (double)arrayspec->location.long_mins / 60.0
              + (double)arrayspec->location.long_secs / 3600.0)/3.819718634205;


    arrayspec->location.lat_radian = ((double)arrayspec->location.lat_degs
              + (double)arrayspec->location.lat_mins / 60.0
              + (double)arrayspec->location.lat_secs / 3600.0)/57.29577951308;
    msg ("Array is '%s', lat/long (radians): %.10g, %.10g", 1, 
	 sitename,arrayspec->location.lat_radian,
	 arrayspec->location.lon_radian);

    return (0);
    }
