/************************************************************************/
/*                                                                      */
/* Based on P. Sherwood's ObservationRead.c, simplified to better fit   */
/* this specialized & relatively static application.                    */
/* Modified to fit visgen architecture.  Adapted for style, to use      */
/* standard utilities, and to store observation parameters in a         */
/* structure instead of in globals.  Also changed to support more       */
/* complex observation specifications, with multiple scans/channels     */
/*                                                                      */
/*      Inputs:         filename        Name of input file              */
/*                                                                      */
/*      Output:         obsparam        Filled observing_param struct   */
/*                      return value    0=OK, else bad                  */
/*                                                                      */
/* Created 2 Feb 2002 by CJL                                            */
/*                                                                      */
/************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <fitsio.h>
#include "novas.h"
#include "utility.h"
#include "observing_param.h"
#include "visgen.h"

#define TRUE 1
#define FALSE 0

#define EMSG(text) {msg (text, 2); return (-1);}

#define MAX_LINE 256
#define GHA_TAG "GHA"
#define NUM_FIELDS 15
                                        // Match with switch statement below
char *fields[NUM_FIELDS] = {"FOV_center_RA",
                    "FOV_center_Dec",
                    "FOV_size_RA",
                    "FOV_size_Dec",
                    "Scan_start",
                    "Scan_duration",
                    "Channel",
                    "Endscan",
                    "Corr_int_time",
                    "Corr_chan_bw",
                    "Elevation_limit",
                    "Time_cells",
                    "Freq_cells",
                    "PNT_center_RA",
                    "PNT_center_Dec",
                    };

// Function to check bandwidths per channel match
int checkIFBandwidth(struct observing_param *obsparam);

int read_obs_spec(char *filename, 
		  struct observing_param *obsparam, 
		  struct vg_param *vgparam, 
		  double array_lon_radian )
{
    char line[MAX_LINE], value[MAX_LINE], fieldname[MAX_LINE], junk[MAX_LINE], *ptr;
    int i, j, flen, nline, nscan, nchan=0, y, d, h, m;
    int match, fovr, fovd, fovsr, fovsd, cint, cbw, elev, error, in_scan;
    int tcells, fcells, timecell, freqcell, status=0;
    float s;
    double dval, freq, bw, diff=0;
    fitsfile *fitsfh=NULL;
    FILE *obsfile;
    time_t  now;
    struct tm gmt;
    char fitsname[128];
    long naxes[16], n_cols, n_rows; /* Dimensions of image plane */ 

                                        // Initialize
    obsparam->point_cent_RA = obsparam->point_cent_DEC = -99;
    obsparam->phase_cent_RA = obsparam->phase_cent_DEC = -99;
    obsparam->fov_RA = obsparam->fov_dec = -1.0;
    obsparam->integ_time = obsparam->corr_chan_bw = -1.0;
//    obsparam->el_limit = 10.0 / 57.29577951308;     // 10 degree default
    obsparam->nscan = 0;
    obsparam->nt_cells = 0; // by default, don't do averaging in sub time/freq intervals in correlator
    obsparam->nf_cells = 0;
    nscan = 0;
    for (i=0; i<MAXSCAN; i++) obsparam->scan[i].nfreq = 0;

                                        // Open the observation file
    if ((obsfile = fopen (filename, "r")) == NULL) {
        msg ("Can't find the observation file '%s'", 2, filename);
        return (1);
    }
                                        // Process the lines of the observation file
    nline = 0;
    fovr = fovd = fovsr = fovsd = cint = cbw = elev = timecell = freqcell = in_scan = FALSE;
    while (TRUE) {
                                        // get line from the file, stop on error or EOF
        if ((fgets (line, MAX_LINE-1, obsfile)) == NULL) break; 
        nline++;
                                        // Remove leading/trailing blanks, tabs, newlines
                                        // and ignore blank lines or lines starting "//"
        rm_whitespace (line);
        msg ("obsfile %3d: '%s'", -1, nline, line);
        if (strlen (line) == 0) continue;
        if (line[0]=='/' && line[1]=='/') continue;
                                        // Lines in the observation file have the format
                                        // <fieldname> = <value>
                                        // The blanks around the "=" are optional
                                        // '=' optional in case of endscan statement
        if ((ptr = strchr (line, '=')) == NULL) {
            ptr = line + strlen (line);
            *(ptr+1) = '\0';            // Avoid problems with value strcpy below
            }
        flen = ptr - line;
        strncpy (fieldname, line, flen); 
        fieldname[flen] = '\0'; 
        rm_whitespace (fieldname);
        strcpy (value, ptr + 1);        // skip the '=' 
        rm_whitespace (value);  
                                        // Identify the field, minmatch
        match = -1;
        for (i=0; i<NUM_FIELDS; i++) {
            if (strncmp (fieldname, fields[i], flen) == 0)
                {
                if (match >= 0)
                    {
                    msg ("Ambiguous keyword '%s'", 2, fieldname);
                    return (-1);
                    }
                match = i;
                }
            }
                                        // Process values
        switch (match) {
            case 0:                     // FOV_center_RA  (hr:min:sec)

                if (in_scan || fovr) EMSG ("FOV_center_RA misplaced or repeated")
                if (sscanf (value, "%d:%d:%f %s", &h, &m, &s, junk) != 3) EMSG ("Bad RA format, should be h:m:s");
                if (h>23 || h<0 || m>59 || m<0 || s>=60.0 || s<0.0) EMSG ("RA specified outside of 0h - 24h range");
                                /* convert to radian now */
                obsparam->phase_cent_RA = (h + m/60.0 + s/3600.0)*(M_PI/12.0);
                fovr = TRUE;
                break;

            case 1:                     // FOV_center_Dec  (deg:min:sec)
                if (in_scan || fovd) EMSG ("FOV_center_RA misplaced or repeated")
                if (sscanf (value, "%d:%d:%f %s", &d, &m, &s, junk) != 3)
                    EMSG ("Bad Dec format, should be d:m:s")
                if (h>90 || h<-90 || m>59 || m<0 || s>=60.0 || s<0.0) 
                    EMSG ("Degrees:minutes:seconds limited to +/- 90")
                                        // Prevent infamous negative dec error
                if (value[0] == '-') {m = -m; s = -s;}
                                /* convert to radian now */
                obsparam->phase_cent_DEC = (d + m/60.0 + s/3600.0)*(M_PI/180.0);
                fovd = TRUE;
                break;

            case 2:                     // FOV_size_RA  (arcsec)
                if (in_scan || fovsr) EMSG ("FOV_size_RA misplaced or repeated")
                if (sscanf (value, "%lf %s", &dval, junk) != 1) 
                    EMSG ("Bad FOV_size_RA field")
                if (dval <= 0.0) EMSG ("FOV_size_RA must be positive")
                obsparam->fov_RA = dval;
                fovsr = TRUE;
                break;

            case 3:                     // FOV_size_Dec  (arcsec)
                if (in_scan || fovsd) EMSG ("FOV_size_Dec misplaced or repeated")
                if (sscanf (value, "%lf %s", &dval, junk) != 1) 
                    EMSG ("Bad FOV_size_Dec field")
                if (dval <= 0.0) EMSG ("FOV_size_Dec must be positive")
                obsparam->fov_dec = dval;
                fovsd = TRUE;
                break;

            case 4:                     // Scan_start  (yyyy:ddd:hh:mm:ss.sss)
                /* special case: use GHA rather than absolute time */
                if (!strncmp(value,GHA_TAG,strlen(GHA_TAG))) {
                    obsparam->scan[nscan].gha_used = 1;
                    obsparam->scan[nscan].gha_start = atof(value+strlen(GHA_TAG));
                    if (obsparam->scan[nscan].gha_start < -12 || obsparam->scan[nscan].gha_start > 12) {
                      msg("invalid GHA: %g. must be between -12 and 12 hours",2,obsparam->scan[nscan].gha_start);
                      return -1;
                    }
                    /* convert to radian */
                    obsparam->scan[nscan].gha_start *= (M_PI/12.0);
     
                    /* set the observing params for now, so that the
                       date is meaningul later as a UVFITS file

                        Needed to correct this to be a valid time for the
                        GHA - otherwise MIRIAD and other tools 
                        get the wrong LST*/ 

                    if (nscan==0) {

                        struct tm faketime;
                        double jd=0,gmst=0,gmst_0=0,dummy,ee=0.0;

                        faketime.tm_year = 2007 - 1900;
                        faketime.tm_mon = 8;
                        faketime.tm_mday = 21;
                        faketime.tm_hour = 0;
                        faketime.tm_min = 0;
                        faketime.tm_sec = 0;

                        /* begin correction */
                        /* we know the desired GHA and we have made a fake observing date.
                           We want to know an absolute (UT) time that corresponds to the GHA
                           on the fake date. Find the GST at midnight on the fake date.
                           Then the GHA at midnight is the GST + RA of phase centre.
                           The difference (gha_start - GHA_midnight)
                           must be converted to seconds and added to the fake time to get the
                           absolute time of the observation */

                        /* first we need the JD */
                        /* remember months are counted from zero but mday is counted from 1 */
                        jd = julian_date(faketime.tm_year+1900,faketime.tm_mon+1,faketime.tm_mday,0.0);

                        /* get the Equation of Equinoxes offset time in seconds */
                        earthtilt(jd, &dummy, &dummy, &ee, &dummy, &dummy);

                        /* Greenwich apparent sidereal time at midnight, between 0 and 24 hours */
                        sidereal_time(jd,0.0,ee,&gmst_0);
                        gmst_0 *= (M_PI/12.0); // convert to radian

                        /* calc the required GST given the phase center RA and desired GHA */
                        gmst = obsparam->phase_cent_RA + obsparam->scan[0].gha_start;
                        /* make it between 0 and 2pi */
                        if (gmst < 0)      gmst += 2*M_PI;
                        if (gmst > 2*M_PI) gmst -= 2*M_PI;

                        // LST offset between desired gmt and gmt of fake date
                        diff = gmst-gmst_0;

                        /* convert to time offset in seconds */
                        diff /= 2.0*M_PI;    // radian to sidereal days
                        diff *= 86164.0905;  // sidereal days to seconds 

                        now = timegm(&faketime);

                        // converting to gmt which doesn't have <1sec res, so add the difference back later.
                        now += floor(diff);

                        printf( "JD: %f, ee (sec): %g, GMST0: %f (hrs), GMST (desired): %f (hrs), RA: %f, UT1_OFF: %.1f (seconds), now: %ld \n",
                                jd, ee, gmst_0*(12.0/M_PI), gmst*(12.0/M_PI), obsparam->phase_cent_RA, diff, (long) now );

                        /* print the adjusted UT date */
                        printf("Reverse engineered observation start time: %s",asctime(gmtime(&now)));

                    }
                    else {
                      /* for multiple scans, set up the time to be last time plus offset */
                      diff = (obsparam->scan[nscan].gha_start-obsparam->scan[nscan-1].gha_start)*(3600.0*12.0/M_PI);
                      printf("gha1: %g, gha2: %g, diff: %g\n",obsparam->scan[nscan].gha_start,
                             obsparam->scan[nscan-1].gha_start,diff);
                      // converting to gmt which doesn't have <1sec res, so add the difference back later.
                      now += floor(diff);
                    }
                    gmtime_r(&now,&gmt);

                    obsparam->scan[nscan].start.year = gmt.tm_year+1900;
                    obsparam->scan[nscan].start.day  = gmt.tm_yday;
                    obsparam->scan[nscan].start.hour = gmt.tm_hour;
                    obsparam->scan[nscan].start.minute = gmt.tm_min;
                    /* add the fractional second time offset from to the start time */
                    obsparam->scan[nscan].start.second = gmt.tm_sec + (diff-floor(diff));

                    msg("Scan: %d. GHA: %g (radian) Year: %d, Day: %d H:M:S (%d:%d:%g)",1,nscan,obsparam->scan[nscan].gha_start,
                        (int)obsparam->scan[nscan].start.year, (int)obsparam->scan[nscan].start.day,
                        (int)obsparam->scan[nscan].start.hour, (int)obsparam->scan[nscan].start.minute,
                        obsparam->scan[nscan].start.second);
                } else {
                    if (sscanf (value, "%d:%d:%d:%d:%f %s", &y, &d, &h, &m, &s, junk) != 5) 
                      EMSG ("Bad date format, should be y:d:h:m:s")
                    if (y<2000 || y>2099 || d<1 || d>366 || h>23 || h<0  || m>59 || m<0 || s>=60.0 || s<0.0)
                      EMSG ("Invalid date specified")
                    obsparam->scan[nscan].gha_used =0;
                    obsparam->scan[nscan].gha_start=0;
                    obsparam->scan[nscan].start.year = y;
                    obsparam->scan[nscan].start.day  = d;
                    obsparam->scan[nscan].start.hour = h;
                    obsparam->scan[nscan].start.minute = m;
                    obsparam->scan[nscan].start.second = s;
                }
                nchan = 0;
                in_scan = TRUE;
                break;

            case 5:                     // Scan_duration  (seconds)
                if (! in_scan) EMSG ("Scan_duration statement outside of scan")
                if (sscanf (value, "%lf %s", &dval, junk) != 1) 
                    EMSG ("Bad Scan_duration field")
                if (dval <= 0.0) EMSG ("Scan_duration must be positive")
                obsparam->scan[nscan].duration = dval;
                break;

            case 6:                     // Channel  (freq:bw  MHz)
                if (! in_scan) EMSG ("Channel statement outside of scan")
                if (sscanf (value, "%lg:%lg %s", &freq, &bw, junk) != 2) EMSG ("Bad channel field")
                if (freq<10) EMSG ("Bad channel data")
                obsparam->scan[nscan].freq[nchan].frequency = freq;
                obsparam->scan[nscan].freq[nchan].bandwidth = bw;
                nchan++;
                break;

            case 7:                     // Endscan
                if (! in_scan) EMSG ("Endscan statement out of place")
                obsparam->scan[nscan].nfreq = nchan;
                nscan++;
                in_scan = FALSE;
                break;

            case 8:                     // Corr_int_time  (seconds)
                if (in_scan || cint) EMSG ("Corr_int_time misplaced or repeated")
                if (sscanf (value, "%lf %s", &dval, junk) != 1) 
                    EMSG ("Bad Corr_int_time field")
                if (dval <= 0.0) EMSG ("Corr_int_time must be positive")
                obsparam->integ_time = dval;
                cint = TRUE;
                break;

            case 9:                     // Corr_chan_bw  (Mhz)
                if (in_scan || cbw) EMSG ("Corr_chan_bw misplaced or repeated")
                if (sscanf (value, "%lf %s", &dval, junk) != 1) 
                    EMSG ("Bad Corr_chan_bw field")
                if (dval <= 0.0) EMSG ("Corr_chan_bw must be positive")
                obsparam->corr_chan_bw = dval;
                cbw = TRUE;
                break;

            case 10:                     // Elevation_limit (degrees)
                if (in_scan || elev) EMSG ("Elevation_limit misplaced or repeated")
                msg("WARNING: el_limit no longer specified in obs_spec file. Use station spec file. This param ignored.\n",3);
                break;

            case 11:                        // Set nt_cells
                if (in_scan || timecell) EMSG ("Time_cells misplaced or repeated")
                if (sscanf (value, "%d %s", &tcells, junk) != 1)
                    EMSG ("Bad Time_cells field");
                if (tcells < 0 || tcells%2 !=0 ) EMSG ("Time_cells must be even and >= 0.")
                obsparam->nt_cells = tcells;
                msg ("Successfully entered Time_cells into structure: tcells = %d", 0, tcells);
                timecell = TRUE;
                break;

            case 12:                        // Set nf_cells
                if (in_scan || freqcell) EMSG ("Freq_cells misplaced or repeated")
                if (sscanf (value, "%d %s", &fcells, junk) != 1)
                    EMSG ("Bad Freq_cells field")
                if (fcells < 0 || fcells%2 != 0 ) EMSG ("Freq_cells must be even and >= 0.")
                obsparam->nf_cells = fcells;
                msg ("Successfully entered Freq_cells into structure: fcells = %d", 0, fcells);
                freqcell = TRUE;
                break;

            case 13:                        // Set point_cent_RA
                if (sscanf (value, "%lf %s", &dval, junk) != 1) EMSG ("Bad PNT_center_RA field");
                if (dval < 0 || dval > 24 ) EMSG ("PNT_cent_RA must be between 0 and 24 ")
                obsparam->point_cent_RA = dval*(M_PI/12.0);
                break;
 
            case 14:                        // Set point_cent_DEC
                if (sscanf (value, "%lf %s", &dval, junk) != 1) EMSG ("Bad PNT_center_DEC field");
                if (dval < -90 || dval > 90 ) EMSG ("PNT_cent_DEC must be between -90 and 90")
                obsparam->point_cent_DEC = dval*(M_PI/180.0);
                break;
 
        default:                    // Not recognized
                msg ("Unrecognized field '%s', ignoring", 2, fieldname);
                break;
	}   /* switch (match) { ... */
    }    /* while (TRUE) { ... */
    obsparam->nscan = nscan;
    fclose(obsfile);


    /*
     * Modification by L. Benkevitch:
     * The values for "FOV_size_RA" and "FOV_size_Dec" now do not 
     * have to be specified in the obs_spec file. Normally, these values
     * are calculated from the fits file attributes CDELT1 and CDELT2,
     * assumed to be in degrees/pixel.
     */
    if (! (fovsr && fovsd)) { 
      strncpy(fitsname, vgparam->simname, VISGEN_SIZE_NAME);
      strncat(fitsname, ".fits", 5); 
      fits_open_image(&fitsfh, fitsname, READONLY, &status);
      if (status !=0) {
	fprintf(stderr, "FATAL: Failed to open file <%s> in FITSio. "	\
		"Error code: %d\n", fitsname, status); exit(1); }
      fits_get_img_size(fitsfh, 16, naxes, &status);
      if (status != 0) {
	fprintf(stderr,"FATAL: Failed to get axis sizes with error code %d\n",
		status); exit(1); }
      n_cols = naxes[0];
      n_rows = naxes[1];
      fits_read_key(fitsfh, TDOUBLE, "CDELT1", &dval, NULL, &status);
      if (status!=0) {
	fprintf(stderr,"FATAL: Error reading CDELT1 from header.\n"); exit(1); }
      obsparam->fov_RA = fabs(dval)*n_cols*3600; /* Express in arcseconds */
      fits_read_key(fitsfh, TDOUBLE, "CDELT2", &dval, NULL, &status);
      if (status!=0) { 
	fprintf(stderr,"FATAL: Error reading CDELT2 from header.\n"); exit(1); }
      obsparam->fov_dec = fabs(dval)*n_rows*3600; /* Express in arcseconds */
      fits_close_file(fitsfh, &status);
      if (status !=0) {
	fprintf(stderr,"WARNING: fits_close_file() failed with code %d. " \
		"Do we care?\n", status); exit(1); }
      printf("From '%s' file: FOV size = %g x %g arcsec\n", fitsname, 
	     obsparam->fov_RA, obsparam->fov_dec);

    }
    /* End of modification by L. Benkevitch */

                                        // Check that information is complete
    error = FALSE;
    /* if (! (fovr && fovd && fovsr && fovsd)) */
    /*     {msg ("Field of view incompletely specified", 2); error = TRUE;} */
    if (! (fovr && fovd))
        {msg ("Field of view center incompletely specified", 2); error = TRUE;}
    if (! (cint && cbw))
        {msg ("Correlator parameters incompletely specified", 2); error = TRUE;}
    if (nscan == 0)
        {msg ("No scans specified", 2); error = TRUE;}
    for (i=0; i<nscan; i++)
        {
        if (obsparam->scan[i].nfreq == 0)
            {msg ("Scan %d incompeletely specified", 2, i); error = TRUE;}
        }
//    if (! elev) msg ("No elevation limit specified, using default %g degrees", 2,obsparam->el_limit*(180.0/M_PI));
    if (! timecell) msg ("No number of time cells per visibility point specified, using default of %d",2,obsparam->nt_cells);
    if (! freqcell) msg ("No number of freq cells per visibility point specified, using default of %d",2,obsparam->nf_cells);
    if (obsparam->point_cent_RA == -99.0) obsparam->point_cent_RA = obsparam->phase_cent_RA;
    if (obsparam->point_cent_DEC == -99.0) obsparam->point_cent_DEC = obsparam->phase_cent_DEC;

	// Check bandwidths for multiple IFs
	if(checkIFBandwidth(obsparam) != 0) return -1;
	
                                        // Diagnostic print
    msg ("Point center (RA,DEC radian)= %g, %g", 2, obsparam->point_cent_RA, obsparam->point_cent_DEC);
    msg ("Phase center (RA,DEC radian)= %g, %g", 2, obsparam->phase_cent_RA, obsparam->phase_cent_DEC);
    msg ("FOV size = %g x %g arcsec", 2, obsparam->fov_RA, obsparam->fov_dec);
    msg ("Corr. params = %g secs, %g MHz", 2, obsparam->integ_time, obsparam->corr_chan_bw);
    msg ("Number of scans = %d", 2, obsparam->nscan);
    for (i=0; i<nscan; i++)
        {
        msg ("Scan %d, start = %d:%d:%d:%d:%g, duration = %g, nfreq = %d", 1, i,
                    obsparam->scan[i].start.year, obsparam->scan[i].start.day,
                    obsparam->scan[i].start.hour, obsparam->scan[i].start.minute,
                    obsparam->scan[i].start.second, obsparam->scan[i].duration,
                    obsparam->scan[i].nfreq);
        for (j=0; j<obsparam->scan[i].nfreq; j++)
            msg ("  chan %d, freq = %g MHz, bw = %g MHz", 1, j,
                    obsparam->scan[i].freq[j].frequency, 
                    obsparam->scan[i].freq[j].bandwidth);
        }

    if (error) return (-1);
    return (0);
    }

// Check that each IF in a particular scan has the same bandwidth, otherwise, can't generate UVFITS
// (the way multiple IF's groups are written sequentially to UVFITS requires that they all be the same lenth,
// since the correlator bandwidth doesn't change, only need to check that the bandwidth is equal)
int checkIFBandwidth(struct observing_param *obsparam) {
	int i, j;
	double bw;
	for(i=0; i<obsparam->nscan; i++) {
		bw=-1.0;
		for(j=0; j<obsparam->scan[i].nfreq; j++) {
			if(bw > 0) {
				if(bw != obsparam->scan[i].freq[j].bandwidth) {
					msg("Bandwidth %g did not match previous bandwidth %g. Channels must have same bandwidth.\n", 2,
						obsparam->scan[i].freq[j].bandwidth, bw);
					return -1;
				}
			}
			
			// Store to compare to next channel
			bw = obsparam->scan[i].freq[j].bandwidth;
		}
	}
	
	// If no clashes, return successfully
	return 0;
}
