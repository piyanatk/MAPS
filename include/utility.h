/************************************************************************/
/*                                                                      */
/* General utility structures                                           */
/*                                                                      */
/* Created Jan 8 2002 by CJL                                            */
/*                                                                      */
/************************************************************************/
#ifndef UTILITY_H                         /* Handle multiple includes */
#define UTILITY_H    0

#include <stdio.h>

#define VLIGHT 299792458.

struct coord
    {
    short       ra_hrs;                 /* Self-explanatory */
    short       ra_mins;
    float       ra_secs;
    short       dec_degs;
    short       dec_mins;
    float       dec_secs;
    };

struct geocoord {
    short       long_hrs;               /* Self-explanatory */
    short       long_mins;
    float       long_secs;
    short       lat_degs;
    short       lat_mins;
    float       lat_secs;
  double      lon_radian;               /* lat and long radian */
  double      lat_radian;
};

struct date
    {
    short       year;                   /* Self-explanatory */
    short       day;
    short       hour;
    short       minute;
    float       second;
    };

typedef struct _skypos {
  double x,y,z;              /* local topocentric horizon coords 
			      * (x==east,y==north,z==up) */
  double alt,az;             /* ditto az measured E of N */
  double ha,dec,lat;         /* sky coords */
} skypos;

#include "observing_param.h"
#include "array_specification.h"

/* public function prototypes */
void rm_whitespace (char *line);
void clear_coord (struct coord *crd);
void msg (char *string,int level, ...);
int  compute_gha (struct observing_param *obsparam, double time, double *gha);
int  compute_uvw (struct stationspec *st1, struct stationspec *st2, 
		  double gha, struct observing_param *obsparam,
                  double *u_len, double *v_len, double *w_len);
void double_to_time (double time, struct date *sim_date);
int  account ( char *segment_name);
void environment();
double time_to_double(struct date);

extern FILE *msg_fp;

#endif
