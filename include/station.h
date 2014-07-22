/************************************************************************/
/*                                                                      */
/* Binary representation of array configuration information, for use in */
/* generating output data files to be read into AIPS++                  */
/* This is not the same as the array specification structures used      */
/* internally in the simulator                                          */
/*                                                                      */
/* Created Jan 7th 2002 by CJL                                          */
/*                                                                      */
/************************************************************************/
#ifndef STATION
#define STATION 0

struct ant_loc
    {
    double  north;              /* Meters w.r.t station center */
    double  east;               /* Meters w.r.t station center */
    double  height;             /* Meters w.r.t station center */
                                /* Future versions may have additional antenna */
                                /* characteristics here ... */
    };

struct station
    {
    int             station_no; /* This is the only station identifier in system */
    int             nant;       /* Number of antennas in this station */
    double          stn_x;      /* X coordinate of station center */
    double          stn_y;      /* Y coordinate of station center */
    double          stn_z;      /* Z coordinate of station center */
    struct ant_loc  ants[1];    /* Variable length array */
    };

#endif

