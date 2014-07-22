/************************************************************************/
/*                                                                      */
/* array_specifications.h                                               */
/*                                                                      */
/* Structures to store configuration information for a test LOFAR array */
/*                                                                      */
/* Created Jan 8th 2002 by CJL                                          */
/* Revised Feb 5th 2002 by CJL                                          */
/* Revised Jun 27th 2002 by CJL                                         */
/*                                                                      */
/************************************************************************/
#ifndef ARRAY_SPECIFICATION
#define ARRAY_SPECIFICATION  1

#define MAX_STATIONS 10000       /* Accommodate any future SKA numbers */
#define MAX_ARRAY_BORDER 1000   /* Max # verices of array convex hull  */
#define MAX_LNAME_SIZE 32

/* feed types */
#define NPOL   4
#define POL_X 'X' /* linear */
#define POL_Y 'Y'
#define POL_R 'R' /* circular */
#define POL_L 'L'
#define POL_I 'I' /* isotropic (ideal and ficticious) */

enum ant_type {
  ANT_ISOTROPIC,
  ANT_GAUSSIAN,
  ANT_SHORT_DIPOLE_ON_GROUNDPLANE,
  ANT_IDEAL_PARABOLOID,
  ANT_CROSSED_DIPOLES_ON_GROUNDPLANE,
  ANT_CROSSED_DIPOLES_HORIZONTAL,
  ANT_IDEAL_PARABOLOID_DUAL_LINEAR,
  ANT_IDEAL_PARABOLOID_DUAL_CIRCULAR
};


struct antenna {
  float           north;          /* Meters, wrt station position */
  float           east;           /* Meters, wrt station position */
  float           height;         /* Meters, wrt station position */
  enum ant_type   type;           /* one of the enumerated types above */
  short           num_receptors;  /* number of receptors (e.g. x,y or r,l 
				   * feeds) must be 1 or 2 */
  short           pol_type[2];    /* the pol type of the feed. 
				   * array the size of num_receptors */
  float           gain[2];        /* overall gain. array the size of 
				   * num_receptors  */
  float           phase[2];       /* overall phase. array the size 
				   * of num_receptors  */
  float           *params;        /* pointer to hold array of params 
				   * depending on the antenna type */
};

struct station_layout {
    char            layout_name[MAX_LNAME_SIZE];   /* Unique identifier */
    int             nant;
    struct antenna  *ants;
};

struct stationspec {
    double          north_offset;   /* Meters, nominal wrt ARRAY center */
    double          east_offset;    /* Meters, nominal wrt ARRAY center */
    double          height_offset;  /* Meters, nominal wrt ARRAY center */
    double          longitude;      /* Radians, east positive (derived) */
    double          latitude;       /* Radians (derived) */
    double          x_coord;        /* Meters, nominal wrt EARTH center 
				     * (derived) */
    double          y_coord;        /* Meters, nominal wrt EARTH center 
				     * (derived) */
    double          z_coord;        /* Meters, nominal wrt EARTH center 
				     * (derived) */
    struct station_layout *layout;
    double          low_el_limit;   /* Radians, lowest el station can observe.*/
    float           high_el_limit;  /* Radian. Highest el station can observe */
    float           sumweight;      /* Effective number of ants used in the 
				     * beamformer */
    float           sefd;           /* system equivalent flux density of 
				     * station (Jy) */
    int             st_id;
};

#include "utility.h"

struct array_specification {
  char                array_name[128];    /* Identifier */
  struct geocoord     location;           /* Geographical array location */
  int                 nst;                /* the number of stations */
  int                 num_pols;           /* the number of polarisation types 
					   * this array measures (1 or 2) */
  struct stationspec  station[MAX_STATIONS];
  int    vlbi; /* flag; added by L. Benkevitch, 10-Dec-2010 */
};

/* public function prototypes */
struct station_layout *get_stn_layout( char *layoutdir, char *layoutname);

#endif
