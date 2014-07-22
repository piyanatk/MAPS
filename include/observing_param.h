/************************************************************************/
/*                                                                      */
/* Information required to initiate visibility generation in make_vis   */
/*                                                                      */
/* Created Jan 8th 2002 by CJL                                          */
/*                                                                      */
/************************************************************************/
#ifndef OBSERVING_PARAM_H
#define OBSERVING_PARAM_H 1

#include "utility.h"

#define MAXSCAN 100   /* Max scans (not correlator integrations!) */
#define MAXFREQ 20    /* Max channels (not correlator channels!) */

struct freqinfo {
  double          frequency;     /* Frequency at bottom edge in MHz */
  double          bandwidth;     /* Channel bandwidth in MHz */
};

struct scaninfo {
  int             gha_used;      /* flag to indicate the scan starts in GHA,
				  * not absolute time. */
  double          gha_start;     /* starting GHA (radian) 
				  * if flag above is set */
  struct date     start;         /* Start time of scan, UT */
  double          duration;      /* Length of scan, seconds */
  int             nfreq;         /* Number of IFs in scan */
  struct freqinfo freq[MAXFREQ]; /* Frequency table for this scan */
};


struct observing_param {
  double          point_cent_RA;  /* Field of view (antenna pointing) 
				   * center coordinates, radian */
  double          point_cent_DEC; 
  double          phase_cent_RA;  /* Field of view (correlator) phase
				   * centre coords, radian */
  double          phase_cent_DEC;
  double          fov_RA;         /* FOV size, arcsec */
  double          fov_dec;
  double          integ_time;     /* Correlator integration time in seconds */
  double          corr_chan_bw;   /* Correlator channel bandwidth MHz */
  //  double          el_limit;     /* Elevation limit in radians - moved to
  //                                 * array_specification per station */
  int             nscan;          /* Number of scans in this observation */
  int             nt_cells;       /* Number of time cells per integ_time */
  int             nf_cells;       /* "    "  freq  "     "  corr_chan_bw */
  struct scaninfo scan[MAXSCAN];  /* Array of scan information */
};

#endif
