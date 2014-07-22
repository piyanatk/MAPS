/************************************************************************/
/*                                                                      */
/* Defines the structures which get mapped into binary visibility       */
/* output files.  These files are read by the AIPS++ filler written by  */
/* ASTRON.                                                              */
/*                                                                      */
/* Created Jan 7 2002 by CJL                                            */
/*                                                                      */
/************************************************************************/
#ifndef VISIBILITY                      /* Handle multiple includes */
#define VISIBILITY 0

#include "utility.h"
#include "type_comp.h"

#define SIZ_FIELDNAME 32
#define SIZ_SIMNAME   36
#define SIZ_PRODNAME  8

struct main_header {
  char            simname[SIZ_SIMNAME];        /* Root name for all files in this sim */
  struct date     reftime;                     /* All times relative to this */
  struct coord    field_center;                /* RA and Dec of field center, J2000 */
  char            field_name[SIZ_FIELDNAME];   /* Arbitrary identifying string */
  char            stokes[SIZ_PRODNAME];        /* Stokes IQUV y/n (e.g. "ynnn") */
  int             ntime_block;                 /* To assist file reader */
};

struct freq_span 
{
    double          frequency;          /* Lower edge of channel in MHz */
    double          bandwidth;          /* Width of channel in MHz */
};

struct time_block_header
{
    float               start_time;     /* Seconds relative to reference time */
    float               duration;       /* Seconds */
    float               intg_time;      /* Seconds */
    int                 num_IFs;        /* Number of IFs */
    struct freq_span    ftable[1];      /* Variable length array, one for each IF. */
};


struct visgroup
{
    int             cchan;              /* Corr. chan no in this frq chan (currently REDUNDANT) */
    int             flag;               /* 0 = unflagged, nonzero meaning TBD */
    float           weight;             /* Applies to all stokes */
    complex         vis[4];             /* Variable length array, 1 per pol product */
};

struct visfreq
{
    int             chan;               /* Index into freq table (currently REDUNDANT) */
    int             ncchan;             /* Number of corr chans in this frq chan */
    struct visgroup visgroup[1];        /* Variable length array */
};

struct visblock
{
  double          u;                  /* meters */
  double          v;                  /* meters */
  double          w;                  /* meters */
  int             intg_no;            /* Integration interval number */
  int             station1;           /* Reference station number */
  int             station2;           /* Remote station number */
  int             nfreq;              /* Number of vis freqs in array */
                                      /* NOTE: this is currently redundant and the same as n_IF
                                         in the time_block_header */
  struct visfreq  visfreq[1];         /* Variable length array */
};


#endif

