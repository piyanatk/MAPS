/************************************************************************/
/*                                                                      */
/* Created July 2007 by Randall Wayth                                   */
/*                                                                      */
/************************************************************************/
#ifndef OOB_H
#define OOB_H

#include "type_comp.h"
#include "observing_param.h"

typedef struct _source_model {
  int type;
  double ra,dec;   /* sky location (decimal hours, decimal degrees) */
  float flux[4];  /* Stokes params I,Q,U,V */
} source_model_t;

/* oob station beams. One of these structs for each station. Each struct holds the beams for all the sources */
typedef struct _oob_beam {
  complex *response[4]; /* one Jones matrix element for each source px,py,qx,qy */
  int valid;
  int valid_phase;
  double  *phase;       /* ionospheric phase for each source, per antenna */
} oob_beam_t;

/* external function prototypes */
void init_oob(const int n_oob, const int n_stations, const int n_pol_prodcuts, oob_beam_t **oob_beams);
int read_oob(char *filename, source_model_t **ooblist, int *noob);
void oob_set_debug(int level);
int oob_integrate(struct observing_param *obsparam, struct array_specification *arrayspec,double gha,double freq,
                  int st1_index,int st2_index,source_model_t *oobs, int num_src, oob_beam_t *beams, complex vis[4]);

#endif
