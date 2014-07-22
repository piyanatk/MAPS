/* public function prototype */
/*
 * 06 Jan 2011 added vnoise parameter to add_noise() (Jeremy Steeger) 
 */

#ifndef VISGEN_ADDNOISE_H
#define VISGEN_ADDNOISE_H

int seed_noise (int seed, int myid);
void addnoise_set_debug(int level);

int add_noise (// st1spec, st2spec, freq, obsparam, visibility)
	       const struct stationspec *st1spec,
	       const struct stationspec *st2spec,
	       const struct stnbeam     *st1beam,
	       const struct stnbeam     *st2beam,
	       const float   bandwidth_Hz,          // in Hz
	       const float   integration_time,      // in seconds.
	       struct uvgrid_parms *uvgrid,
	       const float   global_sefd,           // in Jy
	       const float   system_efficiency,     // dimensionless [0..1]
	       complex *visibility,
	       int n_pols,
	       float *vnoise);

#endif
