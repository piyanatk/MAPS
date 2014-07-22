#ifndef MAKEBEAM_H
#define MAKEBEAM_H

#include "station_beam.h"

int makebeam(struct stationspec *stn, struct stnbeam *beam1, 
	     struct beamgrid *bgrid, double freq, double gha,
	     struct observing_param *obsparam);
int station_response(struct stationspec *stn, skypos *ref_dir, 
		     double hang, double dec, double freq, complex beam[4]);
void MakeSkyPos(double ha, double dec, double lat, skypos *result);
void makebeam_set_debug(int debug_level);
#endif
