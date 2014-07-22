/* public function prototype */

#ifndef COORD_TRANS_H
#define COORD_TRANS_H

void coord_hor_sphcar(double *alt, double *az, double *r, 
		      double *x, double *y, double *z, 
		      int tocartesian, int inradians);
void coord_geo_sphcar(double *lat, double *lon, double *r, 
		      double *x, double *y, double *z, 
		      int tocartesian, int inradians);
void coord_altazr2geo(double alt, double az, double r, 
		      double *x, double *y, double *z, 
		      double latr, double lonr, 
		      int tocartesian, int inradians);
void coord_hor2geo(double *h1, double *h2, double *h3,
		   double *g1, double *g2, double *g3, 
		   double latr, double lonr, 
		   int togeographic, int inradians, int spherical,
		   int units);
void simple_E2H(const double ha,const double dec,const double lat,double *az,double *el);
void sperical2DirectionCosines(double lon, double lat, double vec[3]);
#endif


