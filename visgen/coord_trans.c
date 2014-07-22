/* Functions to carry out coordinate transformations used in visgen 
 * for the ionospheric simulation component
 *   
 * Author: Justin C. Kasper (jkasper@cfa.harvard.edu)
 * 
 * History 
 * 
 *  04 March 2008 - Ported over from IDL code
 * 
 */


/* FUNCTION: coord_hor_sphcar
 * 
 * Convert between spherical and cartesian basis in the
 * horizon coordinate system.
 *
 * Cartesian coordinates (Note that this is _LEFT_ handed system
 * x - is in the plane horizontal to the surface of the Earth, pointing North
 * y - is in the plane horizontal to the surface of the Earth, pointing East
 * z - is normal to the surface of the Earth, pointing up
 *
 * Spherical coordinates
 * alt - altitude is angle out of x-y plane towards z
 * az  - azimuth is angle in x-y plane from x increasing towards y
 * r   - distance
 * 
 * 
 * Additional Inputs
 * 
 * togeographic: Set to 1 if converting from horizon to 
 *               geographic coordinates, set to zero if 
 *               converting from geographic to horizon
 * inradians: Set to 1 if all angle variables are in radians
 *            Set to 0 if all angle variables are in degrees
 * 
 */

#include <math.h>

static const double dtor=0.0174532925;

void coord_hor_sphcar(double *alt, double *az, double *r, 
		      double *x, double *y, double *z, 
		      int tocartesian, int inradians) {

  double altr, azr;

  if (tocartesian==1) {
    if (inradians==1) {
      altr = *alt;
      azr = *az;
    } else {
      altr = *alt * dtor;
      azr = *az * dtor;
    }

    *x = *r * cos(azr) * cos(altr);
    *y = *r * sin(azr) * cos(altr);
    *z = *r * sin(altr);
  } else {
    *r = sqrt((*x)*(*x) + (*y)*(*y) + (*z)*(*z) );
    altr = atan2(*z, sqrt((*x)*(*x) + (*y)*(*y)) );
    azr  = atan2(*y, *x);
    if (inradians==1) {
      *alt = altr;
      *az = azr;
    } else {
      *alt = altr / dtor;
      *az = azr / dtor;
    }
  }

}

/* FUNCTION: geo_sphcar
 * 
 * Convert between spherical and cartesian basis in the
 * geographic coordinate system.
 *
 *  Cartesian Version
 *  g1: x - From the center of the Earth through the prime meridian
 *  g2: y - In equatorial plane 90 degrees east of x
 *  g3: z - Along spin axis of the Earth pointing North
 *  
 *  Spherical Version
 *  g1: lat - latitude from the equator positive North
 *  g2: lon - angle from prime meridian positive East
 *  g3: r   - distance from center of the Earth
 * 
 * Additional Inputs
 * 
 * togeographic: Set to 1 if converting from horizon to 
 *               geographic coordinates, set to zero if 
 *               converting from geographic to horizon
 * inradians: Set to 1 if all angle variables are in radians
 *            Set to 0 if all angle variables are in degrees
 */


void coord_geo_sphcar(double *lat, double *lon, double *r, 
		      double *x, double *y, double *z, 
		      int tocartesian, int inradians) {

  double latr, lonr;

  if (tocartesian==1) {
    if (inradians==1) {
      latr = *lat;
      lonr = *lon;
    } else {
      latr = *lat * dtor;
      lonr = *lon * dtor;
    }

    *x = *r * cos(lonr) * cos(latr);
    *y = *r * sin(lonr) * cos(latr);
    *z = *r * sin(latr);
  } else {
    *r = sqrt((*x)*(*x) + (*y)*(*y) + (*z)*(*z) );
    latr = atan2(*z, sqrt((*x)*(*x) + (*y)*(*y)) );
    lonr  = atan2(*y, *x);
    if (inradians==1) {
      *lat = latr;
      *lon = lonr;
    } else {
      *lat = latr / dtor;
      *lon = lonr / dtor;
    }
  }

}

/* Convert from altazr 2 geographic coordinates */
//
// Input: alt, az, r - location in horizon coordinated
void coord_altazr2geo(double alt, double az, double r, 
		      double *x, double *y, double *z, 
		      double latr, double lonr, 
		      int tocartesian, int inradians) {
  
  const double re=6378.0;

  /* Location of center point in alt/az system*/
  double x0geo=0.0,y0geo=0.0,z0geo=0.0;

  /* Components of horizon unit vectors in geographic coordinates */
  double xxh, xyh, xzh, yxh, yyh, yzh, zxh, zyh, zzh;
  double xh, yh, zh;


  x0geo = re * cos(latr) * cos(lonr);
  y0geo = re * cos(latr) * sin(lonr);
  z0geo = re * sin(latr);
    
  /* Components of horizon unit vectors in geographic coordinates */
  xxh = -sin(latr)*cos(lonr);
  xyh = -sin(latr)*sin(lonr);
  xzh =  cos(latr);

  yxh = -sin(lonr);
  yyh =  cos(lonr);
  yzh =  yyh * 0.0;

  zxh =  cos(latr)*cos(lonr);
  zyh =  cos(latr)*sin(lonr);
  zzh =  sin(latr);

  /* x,y,z horizon system projection of r */
  xh = r*cos(az * dtor)*cos(alt * dtor);
  yh = r*sin(az * dtor)*cos(alt * dtor);
  zh = r*sin(alt * dtor);
  
  /* Locations in Geographic coordinates */
  *x = x0geo + xh * xxh + yh * yxh + zh * zxh;
  *y = y0geo + xh * xyh + yh * yyh + zh * zyh;
  *z = z0geo + xh * xzh + yh * yzh + zh * zzh;
  
}


/* convert between horizon and geographic coordinates
 *
 * Horizon Coordinates (h1, h2, h3)
 * 
 *  Cartesian Version
 *   h1: x - From station in local north direction
 *   h2: y - From station in local east direction
 *   h3: z - From station along radial from center of Earth
 * 
 *  Spherical Version
 *   h1: alt - Angle from horizon increasing towards zenith
 *   h2: az  - Angle from North towards East
 *   h3: r   - Physical distance from station along alt/az
 * 
 * Geographic Coordinates (g1, g2, g3)
 *
 *  Cartesian Version
 *  g1: x - From the center of the Earth through the prime meridian
 *  g2: y - In equatorial plane 90 degrees east of x
 *  g3: z - Along spin axis of the Earth pointing North
 *  
 *  Spherical Version
 *  g1: lat - latitude from the equator positive North
 *  g2: lon - angle from prime meridian positive East
 *  g3: r   - distance from center of the Earth
 * 
 * Additional Inputs
 * 
 * latr: Latitude of the station in geographic coordinates
 *       radians, zero at the equator and increasing North
 * lonr: Longitude of the station in geographic coordinates
 *       radians, zero at the prime medidian and increasing
 *       East
 * togeographic: Set to 1 if converting from horizon to 
 *               geographic coordinates, set to zero if 
 *               converting from geographic to horizon
 * inradians: Set to 1 if all angle variables are in radians
 *            Set to 0 if all angle variables are in degrees
 * spherical: Set to 1 for all variables to be in spherical coordinates
 *            Set to 0 for all variables to be in cartesian coordinates
 * units: Set to 0 if the physical length units are in km
 *        Set to 1 if the physical length units are Earth radius
 *        Set to 2 if the physical length units are m
 */

void coord_hor2geo(double *h1, double *h2, double *h3,
		   double *g1, double *g2, double *g3, 
		   double latr, double lonr, 
		   int togeographic, int inradians, int spherical,
		   int units) {

  double xh, yh, zh;
  double xxh, xyh, xzh, yxh, yyh, yzh, zxh, zyh, zzh;
  double x0geo, y0geo, z0geo;
  double xg, yg, zg;
  double dxg, dyg, dzg;

  double re;
  if (units==0) {
    re = 6378.0;
  }
  if (units==1) {
    re = 1.0;
  }
  if (units==2) {
    re = 6378000.0;
  }
  
  x0geo = re * cos(latr) * cos(lonr);
  y0geo = re * cos(latr) * sin(lonr);
  z0geo = re * sin(latr);

  /* Check if we are converting to geographic coordinates */
  if (togeographic==1) {

    /* Check if input is in spherical coordinates */
    if (spherical==1) {
      coord_hor_sphcar(h1, h2, h3, &xh, &yh, &zh, 1, inradians);
    } else {
      xh = *h1;
      yh = *h2;
      zh = *h3;
    }

    /* Now convert to geographic coordinates */
    /* Components of horizon unit vectors in geographic coordinates */
    xxh = -sin(latr)*cos(lonr);
    xyh = -sin(latr)*sin(lonr);
    xzh =  cos(latr);

    yxh = -sin(lonr);
    yyh =  cos(lonr);
    yzh =  yyh * 0.0;

    zxh =  cos(latr)*cos(lonr);
    zyh =  cos(latr)*sin(lonr);
    zzh =  sin(latr);


    /* Locations in Geographic coordinates */
    xg = x0geo + xh * xxh + yh * yxh + zh * zxh;
    yg = y0geo + xh * xyh + yh * yyh + zh * zyh;
    zg = z0geo + xh * xzh + yh * yzh + zh * zzh;

    /* Convert to spherical coordinates if needed */
    if (spherical==1) {
      coord_geo_sphcar(g1, g2, g3, &xg, &yg, &zg, 0, inradians);
    } else {
      *g1 = xg;
      *g2 = yg;
      *g3 = zg;
    }

  } else {
    /*  Convert from geographic to horizon coordinates */

    /* If geographic coordinates are spherical convert to cartesian */
    if (spherical==1) {
      coord_geo_sphcar(g1, g2, g3, &xh, &yh, &zh, 0, inradians);
    } else {
      xh = *g1;
      yh = *g2;
      zh = *g3;
    }
    
    /* Now convert to horizon coordinates */
    
    /* Components of horizon unit vectors in geographic coordinates */
    xxh = -sin(latr)*cos(lonr);
    xyh = -sin(latr)*sin(lonr);
    xzh =  cos(latr);

    yxh = -sin(lonr);
    yyh =  cos(lonr);
    yzh =  yyh * 0.0;

    zxh =  cos(latr)*cos(lonr);
    zyh =  cos(latr)*sin(lonr);
    zzh =  sin(latr);

    /*  Difference between positions and geographic location of array
	; in geographic coordinates */
    dxg = xg - x0geo;
    dyg = yg - y0geo;
    dzg = zg - z0geo;
    
    /* Project delta in geo onto each component of horizon system */
    xh = dxg*xxh + dyg*xyh + dzg*xzh;
    yh = dxg*yxh + dyg*yyh + dzg*yzh;
    zh = dxg*zxh + dyg*zyh + dzg*zzh;

    /* Convert to spherical coordinates if needed */
    if (spherical==1) {
      coord_hor_sphcar(h1, h2, h3, &xh, &yh, &zh, 0, inradians);
    } else {
      *h1 = xh;
      *h2 = yh;
      *h3 = zh;
    }

  }

}


/* simple transform for equatorial to horizon coords */
/* all angles in radian. This is supposed to be identical to the sla_e2h function */
void simple_E2H(const double ha,const double dec,const double lat,double *az,double *el) {
    double sh,sd,ch,cd,sl,cl;
    double x,y,z,r,a;
    
    sh = sin(ha); ch = cos(ha);
    sd = sin(dec);cd = cos(dec);
    sl = sin(lat);cl = cos(lat);
    
    /* transform to cartesian on unit sphere */
    x = -ch*cd*sl + sd*cl;
    y = -sh*cd;
    z =  ch*cd*cl + sd*sl;

    /* transform to spherical coords */
    r = sqrt(x*x + y*y);
    if (r==0.0) {
        a = 0.0;
    }
    else {
        a=atan2(y,x);
    }
    /* make az betweeen 0 and 2pi */
    if (a < 0.0) a += 2.0*M_PI;
    *az = a;
    *el= atan2(z,r);
}


/* simple conversion of spherical to direction cosines. This is supposed to be
   idential to the sla_cs2c function. Angles in radian. */
void sperical2DirectionCosines(double lon, double lat, double vec[3]){
    double clat;
    
    clat = cos(lat);
    vec[0] = cos(lon)*clat;
    vec[1] = sin(lon)*clat;
    vec[2] = sin(lat);
}


