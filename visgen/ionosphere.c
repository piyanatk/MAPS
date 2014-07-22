/************************************************************************
    Functions associated with computing the TEC and phase from an
    ionosphere model. Collated and updated from functions
    compute_ionoscreen.c,open_ionosphere.c, path.c, interp_t.c
    
    Modification History
    01/01/2008 - Randall Wayth - Original subroutines gathered 
    Previous work - Roger Cappallo, Shep Doeleman, Colin Lonsdale
    
    This module calculates the phase shift due to the ionosphere given
    an ionosphere model, a station location and a look direction. The
    phase shift is calculated as a difference in excess path length
    through the ionosphere between different stations. It does *not*
    actually calculate the *real* refractive offset of a source, just
    the excess path length. This means that any shift seen in the
    location of a source in a map made from interferometric data is
    purely due to a phase shift as seen by the interferometer, rather
    than the true movement of a source due to real
    refraction. Furthermore, the station gains applied to calculating
    the apparent brightness/flux density of a source are those of the
    true direction to the source, not where it ends up in the image.
    
    The phase shift, phi, on a baseline is simply phi =
    2*pi*(delta_L)/lambda = 2*pi*f/c*(delta_L) where delta_L is the
    difference in excess path length between stations in the baseline
    and lambda is the observing wavelength (or f is the frequency).
    Following chapter 13 of Thompson, Moran & Swenson, equation
    13.138, the excess path length is 40.3/(f^2) integral{n_e(h) dh}
    where n_e is the electron density (e m^-3) at height h.  The
    integral is simply the column density or "total electron content"
    (tec) of the ionosphere, as seen by a station in a particular look
    direction.  So, phi = -40.3*2pi/(cf)*[tec_station2 - tec_station1]
    = -8.44e-7/f*[tec_station2 - tec_station1] Note that various
    assumptions go into generating equation 13.138, specifcically that
    the observing freq is >> the plasma freq of the ionosphere and
    gyrorotation freq of electrons in Earth's B field.
    
    TESTING: For an input ionosphere which is a phase ramp across the
    array with gradient=m [e/m^3/km], then for a baseline of length b
    in the direction of the gradient, the delta_TEC will be = b*m*h
    where h is the length through the ionosphere Now equate this to
    the expected phase difference on the same baseline = b*sin(theta)
    (= b*theta for small angles) for a source angular offset of theta
    radian. Thus theta = (-8.44e-7/freq)*m*h

************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <time.h>
#include "fitsio.h"
#include "netcdf.h"
#include "type_comp.h"
#include "comp_func.h"
#include "array_specification.h"
#include "observing_param.h"
#include "convolve.h"
#include "station_beam.h"
#include "visgen.h"
#include "utility.h"
#include "coord_trans.h"
#include "ionosphere.h"
#include "fourc.h"

/* private function prototypes */
int ion_density (struct ion_model *isphere,
                 const double x[3],
                 double *value);
void *matloc(int dim[], size_t cellsize);
void matfree(void *matt, int dim[]);
int tcube_init_path_str(struct path_str *path, double x0, double y0, double z0,
                        double nx, double ny, double nz,
                        double step, int steps);
void tcube_init_tcube_str(struct tcube_str *tcube, double s,
                          int nx, int ny, int nz,
                          double i_scale, double o_scale, double rms);
void tcube_fille_simple(struct tcube_str *tcube );
void tcube_fill_turbulence(struct tcube_str *tcube );
int tcube_path( struct tcube_str *tcube,
                struct path_str *path,
                double *density);
void tcube_print(struct tcube_str *tcube );
double tcube_calc_average( struct tcube_str *tcube );
double tcube_calc_rms( struct tcube_str *tcube );
double tcube_trilinear(struct tcube_str *tcube,
                      double x, double y, double z);
void tcube_write_fits(struct tcube_str *tcube, char *outfilename);
int TranslateCoord(int coord, int dim_length);
int TranslateCoord2(int coord, int dim_length);
double envelope_function(double k, double o_scale, double i_scale);
void translate_array(int x_dim, int y_dim, int z_dim, complex **temp_iono);

/* private global variables */
static int debug=0;

void ionosphere_set_debug(int level) {
    debug = level;
}

int compute_ionoscreen (struct stationspec *stn,
                        int                 st_index,
                        struct beamgrid    *bgrid,
                        double             gha,
                        double             fovra,
                        double             padfact,
                        struct ion_model   *isphere,
                        int                 center_ion,
                        struct iono_screen *ionoscreen)
{
    int i, j,ret,index,fov_range;
    double ra, dec, dra, az, el,lha,phase,xs[3],enh[3];

    // tcube check
    if(debug && isphere->tcube != NULL) {
        printf(" TCUBE dimensions ionoscree X,dx: %d,%g Y,dy: %d,%g Z,dz: %d,%g\n",isphere->tcube->nx,
	    isphere->tcube->sx, isphere->tcube->ny, isphere->tcube->sy, isphere->tcube->nz, isphere->tcube->sz);
    }
    

    // Check to make sure that size>padfactor
    if (ionoscreen->size_x <= padfact || ionoscreen->size_y <= padfact)  {
      msg ("Ionoscreen size %dx%d must exceed padfactor %lf",3,
	   ionoscreen->size_x,ionoscreen->size_y,padfact);
      return (1);
    }
    // precompute station-based latitude and local
    // hour angle
    lha = gha + stn->longitude;
    
    // loop over 2-D output grid. i==Y, j==X.
    for (i=0;i<ionoscreen->size_y;i++) {    
        for (j=0;j<ionoscreen->size_x;j++) {

            // Get index of this beam pixel in 1-D complex array
            index = i * ionoscreen->size_x + j;

            // Only compute path integral over field of view
            fov_range = (ionoscreen->size_x)/(2*padfact);
	

            if ((abs(i-ionoscreen->size_y/2)<=fov_range) && (abs(j-ionoscreen->size_x/2)<=fov_range) &&
                (bgrid->grid[index].ra!=BAD_POINT && bgrid->grid[index].dec!=BAD_POINT) ) {		    
                // Read RA and Dec of this beam sample pixel
                ra  = bgrid->grid[index].ra;
                dra = ra - fovra;       // difference from center of fov
                dec = bgrid->grid[index].dec;
                if (debug) printf("comp_iono: x,y: %d,%d\tra %lf dra %lf dec %lf\n", j,i,ra, dra, dec);
	  
                // Calculate the location of the station in the ENH coordinate system of the 
                // center of the array
                enh[0] = 1.0e-3 * stn->east_offset;
                enh[1] = 1.0e-3 * stn->north_offset;
                enh[2] = 1.0e-3 * stn->height_offset;
	  
                // find az & el towards the source in coordinate system at center of array
                // FIXME: Note that this is only valid for small arrays. This whole section needs to be re-done
                // for VLBI type arrays where there is a different ionosphere over each station.
                simple_E2H(gha + isphere->center_lon - dra,dec,isphere->center_lat,&az,&el);

                // find stn coords wrt ionosphere
                if (center_ion) {               // offset relative to ion. center
                    xs[0] = 1.0e-3 * stn->east_offset  + isphere->origin[Y] + 
                            (isphere->delta[Y] * (isphere->dim[Y] - 1)) / 2;
                    xs[1] = 1.0e-3 * stn->north_offset + isphere->origin[X] + 
                            (isphere->delta[X] * (isphere->dim[X] - 1)) / 2;
                    xs[2] = 0.0;
                    if(debug) printf ("centered coords: %g %g\n", xs[0], xs[1]);
                }
                else {               
                    // use absolute coords. of station and ion.
                    //FIXME: this is wrong. Need to calculate the offset in linear "East" and "North" coords
                    // with respect to the centre of the ionosphere and array. To do that 
                    // propertly requires some algebra
                    // similar to that in the path() function.
                    xs[0] = stn->longitude;
                    xs[1] = stn->latitude;
                    xs[2] = 0.0;
                    if(debug) printf ("absolute coords: %g %g\n", xs[0], xs[1]);
                }
	  
                // find phase delay for this station via path integral through ionosphere model.
                ret = path(isphere, enh, az, el, &phase);

                if (ret) return (ret);
	  
                // phase at 1 MHz, noting that the path integral was in km. (see Thompson et al. eqn. 13.138)
                ionoscreen -> phase_delay[st_index][index] = -8.438e-10 * phase;
            }
	        else ionoscreen -> phase_delay[st_index][index] = 0.0;
	
	        if(debug) printf("st: %d, index: %d, grid: %d,%d, phase_delay: %8.3f\n",
			        st_index,index, j,i,ionoscreen -> phase_delay[st_index][index]);
        }
    }
    ionoscreen->valid[st_index] = 1;
    return (0);
}



//   path - does line integral through the ionosphere
//   curved earth version        rjc 2003.4.10
//   
//   Inputs:
//     isphere  pointer to structure with complete ionosphere model
//     xs       station coords [E long(km), Lat (km), alt (km)]
//     az       azimuth of start of path integral (rad) from station
//     el       elevation "  "    "   "     "       "
//     
//   Outputs:
//     integral integrated electron density (radians of phase at 1 MHz)
//     
//   Function returns:
//     0  nominal return; non-zero on errors
//
/* in the curved Earth version, we interpret the input ionosphere density cube as being
   drooped over the surface of the earth. So x,y,z in the data cube becomes east, north, 
   height (represented by
   xi,eta,zeta). Note that increments in (xi,eta,zeta) are still "linear"- latitude plays no role.
   They are "topocentric spherical coords".

   (x,y,z) is a cartesian coord system centred on the surface of the earth.
   we compute a linear path integral through x,y,z space looking out from the station location.
   For each point on this integral, x,y,z must be transformed into xi,eta,zeta to sample
   the ionosphere density.
   
   If the origin of (x,y,z) marks the point on the earth below the centre of 
   the ionosphere and R is the radius
   of the earth, then at any point in (x,y,z), (R+z)^2 + (x^2+y^2) = (R+zeta)^2, so
   zeta = sqrt((R+z)^2 + x^2 + y^2) - R
   the length of an arc on the earth's surface from the (x,y,z) origin is: rho=R*phi
   where tan(phi) = sqrt(x^2+x^2)/(R+z), so
   rho = R*atan(sqrt(x^2+x^2)/(R+z))
   By the definition of coord system, azimuthal angles are preserved, so
   sin(az) = x/sqrt(x^2+y^2)
   cos(az) = y/sqrt(x^2+y^2)
   so
   xi = rho*sin(az) + xi_station
   eta= rho*cos(az) + eta_station
   where xi_station and eta_station are the location of the station in xi,eta coords.
   
   The limits of the path integral must also be computed. in (x,y,z) 
   the intersection of the ionosphere
   at height zeta for elevation el is
   z = -Rsin^2(el) + Rsin^2(el)*sqrt(1+[zeta^2+2Rzeta/R^2sin^2(el)])
*/
int path (struct ion_model *isphere,
	  double enh[3], 
          double az,
          double el,
          double *integral)
{
    int i,
        nstep,
        rc=0;

    double r, r1, r2,           
      z1, z2,                  // integration height limits in local cartesian coords
      zeta1, zeta2,            // integration height limits in curved ionosphere coords
      xi[3],                   // integration point in ionospheric coords.
      h,                       // integration step size
      coeff,                   // integration weight coeff for numerical integration
      sum,
      density=0.0,
      sel2,                    // sin^2(el)
      cel2;                    // cos^2(el)
    
    double p_geo_x, p_geo_y, p_geo_z, p_geo_lat, p_geo_lon, p_geo_alt, p_geo_r;
    double p_arr_hor_x, p_arr_hor_y, p_arr_hor_z;
    double amin,amax;               // Minimum and maximum altitudes in global ionopshere model [km]
    double rmin,rmax;               // Minimum and maximum distances from station for integration

    //   struct tcube_str *tcube;
    double dn_turb=0.0;                 // density perturbation as fraction due to turbulence
    
    const double RE = 6378.0;        // Earth radius [km].
    const double dtor = M_PI / 180.0;

    // Boole's rule for integration - see A&S 25.4.14 (Note typo in A&S calls it "Bode's" rule.)
    // this is just a 4th degree Newton-Cotes formula. There is a 5th element in c, but it is also 7.0
    // and since we conceptually break the integral up into lots of small integrals, we don't count
    // the multiple 7.0 factors twice.
    double c[4] = {7., 32., 12., 32.};

    // debug - force source directly overhead at all sites
    // el = 90. / (180/M_PI); az = 270. / (180/M_PI);

    if(debug) {
       printf("path: az %f el %f\n", az, el);
       fprintf(stdout, " isphere: center_lat: %g center_lon: %g zlog: %d\n", 
    	    isphere->center_lat/dtor, isphere->center_lon/dtor, isphere->zlog);
    }

    sel2 = sin(el);                // pre-compute sin (elevation) squared
    sel2 *= sel2;
    cel2 = cos(el);                // pre-compute cos (elevation) squared
    cel2 *= cel2;


    // Make sure source is above the horizon
    if (el <= 0.0) {
        printf ("Source below horizon at elevation %f deg; path aborting\n", el * (180.0/M_PI));
        return (1);
    }

    if (isphere->do_turbulence == 1 && debug) {
      fprintf(stdout, "Using turbulence in path, extracting tcube\n");
      // tcube = new_isphere.tcube;
      fprintf(stdout, " TCUBE dimensions in path X,dx: %d,%g Y,dy: %d,%g Z,dz: %d,%g\n",  isphere->tcube->nx,
	      isphere->tcube->sx, isphere->tcube->ny, isphere->tcube->sy, isphere->tcube->nz, isphere->tcube->sz);
      
    }

    // calculate lower and upper bounds of integral in ionosphere coords
    // note calculation of upper bound depends on log/linear sampling in z of model
    zeta1 = isphere->origin[Z];
    amin = isphere->origin[Z];
    if (isphere->zlog == 1) {
      zeta2 = pow(10.0, log10(isphere->origin[Z]) + isphere->delta[Z] * (isphere->dim[Z] - 1));
      amax = pow(10.0, log10(isphere->origin[Z]) + isphere->delta[Z] * (isphere->dim[Z] - 1));
    } else {
      zeta2 = isphere->origin[Z] + isphere->delta[Z] * (isphere->dim[Z] - 1);
      amax = isphere->origin[Z] + isphere->delta[Z] * (isphere->dim[Z] - 1);
    }

    // Translate minimum and maximum altitudes (amin,amx)  in global model into 
    // minimum and maximum distances (rmin,rmax) of line of sight integral
    rmin = sqrt(RE*RE*cel2 + amin*amin + RE*amin) - RE*cos(el);
    rmax = sqrt(RE*RE*cel2 + amax*amax + RE*amax) - RE*cos(el);


    // convert lower/upper bounds in height (zeta) to (x,y,z) cartesian coords (actually, just z)
    z1 = -RE * sel2 + RE * sel2 * sqrt (1 + (zeta1 * zeta1 + 2 * RE * zeta1) / (RE * RE * sel2));
    z2 = -RE * sel2 + RE * sel2 * sqrt (1 + (zeta2 * zeta2 + 2 * RE * zeta2) / (RE * RE * sel2));
    // z1 = zeta1; z2 = zeta2;      // debug - force flat earth

    // now convert the bounds on height to bounds on the path integral distance.
    // we are going to integrate from r1 to r2 in (x,y,z) space.
    r1 = z1 / sin (el);
    r2 = z2 / sin (el);

    // printf ("z1 %12.6lf z2 %12.6lf r1 %12.6lf r2 %12.6lf\n", z1, z2, r1, r2);
    
    // determine number of steps and stepsize.
    // generally, we assume that the ionospheric density is relatively slowly varying between
    // pixels. At the very least it is Nyquist sampled, but for any interpolations to work
    // correctly it should be oversampled by a factor of ~3 anyway. So set the number of steps
    // to be the number of samples in height, or at least 4, since a 4th degree integration is used below.
    // (which requires 5 samples). Using the 4th degree means nstep must be a multiple of 4.
    // note also that if the sampling is log in height then we want resolution set at 
    // the lowest altitudes
    //
    // Also, if you are using turbulence, over-ride this algorithm and use the spatial resolution of
    // the density

    if (isphere->do_turbulence == 1) {
      nstep = (rmax-rmin)/(isphere->tcube->sx);
    } else {
      if (isphere->zlog == 1) {
	nstep = (rmax-rmin)/(r1*pow(10.0,isphere->delta[Z])-r1);
      } else {
	nstep = isphere->dim[Z];  // should be at least 1
      }
    }

    while(nstep%4 != 0) nstep++;
    if (debug) fprintf(stdout, "Number of steps: %d  %lf\n", nstep, isphere->tcube->sx);

    h = (rmax - rmin) / nstep;

    // original version
    // h = (r2 - r1) / nstep;
    
    for (i=0; i<4; i++) c[i] *= h*(2./45.);       // adjust integration coeffs. for stepsize

    r = r1; // initialize to the start of the path integral.
    sum = 0.0;
    for (i=0; i<=nstep; i++) {
      
      
      
      // calc x,y,z along path integral in cartesian horizon system centered on array
      
      p_arr_hor_x = r * cos (az) * cos (el) + enh[1];  // north component
      p_arr_hor_y = r * sin (az) * cos (el) + enh[0];  // east component
      p_arr_hor_z = r * sin (el) + enh[2];             // height component
      
      // fprintf(stderr, " East: %g North: %g Height: %g\n", enh[0], enh[1], enh[2]);

      // fprintf(stderr, " p_arr_hor_x: %g p_arr_hor_y: %g p_arr_hor_z: %g\n", 
      // p_arr_hor_x, p_arr_hor_y, p_arr_hor_z);

      // convert from horizon system to geographic (input should be N, E, H)
      coord_hor2geo(&p_arr_hor_x, &p_arr_hor_y, &p_arr_hor_z, 
		    &p_geo_x, &p_geo_y, &p_geo_z, 
		    isphere->center_lat, isphere->center_lon, 
		    1, 1, 0, 0);

      // fprintf(stderr, " p_geo_x: %g p_geo_y: %g p_geo_z: %g\n", 
      // p_geo_x, p_geo_y, p_geo_z);

      // convert from cartesian geographic to spherical geographic
      coord_geo_sphcar(&p_geo_lat, &p_geo_lon, &p_geo_r, 
		       &p_geo_x, &p_geo_y, &p_geo_z, 
		       0, 1);
      p_geo_alt = p_geo_r - RE;

      if (debug) fprintf(stdout, " p_geo_lat: %g p_geo_lon: %g p_geo_alt: %g\n", p_geo_lat/dtor, p_geo_lon/dtor, p_geo_alt);


      xi[0] = p_geo_lon; 
      xi[1] = p_geo_lat; 
      xi[2] = p_geo_alt; 

      // if (debug) printf ("r: %g, x: %g %g %g, rho %g xy12 %g\n", r, x[0], x[1], x[2],rho,xy12);
      
      coeff = c[i % 4];
      
      // all interior edge points are in two sums, but aren't explicitly counted twice,
      // so make sure they are weighted double
      if (i%4==0 && i!=0 && i!=nstep) coeff *= 2;
      
      // is simple model enabled then do it else run global/turbulence/tid
      if (isphere->do_simple == 1) {
	if ((p_arr_hor_z > isphere->simple_zmin) && (p_arr_hor_z < isphere->simple_zmax)) {
	  density = isphere->simple_dc + isphere->simple_grad*p_arr_hor_x;
	} else 
	  density = 0.0;
      } else {
	// find interpolated ion density at this location 
	if ((p_geo_alt < amin) || (p_geo_alt > amax)) {
	  density = 0.0;
	} else {
	  rc = ion_density (isphere, xi, &density);
	}
	
	// If turbulence is enabled then modify the density value
	if (isphere->do_turbulence == 1) {
	  dn_turb = tcube_trilinear(isphere->tcube, 
				    p_arr_hor_x + isphere->t * isphere->tcube->vx, 
				    p_arr_hor_y + isphere->t * isphere->tcube->vy, 
				    p_arr_hor_z + isphere->t * isphere->tcube->vz);
	  if(debug) fprintf(stdout, "   Turbulence:  %g  %g  %g\n", density, dn_turb, density*(1.0 + dn_turb) );
	  density = density * (1.0 + dn_turb);
	  
	}
	
	
	// If there is a TID enabled and we are in the right altitude range then add it
	if ((isphere->do_tid == 1) && (p_geo_alt > isphere->tid_amin) && (p_geo_alt < isphere->tid_amax)) {
	  density = density * (1.0 + isphere->tid_amp * 
			       sin( 2.0*M_PI*(isphere->tid_vel * isphere->t + p_arr_hor_x) / (isphere->tid_len) ) );
	}
	
	
	// if (debug) printf("step: %d. xi %g %g %g,\tdensity %g\n", i, 
	// xi[0]/dtor, xi[1]/dtor, xi[2], density);
	if (rc) {
	  printf ("ion_density ran into problems for az %lf el %lf , return code %d\n", az, el, rc);
	  return (1);
	}
      
      }
      
      // and add into the integral
      sum += coeff * density;
      
      r += h;                     // update the independent variable
    }
    
    if(debug) printf("sum %g\n", sum);
    *integral = sum;  // return the TEC.
    return (0);
}


/**********************************************************************
 *  interp_t  -  interpolates a function of 4 independent variables   *
 *               to a specific time, yielding a function of the       *
 *               three spatial variables (x, y, z)                    *
 *  Note that interp_t attempts to read the minimum of new data, in   *
 *  order to be as efficient as possible. At entry, the ionosphere    *
 *  model structure is assumed to be pointing at 4 contiguous time    *
 *  slices.                                                           *
 *                                                                    *
 * Inputs:                                                            *
 *   isphere - pointer to an ion_model structure                      *
 *   t       - time at which interpolation is to be done              *
 * Outputs:                                                           *
 *   isphere->slice - points to output 3D field                       *
 *   isphere->ntslice - indices of slices in each of 4 RAM slots      *
 *   isphere->tslice - filled in as necessary to create 4 time points *
 *                     surrounding the requested time                 *
 *   isphere->t - time at which slice is evaluated                    *
 *                                                                    *
 * Returns:                                                           *
 *                                                                    *
 *   interp_t (returned value of function)                            *
 *      0: no error                                                   *
 *      1: out of bounds error (x outside of grid, or too close       *
 *         to its edge)                                               *
 * First created                                rjc  2002.11.22       *
 *                                                                    *
 *********************************************************************/

int interp_t(struct ion_model *isphere, double t) {
  int index,
    offset,
    i,
    it,
    ix,
    iy,
    iz,
    ibeg,
    iend,ret;
  size_t start[4], count[4];    
  float p,c[4],*temp[4];
  long inc[4] = { 1, 1, 1, 1 };
  long fpixel[4]={1,1,1,1};
  long lpixel[4]={1,1,1,1};
  double nulval=0.0;
  int status=0,anynul=0;

  // Save new interpolation time into isphere
  isphere->t = t;
  msg ("Interpolating global ionosphere to %lf seconds UT", 2, isphere->t);

  /* preserve tcube
     fprintf(stdout, " TCUBE dimensions in interp_t X,dx: %d,%g Y,dy: %d,%g Z,dz: %d,%g\n", isphere->tcube->nx,
     isphere->tcube->sx, isphere->tcube->ny, isphere->tcube->sy, isphere->tcube->nz, isphere->tcube->sz);
  */

  //struct tcube_str *tcube;
  //tcube = isphere->tcube;
  
  // initialize variables for later reads
  for (i=X; i<=Z; i++) {
    start[i] = 0;
    count[i] = isphere->dim[i];
  }
  
  count[T] = 1;
  
                                    // first figure out if the desired 4
                                    // time slices are in RAM and which new
                                    // slices need to be read
  index = (t - isphere->origin[T]) / isphere->delta[T];
  msg("Iono times: start: %g, delta: %g, ntimes: %d",-2,
      isphere->origin[T],isphere->delta[T],(int)isphere->dim[T]);

  // check that desired entries are in file
  if (index < 1 || index > isphere->dim[T] - 2) {
    msg("Time %g index %d not within interpolation range",2,t, index);
    return (1);
  }

  offset = index - isphere->ntslice - 1;

  if (offset < -3 || offset > 3)  {
    ibeg = 0;                   // total mismatch, need to read all 4 slices
    iend = 3;
  }
  else if (offset < 0)            // -1..-3
    {
      ibeg = 0;                   // partial mismatch, need to read a subset
      iend = - offset - 1;        // 2, 1, or 0
      
      for (i=0; i<4; i++)         // cyclical shift of slice pointers
        temp[i] = isphere->tslice[i];
      
      for (i=0; i<4; i++)
        isphere->tslice[i] = temp[(i+offset+4)%4];
    }
  else if (offset == 0)            // 0 offset   - do nothing
    {
      ibeg = 0;
      iend = -1;
    }
  else                            //  1..3
    {
      ibeg = 4 - offset;          // partial mismatch, need to read a subset
      iend = 3;
      
      for (i=0; i<4; i++)         // cyclical shift of slice pointers
        temp[i] = isphere->tslice[i];
      
      for (i=0; i<4; i++)
        isphere->tslice[i] = temp[(i+offset)%4];
    }
  // printf ("index %d ntslice %d offset %d ibeg %d iend %d\n",
  //         index, isphere->ntslice, offset, ibeg, iend);
  
  // update index of 1st data record in RAM
  isphere->ntslice = index - 1;
  // read in new slices if necessary
  for (i=ibeg; i<=iend; i++) {
      start[T] = isphere->ntslice + i;

      if(isphere->global_type == 0) {
	// Read from a netCDF file //
	ret = nc_get_vara_float (isphere->file_id, isphere->grid_id, start, count, isphere->tslice[i]);
	if (ret != NC_NOERR) {
	  printf ("Error %d reading time slice %d\n", ret, i);
	  return (2);
	}

      } else {

	// Read from a FIT file
	
	fpixel[3] = 1+start[T];
	fpixel[2] = 1;
	fpixel[1] = 1;
	fpixel[0] = 1;
	
	lpixel[3] = 1+start[T];
	lpixel[2] = isphere->dim[X];
	lpixel[1] = isphere->dim[Y];
	lpixel[0] = isphere->dim[Z];	
	
	fprintf(stdout, "   Reading FITS data for slice %d\n", i);
	fprintf(stdout, "      FPIXEL: %ld %ld %ld %ld\n", 
		fpixel[0], fpixel[1], fpixel[2], fpixel[3]);
	fprintf(stdout, "      LPIXEL: %ld %ld %ld %ld\n", 
		lpixel[0], lpixel[1], lpixel[2], lpixel[3]);
	
	fits_read_subset(isphere->fpfit, TFLOAT, fpixel, lpixel, inc,
			 &nulval, isphere->tslice[i], &anynul, &status);
	if (status!=0) {
	  fprintf(stdout, "error reading subset from fits file\n");

	}
      }
  }

  // now interpolate to specific time
  // first get time as p in range [0..1]
  p = (t - isphere->origin[T]) / isphere->delta[T] - index;
  
  // 4 point Lagrangian interpolation
  // see Abramowitz & Stegun eq. 25.2.13
  c[0] = -p * (p - 1) * (p - 2) / 6;
  c[1] =  (p * p - 1) * (p - 2) / 2;
  c[2] = -p * (p + 1) * (p - 2) / 2;
  c[3] =        p * (p * p - 1) / 6;
  
  // collapse 4D space into 3D via interpolation on time
  for (ix=0; ix<isphere->dim[X]; ix++)
    for (iy=0; iy<isphere->dim[Y]; iy++)
      for (iz=0; iz<isphere->dim[Z]; iz++) {
	i = (ix * isphere->dim[Y] + iy) * isphere->dim[Z] + iz;
	isphere->slice[i] = 0;
	for (it=0; it<4; it++) {
	  isphere->slice[i] +=  *(isphere->tslice[it] + i) * c[it];
	  if (isnan(*(isphere->tslice[it] + i))) fprintf(stdout,"nan in interp_t: i: %d\n",i);
	}
      }
  
  return (0);
}


/**********************************************************************
 *  ion_density - interpolates a function (of 3 independent variables)*
 *                that describes the ionosphereic electron density    *
 *                                                                    *
 * Algorithm:                                                         *
 *   Lagrange 4-point interpolation (Abramowitz & Stegun 25.2.13)     *
 *   simultaneously done in all 3 dimensions                          *
 *   with interpolation point lying in the central segment cube,      *
 *   if possible. If the data array has only one point "exterior" to  *
 *   the desired extrapolation point, then the interpolation will be  *
 *   done anyway, but to slightly lower accuracy.                     *
 * Inputs:                                                            *
 *   isphere - points to an ion_model struct containing the           *
 *             model to be interpolated                               *
 *   x       - 3-vector giving coords. at which interpolated output   *
 *             is desired. Units: (km, km, km)                        *
 * Outputs:                                                           *
 *   value   - pointer to interpolated value of the function at       *
 *             location x[0..2]. (elec/m^3)                           *
 *   ion_density (returned value of function)                         *
 *      0: no error                                                   *
 *      1: out of bounds error (x outside of grid)                    *
 * First created                                rjc  2002.12.3        *
 *                                                                    *
 *********************************************************************/
 int ion_density (struct ion_model *isphere,
                 const double x[3],
                 double *value)
{
    int index,
        i,
        ix, ixoff,
        iy, iyoff,
        iz, izoff;

    double p,
           cx[4],
           cy[4],
           cz[4];


    // Given the ENH location of the point relative to the center of the 
    // ionospheric model calculate the lat/lon/alt for
    // interpolation in the global ionosphere

    

                                    // x coordinate (East longitude)
                                    // calculate indices to be used
    index = (x[0] - isphere->origin[Y]) / isphere->delta[Y];
                                    // put interpolation point in range [0..1]
    p = (x[0] - isphere->origin[Y]) / isphere->delta[Y] - index;

    if (index < 0 || index >= isphere->dim[Y]) {
        printf ("x = %g is not within ionosphere model domain (%g, %g)\n",
                x[0], isphere->origin[Y], 
                isphere->origin[Y] + (isphere->dim[Y] + 1) * isphere->delta[Y] );
        return (1);
    }
    else if (index == 0)            // in left-most segment?
        {
        ixoff = index;              // ixoff is index of 1st of 4 pts to be used
        p -= 1.0;                   // p is between [-1..0]
        }
    else if (index == isphere->dim[Y]-1)
        {
        ixoff = index - 2;          // in right-most segment
        p += 1.0;                   // p is in range [1..2]
        }
    else
        ixoff = index - 1;          // in center segment

    
                                    // 4 point Lagrangian interpolation
                                    // see Abramowitz & Stegun eq. 25.2.13
    cx[0] = -p * (p - 1) * (p - 2) / 6;
    cx[1] =  (p * p - 1) * (p - 2) / 2;
    cx[2] = -p * (p + 1) * (p - 2) / 2;
    cx[3] =        p * (p * p - 1) / 6;

                                    // y coordinate (latitude)
                                    // calculate indices to be used
    index = (x[1] - isphere->origin[X]) / isphere->delta[X];
                                    // put interpolation point in range [0..1]
    p = (x[1] - isphere->origin[X]) / isphere->delta[X] - index;

    if (index < 0 || index >= isphere->dim[X]-1) {
        printf ("y = %g is not within ionosphere model domain (%g, %g)\n",
                x[1], isphere->origin[X], 
                isphere->origin[X] + (isphere->dim[X] + 1) * isphere->delta[X] );
        return (1);
        }
    else if (index == 0)            // in left-most segment?
        {
        iyoff = index;              // iyoff is index of 1st of 4 pts to be used
        p -= 1.0;                   // p is between [-1..0]
        }
    else if (index == isphere->dim[X]-2)
        {
        iyoff = index - 2;          // in right-most segment
        p += 1.0;                   // p is in range [1..2]
        }
    else
        iyoff = index - 1;          // in center segment

    
                                    // 4 point Lagrangian interpolation
                                    // see Abramowitz & Stegun eq. 25.2.13
    cy[0] = -p * (p - 1) * (p - 2) / 6;
    cy[1] =  (p * p - 1) * (p - 2) / 2;
    cy[2] = -p * (p + 1) * (p - 2) / 2;
    cy[3] =        p * (p * p - 1) / 6;

                                    // z coordinate
                                    // calculate indices to be used
    // Note that the action taken for the z-coordinate depends on whether 
    // the sampling in z is linear (isphere->zlog=0) or log (isphere->zlog=1)
    if (isphere->zlog == 1) {

      // Z is log spaced
      index = (log10(x[2]) - log10(isphere->origin[Z])) / isphere->delta[Z];
      // put interpolation point in range [0..1]
      p = (log10(x[2]) - log10(isphere->origin[Z])) / isphere->delta[Z] - index;

    } else {

      // Z is linearly spaced
      index = (x[2] - isphere->origin[Z]) / isphere->delta[Z];
      // put interpolation point in range [0..1]
      p = (x[2] - isphere->origin[Z]) / isphere->delta[Z] - index;

    }


    if (index < 0 || index >= isphere->dim[Z]-1)
      {
	if (isphere->zlog == 1) {
	  printf ("z = %g is not within ionosphere model domain (%g, %g)\n",
		  x[2], isphere->origin[Z], 
		  pow(10.0,log10(isphere->origin[Z]) + (isphere->dim[Z] + 1) * isphere->delta[Z]));
	} else {
	  printf ("z = %g is not within ionosphere model domain (%g, %g)\n",
		  x[2], isphere->origin[Z], 
		  isphere->origin[Z] + (isphere->dim[Z] + 1) * isphere->delta[Z] );
	}
        return (1);
      }
    else if (index == 0)            // in left-most segment?
      {
        izoff = index;              // izoff is index of 1st of 4 pts to be used
        p -= 1.0;                   // p is between [-1..0]
      }
    else if (index == isphere->dim[Z]-2)
      {
        izoff = index - 2;          // in right-most segment
        p += 1.0;                   // p is in range [1..2]
      }
    else
      izoff = index - 1;          // in center segment
    
    
    // 4 point Lagrangian interpolation
                                    // see Abramowitz & Stegun eq. 25.2.13
    cz[0] = -p * (p - 1) * (p - 2) / 6;
    cz[1] =  (p * p - 1) * (p - 2) / 2;
    cz[2] = -p * (p + 1) * (p - 2) / 2;
    cz[3] =        p * (p * p - 1) / 6;


    *value = 0;
    for (ix=0; ix<4; ix++)
        for (iy=0; iy<4; iy++)
            for (iz=0; iz<4; iz++) {
                i = ((ix + ixoff) * isphere->dim[Y] + (iy + iyoff)) * isphere->dim[Z] + iz + izoff;
                *value += cx[ix] * cy[iy] * cz[iz] * isphere->slice[i];
		if (isnan(isphere->slice[i])) fprintf(stdout,"nan for i: x,y,z: %d,%d,%d\n",ix,iy,iz);
            }


    return (0);                     // signify no problems
}


/**********************************
set parameters and allocated memory for station-based ionosphere phase screen
**********************************/
int init_ionoscreen(const int num_stations, 
		    const int gridsize_x, const int gridsize_y, 
		    struct iono_screen *ionoscreen) {

  int i;
  
  /* alloc memory in ionoscreen */
  ionoscreen->num_stations = num_stations;
  ionoscreen->size_x       = gridsize_x;
  ionoscreen->size_y       = gridsize_y;
  ionoscreen->valid = calloc(num_stations,sizeof(int));
  ionoscreen->phase_delay = calloc(num_stations,sizeof(float *));
  if (ionoscreen->valid==NULL || ionoscreen->phase_delay==NULL) return -1;
  for(i=0; i<num_stations; i++) {
    ionoscreen->phase_delay[i] = calloc(gridsize_x*gridsize_y,sizeof(float));
  }
  return 0;
}


/***********************
 import the ionsphere model
 *************************/
/**************************************************************************
 * open_ionosphere opens the ionospheric model file, and saves the        *
 *                 parameters found there in a structure. It also         *
 *                 allocates space for 4 time slices of 3D space, and     *
 *                 reads data from the first 4 times into memory.         *
 *                                                                        *
 * Input and output:                                                      *
 *   isphere  pointer to structure containing information relating to     *
 *            the ionosphere model.                                       *
 *                                                                        *
 * Function value returns non-zero in case of error, zero otherwise.      *
 *                                                                        *
 *                                                                        *
 *  first created                                 rjc  2002.11.21         *
 *************************************************************************/
int open_ionosphere (struct ion_model *isphere)
{
  int time_id, x_id, y_id, z_id, ret, i;
  size_t start[4],count[4];

  int n_hour, n_lat, n_lon, n_height;
  double hr_start, hr_end, lat_cent, lat_wid, lon_cent, lon_wid, height_start, height_end;
  char *comment=NULL;
  FILE *fpset=NULL;
  char *substr=NULL;
  int status=0,anynul=0;
  double dtor = 0.0174532925;
  long inc[4] = { 1, 1, 1, 1 };
  long fpixel[4]={1,1,1,1};
  long lpixel[4]={1,1,1,1};
  char line[IONO_MAX_FNAME], tmpline[IONO_MAX_FNAME];
  char *tag, *pch; 
  
  struct tcube_str *tcube=NULL;

  double nulval=0.0;
  int nret;

  // Commands to force the ionosphere to a particule time in UT
  int do_force=0;     // Do or do not force
  double force_ut=0.0;  // UT in seconds

  int do_turbulence=0,tcube_elements=128;
  double tcube_vx=0.1, tcube_vy=0.1, tcube_vz=0.0,
    tcube_iscale=25.0, tcube_oscale=5.0, 
    tcube_resolution=1.0, tcube_rms=0.01;
 
  int do_tid=0;
  double tid_amp=0.1, tid_amin=100.0, tid_amax=300.0, 
    tid_vel=0.15, tid_len=200.0;

  int do_simple=0;
  double simple_grad=1000.0, simple_dc=1E13, 
    simple_zmin=100.0, simple_zmax=300.0;

  /* Read in options from file */
  msg("opening iono settings file <%s>",2,isphere->settings_fname);
  fpset = fopen(isphere->settings_fname, "r");
  if (fpset == NULL) {
    fprintf(stderr, "Error reading ionosphere settings from file <%s>\n", isphere->settings_fname);
    exit(1);
  }

  /* read info from file */
  while(!feof(fpset)) {
    fgets(line,IONO_MAX_FNAME,fpset);
    strcpy(tmpline, line); /* make a copy of the original line */

    /* ignore comments etc */
    if (line[0]=='\n' || line[0]==';' || line[1]==';') {
      // if(debug ==1) fprintf(stderr," Skipping comment: %s",line);
      continue;
    }
    
    /* Extract the string for a single density from the new line,
       convert it into a float, exponentiate it so it is in units of
       electrons per cm^3, and place it into the list */
    tag = strtok (tmpline," ");

    /* Get the next density from the string for this line */
    pch = strtok(NULL, " ");
    
    if (strcmp(tag, "DO_TURBULENCE") == 0) {
      nret = sscanf(pch,"%d",&do_turbulence); 
      fprintf(stdout, "DO_TURBULENCE: %d\n", do_turbulence);
    }

    if (strcmp(tag, "TCUBE_VX") == 0) {
      nret = sscanf(pch,"%lf",&tcube_vx); 
      fprintf(stdout, "TCUBE_VX: %lf\n", tcube_vx);
    }

    if (strcmp(tag, "TCUBE_VY") == 0) {
      nret = sscanf(pch,"%lf",&tcube_vy); 
      fprintf(stdout, "TCUBE_VY: %lf\n", tcube_vy);
    }

    if (strcmp(tag, "TCUBE_VZ") == 0) {
      nret = sscanf(pch,"%lf",&tcube_vz); 
      fprintf(stdout, "TCUBE_VZ: %lf\n", tcube_vz);
    }

    if (strcmp(tag, "TCUBE_ISCALE") == 0) {
      nret = sscanf(pch,"%lf",&tcube_iscale); 
      fprintf(stdout, "TCUBE_ISCALE: %lf\n", tcube_iscale);
    }

    if (strcmp(tag, "TCUBE_OSCALE") == 0) {
      nret = sscanf(pch,"%lf",&tcube_oscale); 
      fprintf(stdout, "TCUBE_OSCALE: %lf\n", tcube_oscale);
    }

    if (strcmp(tag, "TCUBE_ELEMENTS") == 0) {
      nret = sscanf(pch,"%d",&tcube_elements); 
      fprintf(stdout, "TCUBE_ELEMENTS: %d\n", tcube_elements);
    }

    if (strcmp(tag, "TCUBE_RESOLUTION") == 0) {
      nret = sscanf(pch,"%lf",&tcube_resolution); 
      fprintf(stdout, "TCUBE_RESOLUTION: %lf\n", tcube_resolution);
    }

    if (strcmp(tag, "TCUBE_RMS") == 0) {
      nret = sscanf(pch,"%lf",&tcube_rms); 
      fprintf(stdout, "TCUBE_RMS: %lf\n", tcube_rms);
    }

    if (strcmp(tag, "DO_TID") == 0) {
      nret = sscanf(pch,"%d",&do_tid); 
      fprintf(stdout, "DO_TID: %d\n", do_tid);
    }

    if (strcmp(tag, "TID_AMP") == 0) {
      nret = sscanf(pch,"%lf",&tid_amp); 
      fprintf(stdout, "TID_AMP: %lf\n", tid_amp);
    }

    if (strcmp(tag, "TID_AMIN") == 0) {
      nret = sscanf(pch,"%lf",&tid_amin); 
      fprintf(stdout, "TID_AMIN: %lf\n", tid_amin);
    }

    if (strcmp(tag, "TID_AMAX") == 0) {
      nret = sscanf(pch,"%lf",&tid_amax); 
      fprintf(stdout, "TID_AMAX: %lf\n", tid_amax);
    }

    if (strcmp(tag, "TID_VEL") == 0) {
      nret = sscanf(pch,"%lf",&tid_vel); 
      fprintf(stdout, "TID_VEL: %lf\n", tid_vel);
    }

    if (strcmp(tag, "TID_LEN") == 0) {
      nret = sscanf(pch,"%lf",&tid_len); 
      fprintf(stdout, "TID_LEN: %lf\n", tid_len);
    }

    if (strcmp(tag, "DO_SIMPLE") == 0) {
      nret = sscanf(pch,"%d",&do_simple); 
      fprintf(stdout, "DO_SIMPLE: %d\n", do_simple);
    }

    if (strcmp(tag, "SIMPLE_GRAD") == 0) {
      nret = sscanf(pch,"%lf",&simple_grad); 
      fprintf(stdout, "SIMPLE_GRAD: %lf\n", simple_grad);
    }

    if (strcmp(tag, "SIMPLE_DC") == 0) {
      nret = sscanf(pch,"%lf",&simple_dc); 
      fprintf(stdout, "SIMPLE_DC: %lf\n", simple_dc);
    }

    if (strcmp(tag, "SIMPLE_ZMIN") == 0) {
      nret = sscanf(pch,"%lf",&simple_zmin); 
      fprintf(stdout, "SIMPLE_ZMIN: %lf\n", simple_zmin);
    }

    if (strcmp(tag, "SIMPLE_ZMAX") == 0) {
      nret = sscanf(pch,"%lf",&simple_zmax); 
      fprintf(stdout, "SIMPLE_ZMAX: %lf\n", simple_zmax);
    }

    if (strcmp(tag, "DO_FORCE") == 0) {
      nret = sscanf(pch,"%d",&do_force); 
      fprintf(stdout, "DO_FORCE: %d\n", do_force);
    }

    if (strcmp(tag, "FORCE_UT") == 0) {
      nret = sscanf(pch,"%lf",&force_ut); 
      fprintf(stdout, "FORCE_UT: %lf\n", force_ut);
    }
  }
  fclose(fpset);

  // Forcing time in UT
  isphere->do_force = do_force;
  isphere->force_ut = force_ut;

  // Place simple ionosphere model into isphere
  isphere->do_simple = do_simple;
  isphere->simple_grad = simple_grad;
  isphere->simple_dc = simple_dc;
  isphere->simple_zmin = simple_zmin;
  isphere->simple_zmax = simple_zmax;

  // Place TID info into isphere
  isphere->do_tid = do_tid;
  isphere->tid_amp = tid_amp;
  isphere->tid_amin = tid_amin;
  isphere->tid_amax = tid_amax;
  isphere->tid_vel = tid_vel;
  isphere->tid_len = tid_len;
  
  if (do_turbulence == 1) {

    isphere->do_turbulence = 1;

    /* allocate a turbulence cube data structure */
    isphere->tcube = calloc(1,sizeof(struct tcube_str));
    if( isphere->tcube==NULL) {
        fprintf(stderr, "ERROR: no malloc for turbulence cube\n");
        return -1;
    }
    tcube = isphere->tcube;

    /* Create tcube based on input values */
    if (debug) {
      fprintf(stdout, "Initializing the turbulence cube\n");
    }
    tcube_init_tcube_str( tcube, tcube_resolution,
			  tcube_elements, tcube_elements, tcube_elements, 
			  tcube_iscale, tcube_oscale, tcube_rms);
    
    // fprintf(stderr, "Dimensions X: %d Y: %d Z: %d\n", tcube->nx, 
    // tcube->ny, tcube->nz);
    
    /* Fill the turbulence cube */
    if (debug) {
      fprintf(stdout, "Fill the turbulence cube\n");
    }
    tcube_fill_turbulence( tcube );
    
    /* Assign the tcube to the isphere structure */
    isphere->tcube = tcube;
  }


  /* Figure out if the input ionospheric file is a netCDF or a FITS file 
     based on the file extension */
  substr = strstr(isphere->iono_fname,"fit");
  if ( substr !=NULL ) {
    fprintf(stdout, "The global ionosphere is an IRI model in a FIT file\n");
    isphere->global_type=1;
    isphere->zlog=0;
  } else {
    fprintf(stdout, "This global ionosphere is in a NetCDF file\n");
    isphere->global_type=0;
    isphere->zlog=1;  // Assume log z heights but double check when loading in
  }    
  
  // Read in the model using rhe appropriate routines depending on 
  // whether the file is netCDF or FITS
  if(isphere->global_type == 0) {
	/* open the file. It must remain open for the duration of the program
       since different time slices may be read in the future */
    ret = nc_open (isphere->iono_fname, NC_NOWRITE, &isphere -> file_id);
    if (ret != NC_NOERR)  {
      msg("Error in nc_open, ret = %d",2, ret);
      return (-1);
    }
    msg("Ionosphere data file %s opened.",1, isphere->iono_fname);
    // Retrieve all global attributes
    ret = nc_get_att_double (isphere -> file_id, NC_GLOBAL, "start_time", &isphere-> origin[T])
      | nc_get_att_double (isphere -> file_id, NC_GLOBAL, "delta_time", &isphere-> delta[T])
      | nc_get_att_double (isphere -> file_id, NC_GLOBAL, "start_lat",  &isphere-> origin[X])
      | nc_get_att_double (isphere -> file_id, NC_GLOBAL, "delta_lat",  &isphere-> delta[X])
      | nc_get_att_double (isphere -> file_id, NC_GLOBAL, "start_lon",  &isphere-> origin[Y])
      | nc_get_att_double (isphere -> file_id, NC_GLOBAL, "delta_lon",  &isphere-> delta[Y])
      | nc_get_att_double (isphere -> file_id, NC_GLOBAL, "start_alt",  &isphere-> origin[Z])
      | nc_get_att_double (isphere -> file_id, NC_GLOBAL, "delta_alt",  &isphere-> delta[Z]);
    
    if (ret != NC_NOERR) {
      msg("Error %d getting global attributes",2, ret);
      return (-1);
    }
    
    
    // Get variable ID and dimension IDs
    ret = nc_inq_varid (isphere -> file_id, "ionosphere_grid", &isphere -> grid_id)
      | nc_inq_dimid (isphere -> file_id, "t", &time_id)
      | nc_inq_dimid (isphere -> file_id, "x", &x_id)
      | nc_inq_dimid (isphere -> file_id, "y", &y_id)
      | nc_inq_dimid (isphere -> file_id, "z", &z_id);
    
    if (ret != NC_NOERR)
      {
        msg("Error %d getting dimension IDs", 2,ret);
        return (-1);
      }
    // Get dimension lengths
    ret = nc_inq_dimlen (isphere -> file_id, time_id, &isphere -> dim[T])
      | nc_inq_dimlen (isphere -> file_id, x_id,    &isphere -> dim[X])
      | nc_inq_dimlen (isphere -> file_id, y_id,    &isphere -> dim[Y])
      | nc_inq_dimlen (isphere -> file_id, z_id,    &isphere -> dim[Z]);
    
    if (ret != NC_NOERR) {
      msg("Error %d getting dimension lengths", 2,ret);
        return (-1);
    }
    

    

    // Allocate space and read first 4 time slices
    for (i=X; i<=Z; i++) {
      start[i] = 0;
      count[i] = isphere->dim[i];
    }
    
    count[T] = 1;
    
    isphere->slice = calloc (isphere->dim[X] * isphere->dim[Y] * isphere->dim[Z], sizeof (float));
    if (isphere->slice == NULL) {
      msg("Error allocating memory for timeslice",2);
      return (-1);
    }
    
    for (i=0; i<4; i++) {
        isphere->tslice[i] = calloc (isphere->dim[X] * isphere->dim[Y] * isphere->dim[Z], sizeof (float));
        if (isphere->tslice[i] == NULL) {
            msg("Error allocating memory for timeslice %d",2, i);
            return (-1);
	}
	// nor read data for this slice
        start[T] = i;
        ret = nc_get_vara_float(isphere->file_id, isphere->grid_id, start, count, isphere->tslice[i]);
        if (ret != NC_NOERR) {
            msg("Error %d reading time slice %d",2, ret, i);
            return (-1);
	}
    }
  
  } else {
    /* input is a FIT file */

	/* open the file. It must remain open for the duration of the program
       since different time slices may be read in the future */
    fits_open_file(&isphere->fpfit, isphere->iono_fname, 0, &status);
 
    msg("Ionosphere data file %s opened.",1, isphere->iono_fname);
    // Retrieve all global attributes
    
    fits_read_key(isphere->fpfit, TINT, "NHOUR" , &n_hour, comment, &status);
    fits_read_key(isphere->fpfit, TDOUBLE, "HRSTART" , &hr_start, comment, &status);
    fits_read_key(isphere->fpfit, TDOUBLE, "HREND" , &hr_end, comment, &status);
    
    fits_read_key(isphere->fpfit, TINT, "NLAT" , &n_lat, comment, &status);
    fits_read_key(isphere->fpfit, TDOUBLE, "LATCEN" , &lat_cent, comment, &status);
    fits_read_key(isphere->fpfit, TDOUBLE, "LATWID" , &lat_wid, comment, &status);
    
    fits_read_key(isphere->fpfit, TINT, "NLON" , &n_lon, comment, &status);
    fits_read_key(isphere->fpfit, TDOUBLE, "LONCEN" , &lon_cent, comment, &status);
    fits_read_key(isphere->fpfit, TDOUBLE, "LONWID" , &lon_wid, comment, &status);
    
    fits_read_key(isphere->fpfit, TINT, "NHEIGHT" , &n_height, comment, &status);
    fits_read_key(isphere->fpfit, TDOUBLE, "HSTART" , &height_start, comment, &status);
    fits_read_key(isphere->fpfit, TDOUBLE, "HEND" , &height_end, comment, &status);
    fits_read_key(isphere->fpfit, TINT, "ZLOG" , &isphere->zlog, comment, &status);
    
    fprintf(stdout, "Hour Num/Start/End: %d %lf %lf\n", n_hour, hr_start, hr_end);
    fprintf(stdout, "Lat Num/Cen/Wind: %d %lf %lf\n", n_lat, lat_cent, lat_wid);
    fprintf(stdout, "Lon Num/Cen/Wind: %d %lf %lf\n", n_lon, lon_cent, lon_wid);
    fprintf(stdout, "Height Num/Start/End: %d %lf %lf\n", n_height, height_start, height_end);

    // Starting time in seconds of day, time step in seconds of day
    isphere->origin[T] = hr_start*3600.0;
    isphere->delta[T] = (hr_end-hr_start)*3600.0/(n_hour-1.0);

    // Starting latitude in radians, lat step in radians
    isphere->origin[X] = (lat_cent-lat_wid/2.0)*dtor;
    isphere->delta[X] = lat_wid*dtor/(n_lat-1.0);
    isphere->center_lat = lat_cent*dtor;

    // Starting longitude in radians, lon step in radians
    isphere->origin[Y] = (lon_cent-lon_wid/2.0)*dtor;
    isphere->delta[Y] = lon_wid*dtor/(n_lon-1.0);
    isphere->center_lon = lon_cent*dtor;

    // Starting height in km, height step in LOG10
    isphere->origin[Z] = height_start;
    if (isphere->zlog) {
      isphere->delta[Z] = (log10(height_end)-log10(height_start))/(n_height-1.0);
    } else {
      isphere->delta[Z] = (height_end-height_start)/(n_height-1.0);
    }

    // Set dimension lengths
    isphere->dim[T] = n_hour;
    isphere->dim[X] = n_lat;
    isphere->dim[Y] = n_lon;
    isphere->dim[Z] = n_height;
    

    // Allocate space and read first 4 time slices
    for (i=X; i<=Z; i++) {
      start[i] = 0;
      count[i] = isphere->dim[i];
    }
    
    count[T] = 1;
    
    fprintf(stdout, "Allocating space for timeslice\n");
    isphere->slice = calloc (isphere->dim[X] * isphere->dim[Y] * isphere->dim[Z], sizeof (float));
    if (isphere->slice == NULL) {
      msg("Error allocating memory for timeslice",2);
      return (-1);
    }
    
    fprintf(stdout, "Loading in timeslices\n");
    for (i=0; i<4; i++)                 // i is a time-like index here
      {

	fprintf(stdout, "   Creating space for slice %d with %d elements\n", i, 
		(int)(isphere->dim[X] * isphere->dim[Y] * isphere->dim[Z]));

	//        isphere->tslice[i] = (double *) 
	//  calloc (isphere->dim[X] * isphere->dim[Y] * isphere->dim[Z], sizeof (double));

	isphere->tslice[i] = calloc(isphere->dim[X] * isphere->dim[Y] * isphere->dim[Z], sizeof (double));


        if (isphere->tslice[i] == NULL)
	  {
            msg("Error allocating memory for timeslice %d",2, i);
            return (-1);
	  }
	// nor read data for this slice
        start[T] = i;

	fpixel[3] = 1+i;
	fpixel[2] = 1;
	fpixel[1] = 1;
	fpixel[0] = 1;
	
	lpixel[3] = 1+i;
	lpixel[2] = n_lat;
	lpixel[1] = n_lon;
	lpixel[0] = n_height;	
	
	//fprintf(stderr, "   Reading FITS data for slice %d\n", i);
	//fprintf(stderr, "      FPIXEL: %ld %ld %ld %ld\n", 
	//	fpixel[0], fpixel[1], fpixel[2], fpixel[3]);
	//fprintf(stderr, "      LPIXEL: %ld %ld %ld %ld\n", 
	//	lpixel[0], lpixel[1], lpixel[2], lpixel[3]);

	fits_read_subset(isphere->fpfit, TFLOAT, fpixel, lpixel, inc,
			 &nulval, isphere->tslice[i], &anynul, &status);
	if (status!=0) {
	  fprintf(stderr, "error reading subset from fits file\n");
	}
	
	fprintf(stdout, "   FITS read complete %d\n", i);

    }
  }


  // Center of global ionosphere model
  isphere->center_lat = isphere->origin[X]+isphere->delta[X]*0.5*(isphere->dim[X]-1.0);
  isphere->center_lon = isphere->origin[Y]+isphere->delta[Y]*0.5*(isphere->dim[Y]-1.0);
    

  
  // Report summary of global ionosphere
  msg("Ionosphere start_time: %g,\tdelta_time: %g\t dim_size: %d",1,
      isphere-> origin[T],isphere-> delta[T],isphere -> dim[T]);
  msg("Ionosphere start_lat:  %g,\tdelta_lat:  %g\t dim_size: %d",1,
      isphere-> origin[X]/dtor,isphere-> delta[X]/dtor,isphere -> dim[X]);
  msg("Ionosphere start_lon:  %g,\tdelta_lon:  %g\t dim_size: %d",1,
      isphere-> origin[Y]/dtor,isphere-> delta[Y]/dtor,isphere -> dim[Y]);
  msg("Ionosphere start_alt:  %g,\tdelta_alt:  %g\t dim_size: %d",1,
      isphere-> origin[Z],isphere-> delta[Z],isphere -> dim[Z]);
  msg("Ionosphere center latitude: %g  longitude:  %g",1,
      isphere->center_lat/dtor, isphere->center_lon/dtor);
  
  if (isphere->zlog) {
    for (i=T; i<=Y; i++)
      msg("domain[%d] %.3f %.3f",0, i, isphere->origin[i], 
	  isphere->origin[i] + (isphere->dim[i] - 1.0) * isphere-> delta[i]);
    msg("domain[%d] %.3f %.3f",0, Z, isphere->origin[Z], 
	pow(10.0, log10(isphere->origin[Z]) + 
	    (isphere->dim[Z] - 1.0) * isphere-> delta[Z]));
  } else {
    for (i=T; i<=Z; i++)
      msg("domain[%d] %.3f %.3f",0, i, isphere->origin[i], 
	  isphere->origin[i] + (isphere->dim[i] - 1.0) * isphere-> delta[i]);
  }

  
  isphere -> ntslice = 0;             // First slice in RAM has index of 0

  // Done!
  return (0);

}


  


int tcube_init_path_str(struct path_str *path, double x0, double y0, double z0, 
		  double nx, double ny, double nz, 
		  double step, int steps) {

  int status = 0;

  path->x0 = x0;
  path->y0 = y0;
  path->z0 = z0;

  path->nx = nx;
  path->ny = ny;
  path->nz = nz;

  path->steps = steps;
  path->step = step;

  return status;

}

void tcube_init_tcube_str(struct tcube_str *tcube, double s, 
			  int nx, int ny, int nz,
			  double i_scale, double o_scale, double rms) {

  /* A turbulent cube structure */
  /* tcube_str thiscube; */

  /* Number of elements in each dimension */
  /* int nx=128,ny=128,nz=128; */
  
  /* Step size in km */
  /* double s=1.0; */
  
  /* The big matrix */
  int bigmatrixsizes[4]={1,1,1,0}; // 0 terminated array of dimension size 
  
  fprintf(stdout, "Dimensions X: %d Y: %d Z: %d\n", nx, 
	  ny, nz);

  fprintf(stdout, "Step size X: %lf Y: %lf Z: %lf\n", s, s, s);

  /* Populate this cube */
  tcube->nx = nx;
  tcube->ny = ny;
  tcube->nz = nz;

  tcube->sx = s;
  tcube->sy = s;
  tcube->sz = s;

  tcube->i_scale = i_scale;
  tcube->o_scale = o_scale;
  tcube->rms = rms;

  fprintf(stdout, "Filled params\n");

  fprintf(stdout, "Dimensions X: %d Y: %d Z: %d\n", tcube->nx, 
	  tcube->ny, tcube->nz);

  /* Create cube */
  fprintf(stdout, "Creating matrix\n");
  
  bigmatrixsizes[0] = nx;
  bigmatrixsizes[1] = ny;
  bigmatrixsizes[2] = nz;
  
  tcube->cube = matloc(bigmatrixsizes,sizeof(float));

  fprintf(stdout, "Cube is initialized\n");

}

/* Fill this cube with a simple function  */
void tcube_fill_simple(struct tcube_str *tcube ) {

  int i, j, k;
  double new_value;

  for (i=0;i<tcube->nx;i++) {
    for (j=0;j<tcube->ny;j++) {
      for (k=0;k<tcube->nz;k++) {
	
	new_value = k;
	fprintf(stdout, " %d %d %d: %f\n", i, j, k, new_value);
	
	tcube->cube[i][j][k] = new_value;

      }
    }
  }

}


/* Fill this cube with 3D isotropic turbulence  */
void tcube_fill_turbulence(struct tcube_str *tcube) {

  int x_dim, y_dim, z_dim;
  double c, i_scale, o_scale;
  double dx, dy, dz;
  double rms_desired=0.0, rms_calc=0.0;

  complex *temp_iono;
  int i, j, k, index, tIndex;
  int dims[3];
  double dist;
  double amplitude, phase;

  /* fprintf(stderr, "Filling cube with turbulence model\n"); */

  x_dim = tcube->nx;
  y_dim = tcube->ny;
  z_dim = tcube->nz;
  dx = tcube->sx;
  dy = tcube->sy;
  dz = tcube->sz;

  o_scale = tcube->o_scale;
  i_scale = tcube->i_scale;

  c = 0.0045;

  temp_iono = calloc(x_dim*y_dim*z_dim, sizeof(complex)); 

  time_t t1;
  (void) time(&t1);
  srandom(-(long) t1);
  srandom(-1); /* keep random seed fixed for easier comparison */
  
  dims[0] = x_dim; dims[1] = y_dim; dims[2] = z_dim;
  for(i = 0; i < x_dim; i++) 
    for(j = 0; j < y_dim; j++)
      for(k = 0; k < z_dim; k++) {
        index = i*y_dim*z_dim + j*z_dim + k;
        tIndex = TranslateCoord2(i, x_dim)*y_dim*z_dim + 
	  TranslateCoord2(j, y_dim)*z_dim + TranslateCoord2(k, z_dim);
	
        if(index >= x_dim*y_dim*z_dim/2 + y_dim*z_dim/2 + z_dim/2) break;
	
        dist = pow(pow(x_dim/2-i, 2) + 
		   dx*x_dim/(dy*y_dim)*pow(y_dim/2-j, 2) + 
		   dx*x_dim/(dz*z_dim)*pow(z_dim/2-k, 2), .5);  
	
        amplitude = envelope_function(dist, o_scale, i_scale) * 
	  c * (pow(dist, -11.0/3.0));
        phase = 2.0*M_PI*((double)random()/(double)RAND_MAX); 
	temp_iono[index] = rect(amplitude, phase);
	temp_iono[tIndex] = c_conj(temp_iono[index]);
      }
  translate_array(x_dim, y_dim, z_dim, &temp_iono);
  fourc(temp_iono,dims, 3, -1);
  translate_array(x_dim, y_dim, z_dim, &temp_iono);
  
  for(i = 0; i < x_dim; i++) 
    for(j = 0; j < y_dim; j++)
      for(k = 0; k < z_dim; k++) {
 
	index = x_dim*y_dim*z_dim + i*y_dim*z_dim + j*z_dim + k;
	tcube->cube[i][j][k] = c_real(temp_iono[i*y_dim*z_dim + j*z_dim + k]);
	
      }
  
  
  free(temp_iono);

  /* Now set RMS to requested value */
  rms_desired = tcube->rms;
  rms_calc = tcube_calc_rms(tcube);
  fprintf(stdout, "RMS calculated: %f desired: %f\n", 
	  rms_calc, rms_desired);
  for(i = 0; i < x_dim; i++) 
    for(j = 0; j < y_dim; j++)
      for(k = 0; k < z_dim; k++) {
	tcube->cube[i][j][k] *= (rms_desired/rms_calc);
      }
  
  /* fprintf(stderr, "     Turbulence model complete\n"); */

}

/* Determine values of the density along a path */
int tcube_path( struct tcube_str *tcube, 
		struct path_str *path, 
		double *density) {
  
  int status=0;
  int i=0;
  int steps=0;
  double step=0.0;
  double x=0.0, y=0.0, z=0.0;
  double nx=0.0, ny=0.0, nz=0.0;

  /* Number of steps and step size */
  steps = path->steps;
  step = path->step;

  /* Set starting position */
  x = path->x0;
  y = path->y0;
  z = path->z0;

  /* Direction */
  nx = path->nx;
  ny = path->ny;
  nz = path->nz;

  for (i=0;i<steps;i++) {

    x += step*nx;
    y += step*ny;
    z += step*nz;

    density[i] = tcube_trilinear(tcube, x, y, z);
    fprintf(stdout, "%f %f %f %f\n", x, y, z, density[i]);

  }
  
  status = 1;

  return status;
}


/* Calculate the average value of the turbulence */
double tcube_calc_average( struct tcube_str *tcube ) {

  int i, j, k;
  int nx, ny, nz;
  double average=0.0;

  nx = tcube->nx;
  ny = tcube->ny;
  nz = tcube->nz;

  for (i=0;i<nx;i++) {
    for (j=0;j<ny;j++) {
      for (k=0;k<nz;k++) {
	average = average + tcube->cube[i][j][k];
      }
    }
  }

  return average/(nx*ny*nz);
}

/* Calculate the RMS value of the turbulence */
double tcube_calc_rms( struct tcube_str *tcube ) {

  int i, j, k;
  int nx, ny, nz;
  double average=0.0, rms=0.0;

  /* First get the average */
  average = tcube_calc_average(tcube);

  nx = tcube->nx;
  ny = tcube->ny;
  nz = tcube->nz;

  for (i=0;i<nx;i++) {
    for (j=0;j<ny;j++) {
      for (k=0;k<nz;k++) {
	rms = rms + pow(tcube->cube[i][j][k]-average,2);
      }
    }
  }

  return sqrt(rms/(nx*ny*nz));
}

/* Print this cube */
void tcube_print( struct tcube_str *tcube ) {

  int i, j, k;
  int print_contents=0;

  fprintf(stdout, "tcube summary\n");
  fprintf(stdout, " Dimensions X: %d Y: %d Z: %d\n", tcube->nx, tcube->ny, tcube->nz);
  fprintf(stdout, " Step sizes X: %f Y: %f Z: %f\n", tcube->sx, tcube->sy, tcube->sz);
  fprintf(stdout, " inner scale: %f   outer scale %f\n", tcube->i_scale, tcube->o_scale);


  if (print_contents == 1) {
    for (i=0;i<tcube->nx;i++) {
      
      fprintf(stdout, "X: %d\n", i);
      
      for (j=0;j<tcube->ny;j++) {
	
	fprintf(stdout, " Y: %d    ", j);
	
	for (k=0;k<tcube->nz;k++) {
	  
	  fprintf(stdout, "%f ", tcube->cube[i][j][k]);
	  
	}
	
	fprintf(stdout, "\n");
      }
    }
  }

}

void matfree(void *matt, int dim[])
{
  int bob;
  if (*dim == 0) return;
  if (dim[1]==0) free(matt);
  else {
    for(bob=0;bob<*dim;bob++)
      matfree(((void**)matt)[bob],dim+1);
    free(matt);
  }
}


/* recursively allocate a 3D array */
void *matloc(int dim[],size_t cellsize) {
  void **matt; int bob;

  /* nothing to do */
  if (dim[0]==0) return NULL;

  /* return a 1-D array if the next dimesion is zero */
  if (dim[1]==0) return (void*)malloc(dim[0]*cellsize);

  /* allocate pointers for the current dimension */
  matt = (void**)malloc(dim[0]*sizeof(void*));
  if (matt==NULL) return matt;

  /* for each element, go and allocate an array for the remainnig dimension */
  for(bob=0; bob<dim[0]; bob++)
    if ((matt[bob]=matloc(dim+1,cellsize)) == NULL) {
	/* something bad happened. Free and exit */
	while (bob--) matfree(matt[bob],dim+1);
	free(matt);
	return 0;
    }
  return matt;
} 


/* Perform trilinear interpolation using tcube */
double tcube_trilinear(struct tcube_str *tcube, 
		      double x, double y, double z) {

  double dx, dy, dz;
  double sx, sy, sz;
  int nx, ny, nz;

  
  /* Value of the interpolated function at each vertex */
  double V000, V100, V010, V001, V101, V011, V110, V111;

  /* Interpolated result */
  double result;

  /* Indices for the corners of the vertex */
  int i0, i1, j0, j1, k0, k1;
  
  /* extract parameters from tcube */
  sx = tcube->sx;
  sy = tcube->sy;
  sz = tcube->sz;

  nx = tcube->nx;
  ny = tcube->ny;
  nz = tcube->nz;

  i0 = floor(x/sx);
  dx = x - sx*i0;

  i0 = i0 % nx;
  if (i0 < 0) {
    i0 = (i0 + 2*nx) % nx;
  }
  i1 = (i0 + 1) % nx;

  j0 = floor(y/sy);
  dy = y - sy*j0;
  j0 = j0 % ny;
  if (j0 < 0) {
    j0 = (j0 + 2*nx) % nx;
  }
  j1 = (j0 + 1) % ny;

  k0 = floor(z/sz);
  dz = z - sz*k0;
  k0 = k0 % nz;
  if (k0 < 0) {
    k0 = (k0 + 2*nx) % nx;
  }
  k1 = (k0 + 1) % nz;


  V000 = tcube->cube[i0][j0][k0];
  V100 = tcube->cube[i1][j0][k0];
  V010 = tcube->cube[i0][j1][k0];
  V001 = tcube->cube[i0][j0][k1];
  V101 = tcube->cube[i1][j0][k1];
  V011 = tcube->cube[i0][j1][k1];
  V110 = tcube->cube[i1][j1][k0];
  V111 = tcube->cube[i1][j1][k1];

  result =  V000 * (1.0 - dx) * (1.0 - dy) * (1.0 - dz) +
    V100 * dx * (1.0 - dy) * (1.0 - dz) + 
    V010 * (1.0 - dx) * dy * (1.0 - dz) + 
    V001 * (1.0 - dx) * (1.0 - dy) * dz +
    V101 * dx * (1.0 - dy) * dz + 
    V011 * (1.0 - dx) * dy * dz + 
    V110 * dx * dy * (1.0 - dz) + 
    V111 * dx * dy * dz;
  
  /* 
     fprintf(stderr, "Location: (%f, %f, %f)\n", x, y, z);
     fprintf(stderr, " i0: %d i1: %d dx: %f\n", i0, i1, dx);
     fprintf(stderr, " j0: %d j1: %d dy: %f\n", j0, j1, dy);
     fprintf(stderr, " k0: %d k1: %d dz: %f\n", k0, k1, dz);
  */

  return result;

}



/* write a turbulent cube out as a fits file */
void tcube_write_fits(struct tcube_str *tcube, char *outfilename) {

  int status=0,verbose=0;
  
  fitsfile *fpout =NULL;
  long naxis=3;
  long naxes[3]={0,0,0},fpixel[3]={1,1,1},lpixel[3]={1,1,1};
  
  /* char *outfilename="tcube.fit"; */
  
  int i,j,k;
  int nx, ny, nz;

  double sx, sy, sz, i_scale, o_scale, rms;

  double *density=NULL;
  
  /* create the output FITS file */
  remove(outfilename); /* just overwrite any existing file */
  fits_create_file(&fpout, outfilename, &status);
  if (status!=0) {
    fprintf(stderr, "failed to create fits file %s\n",outfilename);
    exit(1);
  }

  /* Extract parameters from the tcube structure */

  nx = tcube->nx;
  ny = tcube->ny;
  nz = tcube->nz;

  sx = tcube->sx;
  sy = tcube->sy;
  sz = tcube->sz;
  
  i_scale = tcube->i_scale;
  o_scale = tcube->o_scale;
  rms = tcube->rms;

  /* Set size of image */
  naxes[2] = nx;
  naxes[1] = ny;
  naxes[0] = nz;

  
  /* Create the fits file */
  if (verbose == 1) {
    fprintf(stderr,"Creating FITS image file\n");
  }
  
  fits_create_img( fpout, DOUBLE_IMG, naxis, naxes, &status);
  if (status!=0) {
    fprintf(stderr, "failed to initialize DOUBLE_IMG in file %s\n",
	    outfilename);
    exit(1);
  }
  
  /* create one dimensional array to hold density as a function of height */
  density = calloc(nz,sizeof(double));
    
  for(i = 0; i < nx; i++) {
    for(j = 0; j < ny; j++) {
      
      for(k = 0; k < nz; k++) {
	density[k] = tcube->cube[i][j][k];

      }

      /* Write the line of heights to the FITS file */
      fpixel[2] = 1+i;
      fpixel[1] = 1+j;
      fpixel[0] = 1;
      
      lpixel[2] = 1+i;
      lpixel[1] = 1+j;
      lpixel[0] = nz;	

      fits_write_subset(fpout, TDOUBLE, fpixel, lpixel, density, &status);
      if (status!=0) {
	  fprintf(stderr,"unable to write image\n");    
	  exit(1);
      }

    }
  }

 /* Write metadata to fits file */
  if (verbose == 1) {
    fprintf(stderr, "Writing metadata to FITS file\n");
  }
  
  fits_update_key(fpout, TINT, "NX" , &nx, 
		  "[#] Number of samples in X", &status);
  fits_update_key(fpout, TDOUBLE, "DX" , &sx, 
		  "[km] Step size in X", &status);

  fits_update_key(fpout, TINT, "NY" , &ny, 
		  "[#] Number of samples in Y", &status);
  fits_update_key(fpout, TDOUBLE, "DY" , &sy, 
		  "[km] Step size in Y", &status);

  fits_update_key(fpout, TINT, "NZ" , &nz, 
		  "[#] Number of samples in Z", &status);
  fits_update_key(fpout, TDOUBLE, "DZ" , &sz, 
		  "[km] Step size in Z", &status);

  fits_update_key(fpout, TDOUBLE, "ISCALE" , &sz, 
		  "[TBD] Inner scale of turbulence", &status);
  fits_update_key(fpout, TDOUBLE, "OSCALE" , &sz, 
		  "[TBD] Outer scale of turbulence", &status);
  fits_update_key(fpout, TDOUBLE, "RMS" , &sz, 
		  "[%] Root mean square amplitude", &status);



  if (status!=0) {
    fprintf(stderr,"error: unable to update keyword\n");
    exit(1);
 }
  
  /* clean up */
  if (verbose == 1) {
    fprintf(stderr, "Closing FITS file\n");
  }
  fits_close_file(fpout, &status);
  if (status!=0) {
    fprintf(stderr,"error: status %d on closefile", status);
    exit(1);
  }
  
  
  /* Done! */
  if(density!=NULL) free(density);
}


/* The following few functions have been stolen and modified from the code in
   "generate_iono_turbulent.c" in the maps_ionosphere_generation folder. It is
   not intended to keep the code in sync with the other */
int TranslateCoord(int coord, int dim_length)
{
  coord -= floor(dim_length/2);
  if(coord < 0) coord += dim_length;
  return coord;
}

int TranslateCoord2(int coord, int dim_length)
{
  coord = dim_length - coord;
  if(coord == dim_length) coord = 0;
  return coord;
}


/*
 * The envelope function is multiplied with the general function for
 * the Kolmogorov Spectra in k-space to reduce the effect of waves
 * that are both too small and too large where these sizes are specified
 * by the o_scale (large wave limit) and i_scale (small wave limit)
 */
double envelope_function(double k, double o_scale, double i_scale)
{
  if(k < o_scale) return exp(k-o_scale);
  if(k > i_scale) return 0; 
  return 1; 
} 


/*
 * Translates the values in the passed complex array so that the origin, which is considered
 * to lie at the element with coordinates (x_dim/2, y_dim/2, z_dim/2) lies at the origin.
 * This is done so that the array can be passed to the fourc function.
 */
void translate_array(int x_dim, int y_dim, int z_dim, complex **temp_iono)
{
  complex *temp, *new_iono = malloc(sizeof(complex)*x_dim*y_dim*z_dim);
  int i, j, k;
  int x, y, z;
  int index, nIndex;

  for(i = 0; i < x_dim; i++) 
    for(j = 0; j < y_dim; j++) 
      for(k = 0; k < z_dim; k++) {
        index = i*y_dim*z_dim + j*z_dim + k;
        nIndex = index;
        x = floor(((double) nIndex)/(y_dim*z_dim));  nIndex -= x*y_dim*z_dim;
        y = nIndex/z_dim; nIndex -= y*z_dim;
        z = nIndex;
        x = TranslateCoord(x, x_dim); y = TranslateCoord(y, y_dim); z = TranslateCoord(z, z_dim);
        nIndex = x*y_dim*z_dim + y*z_dim + z;
        new_iono[nIndex] = (*temp_iono)[index];
      }
  temp = *temp_iono;
  *temp_iono = new_iono;
  free(temp);
}

