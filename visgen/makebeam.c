/******************************************************************************/
/*                                                                            */
/* This program calculates the complex beam formed by phasing up all elements */
/* in a station.  This routine needs to know the FOV center, the observing    */
/* time, the station position (long,lat), relative antenna positions within   */
/* the station, antenna weights, output array size, frequency of interest.    */
/* Output is a filled in structure consisting of beam values in two arrays,   */
/* one real, one complex.                                                     */
/*                                                                            */
/* Output : filled in stnbeam structure                                       */
/* Created March 2002 SSD                                                     */
/* Modified considerably : Oct 16 2002					      */
/* major modifications 2007 by Randall Wayth:                                 */
/*   - polarised response                                                     */
/*   - generic antenans                                                       */
/*   - separate functions to calculate the station beam vs fill the beamgrid  */
/******************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "fitsio.h"
#include "type_comp.h"
#include "comp_func.h"
#include "station_beam.h"
#include "utility.h"
#include "array_specification.h"
#include "observing_param.h"
#include "makebeam.h"


/* private function prototypes */
int primary_beam_response(struct antenna *ant,int receptor, int skypol_index, 
			  skypos *ref_dir,skypos *look_dir, float lambda,
			  complex *response );
double AngularDist(skypos *ref_dir,skypos *look_dir);
float bessj1(float x);
float end_fire_2element(float d, float th);
float jinc(float x);
int DumpStationbeamToFITS(struct   stnbeam *beam, 
			  struct stationspec *stn, int npols,double freq,
			  double gha,double dec_ref,double ra_ref,
			  float  crval_ra, float crval_dec,
			  float  cdelt_ra, float cdelt_dec);

/* private global vars */
static int debug_makebeam=0;

/* set module level debugging */
void makebeam_set_debug(int debug_level) {
    debug_makebeam = debug_level;
}

/******************************
 * Calculate the station beam Jones matrices for station stn for all the 
 * points on the sky specified in bgrid.
 * The result goes in beam.
 * stn (in): pointer to stationspec
 * beam (out): resulting station voltage response for all pol products 
 *             goes here
 * bgrid (in): this is a grid of points on the sky to sample the beam. Should 
 *             be sufficiently close to Nyquist sample the station beam, but 
 *             there is currently no check for this.
 * freq (in): observing freq in MHz
 * gha (in):  the Greenwhich HA (radian)
 * ra_ref,dec_ref (in): the station beam pointing direction (radian)
 * ra_ref_correlator:   the ra of the correlator phase reference poi
 * NOTE: currently bgrid->size is fixed at NUM_SAMPS=32 defined in 
 *       fill_beamgrid.c 
 ******************************/
int makebeam(/* stn, beam, bgrid, freq, gha, dec_ref, ra_ref) */
	     struct stationspec *stn,
	     struct stnbeam *beam,
	     struct beamgrid *bgrid,
	     double freq,
	     double gha,
	     struct observing_param *obsparam)
{ 
  int i,j,receptor,skypol, index, polindex=1, size, result=0,n_pols=0;
  double hang_ref, dec, ra;
  complex beam_res[4];
  /* unit vec and alt/az and ha/dec/lat of sky positions  *
   * for pointing ref direction.                          */
  skypos ref_dir;    


  size = bgrid->size;
  beam->size = bgrid->size; 
  beam->xref_pixel = bgrid->xref_pixel;
  beam->yref_pixel = bgrid->yref_pixel;
  n_pols = stn->layout->ants[0].num_receptors;

  /* make space for the beam grid(s). There will be either 1 or 4 
   * depending on whether we are using dual pol
   * or single pol receptors */
  if (n_pols < 1 || n_pols > 2) {
    msg("ERROR: makebeam: station has illegal number of pols. must be 1 or 2",
	3);
    return 1;
  }
  for (i=0; i< n_pols*n_pols; i++) {
    /* these get re-used, so only allocate once */
    if (beam->beam[i] ==NULL) beam->beam[i] = calloc(size*size,sizeof(complex));
    if (beam->beam[i] ==NULL){
      msg("ERROR: makebeam: no malloc",3);
      return -1;
    }
  }

				/* Figure out some useful quantities: 
				 * hour angle of FOV center, size of 
				 * beam sampling grid.		*/

  /* set up the reference dir for the pointing center of the station. 
   * Note that may be different from the phase center of the array 
   * (correlator). Each station will have a different location and since 
   * the stations follow the curvature of the Earth, have a slightly 
   * different HA.
   * The GHA is for the correlator phase center. We need to know the HA 
   * of the pointing center.
   * The HA(point) = LST - RA(point). We know the RA of the pointing 
   * center and correlator phase center.
   * LST = GHA(correlator) + lon + RA(correlator), so we use this to 
   * calculate the HA of the pointing center below */
  hang_ref = gha + stn->longitude + obsparam->phase_cent_RA - 
    obsparam->point_cent_RA;
  MakeSkyPos(hang_ref, obsparam->point_cent_DEC, stn->latitude, &ref_dir);

  /* loop over beam sampling grid, fill in beam grid */
  msg("Makebeam for station %d. There are %d antennas in the station. " \
      "HA ref: %g. GHA: %g, long: %g. Ref (xyz): (%g,%g,%g), alt/az: %g,%g",
      -1, beam->station, stn->layout->nant, hang_ref, gha,stn->longitude,
      ref_dir.x, ref_dir.y, ref_dir.z, ref_dir.alt, ref_dir.az);

  for (i=0; i<size; i++){
      for (j=0; j<size; j++)  {

        /* Get index of this beam pixel in 1-D complex array 	*/
        index = i*size + j;

        /* Read RA and DEC of this beam sample pixel	*/
        ra = bgrid->grid[index].ra;
        dec = bgrid->grid[index].dec;

        /* if this is not a valid sky coord, set result to zero */
        if (ra == BAD_POINT || dec ==BAD_POINT) {
          for (polindex=0; polindex < n_pols*n_pols; polindex++) {
            beam->beam[polindex][index] = c_zero();
          }
          continue;
        }

        /* Need hour angle at this ra, but don't have to */
        /* recompute gha to get this, just remember gha = gst - ra, */
        /* so take (gha+ra_ref-ra) to get gha at the offset point. */
        /* Or hang = hang_ref + (ra_ref - ra).  It's easy. */
        beam_res[0] = beam_res[1] = beam_res[2] = beam_res[3] = c_zero();
        //result = station_response (stn, &ref_dir, hang_ref + 
	//           (obsparam->phase_cent_RA - ra), dec, freq, beam_res);
        result = station_response(stn, &ref_dir, hang_ref + 
				  (obsparam->point_cent_RA - ra), dec, freq, 
				  beam_res);
        if (result != 0) return result;

        /* copy data into beamgrid. normalise by the array size 
	 * (for later FFT) */
        for (receptor=0; receptor< n_pols; receptor++) {
          for (skypol=0; skypol<n_pols; skypol++) {
            polindex = receptor*n_pols + skypol;
            beam->beam[polindex][index] = s_mult(beam_res[polindex],1.0/size);
          }
        }
      }
  }
  
  if(debug_makebeam && beam->station==0) {
    int bgrid_ind;
    float delta_ra=0,delta_dec=0;
    bgrid_ind = bgrid->xref_pixel + bgrid->size*bgrid->yref_pixel;
    if (bgrid->size>1) {
        /* reverse-engineer the pixel angular size. 
	 * Only works for square fields */
        delta_ra = delta_dec = asin(bgrid->grid[bgrid_ind+bgrid->size].dec - 
				    bgrid->grid[bgrid_ind].dec);
    }
    DumpStationbeamToFITS(beam, stn, n_pols*n_pols, freq, gha,
			  obsparam->phase_cent_DEC, obsparam->phase_cent_RA,
			  bgrid->grid[bgrid_ind].ra*(180.0/M_PI),
			  bgrid->grid[bgrid_ind].dec*(180.0/M_PI),
			  delta_ra*(180.0/M_PI),delta_dec*(180.0/M_PI));
  }

  beam->valid_flag = 1;
  return 0;
}


/******************************
 * Calculate the combined phased-array response from all the elements 
 * in a station beam.
 * stn: input pointer to the station
 * ref_dir: pointer to skypos with the phase reference direction for 
 *          the station (not of the correlator)
 * hang: the HA of the look direction
 * dec: the DEC of the look direction
 * freq: observation freq in MHz
 * beam: resulting complex Jones matrix response.
 ******************************/
int station_response(struct stationspec *stn, skypos *ref_dir, double hang, 
		     double dec, double freq, complex beam[4]) {

  skypos look_dir;   /* unit vec and alt/az and ha/dec/lat of sky positions */
  double trigarg,wav;
  double dx,dy,dz;
  double wtsum[2];
  complex *phases=NULL,ant;
  int k,receptor,skypol,n_pols,polindex,result;

  n_pols = stn->layout->ants[0].num_receptors;

  MakeSkyPos(hang,dec, stn->latitude, &look_dir);

  /* offset vector between reference position and look direction */
  dx = look_dir.x - ref_dir->x;
  dy = look_dir.y - ref_dir->y;
  dz = look_dir.z - ref_dir->z;
  if (debug_makebeam) {
    printf("Look dir: HA,DEC: %g, %g. ENU: %g,%g,%g. Ref dir: %g,%g,%g\n",
	   hang, dec, look_dir.x, look_dir.y, look_dir.z,
	   ref_dir->x,ref_dir->y,ref_dir->z);
    printf("HA,dec,lat: %.10g,%.10g,%.10g alt: %.10g, az: %.10g, " \
	   "(x,y,z): (%.10g,%.10g,%.10g)\n",
	   hang, dec, stn->latitude, look_dir.alt, look_dir.az, 
	   look_dir.x, look_dir.y, look_dir.z);
  }

  /* pre-calculate phases and weights */
  wav=freq*(2.0L*M_PI*1.0e6/VLIGHT);        /* Wavenumber in units of m^-1 */
  phases = calloc(stn->layout->nant,sizeof(complex));
  if (phases==NULL) return -1;
  wtsum[0] = wtsum[1]=0.0;

  for (k=0; k<stn->layout->nant; k++){
    /* dot product of the baseline vector with the 
     * difference between the phase centre and the
     * look direction unit vectors */
    trigarg = -wav* (stn->layout->ants[k].east*(dx) +
                     stn->layout->ants[k].north*(dy) +
                     stn->layout->ants[k].height*(dz));
    /* phase contribution of the antenna in the phased array */
    phases[k] = c_exp(trigarg);

    /* accumulate the total weight for each receptor over the station */
    /* FIXME: summing gains is probably not the right thing to do. 
     * These weights should really reflect relative antenna
     * efficiency (e.g. A_eff/T_sys or something) rather than the gain. 
     * If all antennas in a station are the same type,
     * then a good approximation is just the total number of antennas 
     * in the station */
    /*
    wtsum[0] += stn->layout->ants[i].gain[0];
    wtsum[1] += stn->layout->ants[i].gain[1];
    */
    wtsum[0]++;
    wtsum[1]++;
  }

  /* protect against div by zero below */
  if (wtsum[0] ==0.0) wtsum[0] = 1.0;
  if (wtsum[1] ==0.0) wtsum[1] = 1.0;

  /* Loop over receptor, sky polarisation vector combinations to 
   * calculate Jones matrix components.
   * For one pol, we just receive total intensity. For two pols, 
   * we receive sky e-field x and y for each receptor, p and q. 
   * So the order of the index is px, py, qx, qy */
  for (receptor=0; receptor<n_pols; receptor++) {
    for (skypol=0; skypol<n_pols; skypol++) {

      polindex = receptor*n_pols + skypol;

      if (look_dir.z <= 0) {
        /* this point is below the horizon, so set its response to zero */
        beam[polindex] = c_zero();
      }
      else {
        /* Now compute complex beam at this ra,dec	*/
        ant = c_zero();
        /* loop over individual antennas in a station */
        for (k=0; k<stn->layout->nant; k++){
          complex beam_voltage;
          
          /* calculate the complex voltage primary beam response */
          result = primary_beam_response(stn->layout->ants+k ,receptor,
					 skypol, ref_dir, &look_dir,
                                         VLIGHT/(freq*1e6), &beam_voltage);
          if (result !=0) return result;


          /* add to total for phased-array station  */
          ant = c_add(ant,c_mult(phases[k],beam_voltage));
          if (debug_makebeam) printf("ant %d. polindex: %d, beam voltage" \
				     " (%g,%gi), phase: (%g,%gi)\n", 
				     k, polindex, 
				     c_real(beam_voltage), 
				     c_imag(beam_voltage),
				     c_real(phases[k]),
				     c_imag(phases[k]) );
        }

        /* move data into result including scaling by sum of weights */
        beam[polindex] = s_mult(ant,1.0/wtsum[receptor]);
        if (debug_makebeam) {
          printf("Total real: %g, imag: %g\n", 
		 c_real(beam[polindex]),c_imag(beam[polindex]));
        }
      }
    }
  }
  if (phases !=NULL) free(phases);
  return 0;
}

/* Calculate the complex primary beam response of an antenna for a given
 * look direction and antenna pointing direction. This is the complex
 * voltage pattern. The power pattern gets created when two of these 
 * are multiplied later in compute_beamconvol().
 * This function calculates and uses unit vectors in local topocentric coords. 
 * In all cases, the "x" direction is east, "y" is north, and "z" is up 
 * towards the zenith. */
int primary_beam_response(
	 struct  antenna *ant,    /* the antenna to calculate PB */
	 int     receptor,        /* index of the receptor (0 or 1) */
	 int     skypol_index,    /* x or y pol in look direction (0 or 1) */
	 skypos  *ref_dir,        /* unit vectors in local topocentric coords 
				  * (ENH) to pointing cent */
	 skypos  *look_dir,       /* unit vec to look direction */
	 float   lambda,           /* wavelength of obs in m */
	 complex *response       /* result goes here */
)
{
  int result=0;

  *response = c_zero();

  /* switch between antenna types. All the details have been read 
   * in read_stn_beam(). New antenna types that require parameters will 
   * need to have read_stn_beam also updated to extract the antenna 
   * params for that type */
  switch (ant->type) {
  case ANT_ISOTROPIC:
    *response = c_make(1.0,0.0);
    break;
  case ANT_CROSSED_DIPOLES_ON_GROUNDPLANE:
    {
      /* retrieve parameters: PA, height,  PA2, height2*/
      float pa, height,d_dot_p,dipole_x,dipole_y,gp_normaliser;
      float skypol_x,skypol_y,skypol_z,sd,sh,sl,cd,ch,cl;

      if (receptor==0){
        pa = ant->params[0]*M_PI/180.0;
        height = ant->params[1];
      }
      else {
        pa = ant->params[2]*M_PI/180.0;
        height = ant->params[3];
      }

      /* calculate the component of the Jones matrix. One component will be 
       *returned per call to this function. If R and D are unit vectors in 
       * the -RA and DEC direction for the given look direction
       * and p and q are unit vectors in the direction of the receptors 
       * (North and East in ideal case), then
       * the Jones matrix for the antenna will be
       *  J = | p dot D     p dot R |
       *      | q dot D     q dot R |
       * So here, p will be for receptor 0 and q will be for receptor 1 */

      /* first, calculate the unit vector for R or D in look direction */
      sd = sin(look_dir->dec);
      sh = sin(look_dir->ha);
      sl = sin(look_dir->lat);
      cd = cos(look_dir->dec);
      ch = cos(look_dir->ha);
      cl = cos(look_dir->lat);

      if (skypol_index==0) {
        /* unit vec in DEC direction */
        skypol_x = sd*sh;
        skypol_y = cl*cd + sl*sd*ch;
        skypol_z = sl*cd - cl*sd*ch;
      }
      else {
        /* unit vec in -RA direction */
        skypol_x = ch;
        skypol_y = -sl*sh;
        skypol_z = cl*sh;
      }

      /* dipole unit vec (z=0 for horizontal dipoles) */
      dipole_y = cos(pa);
      dipole_x = sin(pa);

      /* projection of vectors (dot prod) gives response 
       * of this dipole to the pol direction */
      d_dot_p = skypol_x*dipole_x + skypol_y*dipole_y;

      gp_normaliser = end_fire_2element(2*height/lambda,0.0);
      *response = c_make(d_dot_p*end_fire_2element(2*height/lambda,
		           M_PI/2.0 - look_dir->alt)/gp_normaliser, 0.0);

      /*      printf("receptor: %d, pa: %f skypol: %d, az: %g, " \ 
       *      "pi/2-alt: %g, endfire: %g, dip(%g,%g,%g), pol(%g,%g,%g), " \
       *      "projection: %g\n",
       *      receptor,pa,skypol_index,look_dir->az,M_PI/2.0-look_dir->alt,
       *     end_fire_2element(2*height/lambda,M_PI/2.0-look_dir->alt),
       *     dipole_x,dipole_y,0.0,skypol_x,skypol_y,skypol_z,
       *     d_dot_p); */
      break;
    }

  case ANT_CROSSED_DIPOLES_HORIZONTAL:
    {
      /* retrieve parameters: PA, PA2 */
      double pa, d_dot_p,dipole_x,dipole_y;
      double skypol_x,skypol_y,skypol_z,sd,sh,sl,cd,ch,cl;

      if (receptor==0){
        pa = ant->params[0]*(M_PI/180.0);
      }
      else {
        pa = ant->params[1]*(M_PI/180.0);
      }

      /* calculate the component of the Jones matrix. One component 
       * will be returned per call to this function. If R and D are 
       * unit vectors in the -RA and DEC direction for the given look 
       * direction  and p and q are unit vectors in the direction of 
       * the receptors (North and East in ideal case), then
       * the Jones matrix for the antenna will be
       *  J = | p dot D     p dot R |
       *      | q dot D     q dot R |
       * So here, p will be for receptor 0 and q will be for receptor 1 */

      /* first, calculate the unit vector for R or D in look direction */
      sd = sin(look_dir->dec);
      sh = sin(look_dir->ha);
      sl = sin(look_dir->lat);
      cd = cos(look_dir->dec);
      ch = cos(look_dir->ha);
      cl = cos(look_dir->lat);

      if (skypol_index==0) {
        /* unit vec in DEC direction */
        skypol_x = sd*sh;
        skypol_y = cl*cd + sl*sd*ch;
        skypol_z = sl*cd - cl*sd*ch;
      }
      else {
        /* unit vec in -RA direction */
        skypol_x = ch;
        skypol_y = -sl*sh;
        skypol_z = cl*sh;
      }

      /* dipole unit vec (z=0 for horizontal dipoles) */
      dipole_y = cos(pa);
      dipole_x = sin(pa);

      /* projection of vectors (dot prod) gives response of this dipole 
       * to the pol direction */
      d_dot_p = skypol_x*dipole_x + skypol_y*dipole_y;

      *response = c_make(d_dot_p,0.0);

      break;
    }

  case ANT_SHORT_DIPOLE_ON_GROUNDPLANE: { /* params: PA (E thru N, radian), 
					   * height above GP (m) */
    /* retrieve parameters: PA, height */
    float pa, height,d_dot_s,gp_normaliser;

    pa = ((float *)ant->params)[0]*M_PI/180.0;
    height = ((float *)ant->params)[1];

    /* calculate dot product of dipole unit vec with look direction (z=0) */
    d_dot_s = look_dir->x*cos(pa) + look_dir->y*sin(pa);
    gp_normaliser = end_fire_2element(2*height/lambda,0.0);
    *response = c_make(end_fire_2element(2*height/lambda,
		         M_PI/2.0 - look_dir->alt)/gp_normaliser*
		         sqrt(1.0 - d_dot_s*d_dot_s), 0.0);
    /*
    //	printf("lambda: %g, 2h/l: %g, pi/2-alt: %g, endfire: %g, " \ 
    //         "projection: %g\n", lambda, 
    //         2*height/lambda,M_PI/2.0 - look_dir->alt, 
    //         end_fire_2element(2*height/lambda,M_PI/2.0 - look_dir->alt),
    //         sqrt(1.0 - d_dot_s*d_dot_s));
    */
    break;
  }
  case ANT_IDEAL_PARABOLOID:/* params: scale length */ {
    float diameter,dist;

    diameter = ((float *)ant->params)[0]; /* diameter of dish in meters */
    /* calculate the angular distance between the antenna pointing direction
       and the direction of interest ("look direction") */
    dist = AngularDist(ref_dir,look_dir);
    *response = c_make(jinc(M_PI*diameter/lambda*sin(dist)),0.0);
    break;
  }
  case ANT_GAUSSIAN: /* params: scale length */ {
    float x,dist;

    dist = AngularDist(ref_dir,look_dir);
    x = dist/((float *)ant->params)[0]; /* HWHM of beam in radian */
    *response = c_make(exp((x*x)/(-2.0)),0.0);
    break;
  }
  case ANT_IDEAL_PARABOLOID_DUAL_LINEAR:/* params: PA1, diam1, PA2, diam2 */ {
    float diameter,dist,pa;
    float skypol_x,skypol_y,skypol_z,sd,sh,sl,cd,ch,cl;

    if (receptor==0){
        pa = ant->params[0]*(M_PI/180.0);
        diameter = ant->params[1];
    }
    else {
        pa = ant->params[2]*(M_PI/180.0);
        diameter = ant->params[3];
    }
    // first, calculate the unit vector for -RA or DEC in look direction
    sd = sin(look_dir->dec); cd = cos(look_dir->dec);
    sh = sin(look_dir->ha);  ch = cos(look_dir->ha);
    sl = sin(look_dir->lat); cl = cos(look_dir->lat);

    if (skypol_index==0) {
        // unit vec in DEC direction
        skypol_x = sd*sh;
        skypol_y = cl*cd + sl*sd*ch;
        skypol_z = sl*cd - cl*sd*ch;
    }
    else {
        // unit vec in -RA direction
        skypol_x = ch;
        skypol_y = -sl*sh;
        skypol_z = cl*sh;
    }
    
    // *** NOTE UNFINISHED ***
    // projection of sky onto dipoles needs to be completed.
    
    /* calculate the angular distance between the antenna pointing direction
       and the direction of interest ("look direction") */
    dist = AngularDist(ref_dir,look_dir);
    *response = c_make(jinc(M_PI*diameter/lambda*sin(dist)),0.0);
    break;
  }
  default:
    fprintf(stderr,
	    "ERROR: primary_beam_response(): unknown antenna type %d\n",
	    ant->type);
    result = -1;
    break;
  }

  /* apply the overall antenna gain and phase */
  *response = s_mult(c_mult(*response, c_exp(ant->phase[receptor])),
		     ant->gain[receptor]);

  if (debug_makebeam) {
    printf("receptor: %d, skypol: %d, response: %g,%g gain: %g\n",
	   receptor, skypol_index, c_real(*response), c_imag(*response),
	   ant->gain[receptor]); 
  }

  return result;
}


/******************************
 calculate the response of a 2-element end-fire antenna, equivalent to
   one antenna above a perfect groundplane
   d = distance between antennas (=2* height above GP) in wavelengths
   th = zenith angle (radian)
   Sep 2007: normalise this so that unit response is at zenith.
            this might turn out to be not the right thing to do.
*******************************/
float end_fire_2element(float d, float th) {
  return sin(M_PI*d*cos(th))/sin(M_PI*d);
}


/*******************************
 * calculate the "Jinc" function for the 
 * voltage response of a circular aperture
 ***************************** */
float jinc(float x) {
  if(fabs(x) < 1e-6) x = 1e-6; /* avoid dealing with very small numbers */
  return 2.0*bessj1(x)/x;
}


/*********************************
 * calculate the Bessel function J1
 * rip-off from NumRec
 **********************************/
float bessj1(float x)
/*  Returns the Bessel function J1 (x) for any real x. */
{
  float ax,z;
  /*  Accumulate polynomials in double precision. */
  double xx,y,ans,ans1,ans2;
    /* Direct rational approximation. */
  if ((ax=fabs(x)) < 8.0) {
    y = x*x;
    ans1 = x*(72362614232.0 + y*(-7895059235.0 + y*(242396853.1 + 
	       y*(-2972611.439 + y*(15704.48260 + y*(-30.16036606))))));
    ans2 = 144725228442.0 + y*(2300535178.0 + y*(18583304.74 + 
               y*(99447.43394 + y*(376.9991397 + y*1.0))));
    ans=ans1/ans2;
    
  } 
  else {
    z = 8.0/ax;
    y = z*z;
    xx = ax - 2.356194491;
    ans1 = 1.0 + y*(0.183105e-2 + y*(-0.3516396496e-4 + 
	             y*(0.2457520174e-5 + y*(-0.240337019e-6))));
    ans2 = 0.04687499995 + y*(-0.2002690873e-3 + y*(0.8449199096e-5 + 
			       y*(-0.88228987e-6 + y*0.105787412e-6)));
    ans = sqrt(0.636619772/ax)*(cos(xx)*ans1 - z*sin(xx)*ans2);
    if (x < 0.0) ans = -ans;
  }
  return ans;
}


/********************************
 * anglular distance in radian between two points on unit sphere
 *******************************/
double AngularDist( skypos *ref_dir, skypos *look_dir
			  /* calculate the angular distance between 
			   * the antenna pointing direction
			   * and the direction of interest 
			   * ("look direction") */
) {
  double dist;

  /* try dot product. Angle is acos(product). 
   * Bad for small angles or product close to 1 */
  dist = look_dir->x*ref_dir->x + look_dir->y*ref_dir->y + 
    look_dir->z*ref_dir->z;
  if (fabs(dist) < 0.98) {
    dist = acos(dist);
  }
  else {
    /* use the cross product. Angle is asin(mag of vector). 
     * Bad for angles close to pi/2 */
    double x,y,z;

    x = ref_dir->y*look_dir->z - ref_dir->z*look_dir->y;
    y = ref_dir->z*look_dir->x - ref_dir->x*look_dir->z;
    z = ref_dir->x*look_dir->y - ref_dir->y*look_dir->x;
    dist = sqrt(x*x + y*y + z*z);
    dist = asin(dist);
  }
  return dist;
}


/*****************************
 * Convert an HA, DEC and station latitude (all radian) 
 * into a unit vector on the unit sphere in local
 * topocentric coords. x== east, y==north, z==up. 
 * Calculate the alt and az of the vector too.
******************************/
void MakeSkyPos(double ha, double dec, double lat, skypos *result) {

  double sl,cl,sh,ch,sd,cd;

  cl = cos(lat);
  sl = sin(lat);
  cd = cos(dec);
  sd = sin(dec);
  ch = cos(ha);
  sh = sin(ha);

  result->ha = ha;
  result->dec = dec;
                   /* x,y,z location on unit sphere (i.e. unit vector) 
		    * in local topocentric coords (x==east,y==north,z==up) 
		    * of the phase reference of the phased-array beam
		    * (i.e. pointing centre) */
  result->x = -sh*cd;
  result->y = sd*cl - ch*cd*sl;
  result->z = ch*cd*cl + sd*sl;

  /* calc alt,az for later reference */
  result->az  = atan2(result->x, result->y);
  result->alt = asin(result->z);
  result->lat = lat;
}


/*********************************
 * Dump the station beam (all pols) to a FITS file. 
 * This is for debugging and testing.
 * The beams themselves are complex, but FITS doesn't 
 * currently support complex, so just write
 * out the real part.
 *********************************/
int DumpStationbeamToFITS(struct   stnbeam *beam,
                            struct stationspec *stn,
                            int    npols,
                            double freq,
                            double gha,
                            double dec_ref,
                            double ra_ref,
                            float  crval_ra, float crval_dec,
                            float  cdelt_ra, float cdelt_dec) {
    char filename[FILENAME_MAX];
    fitsfile *fp=NULL;
    int status=0,pol=0,naxis=2,i;
    float *data=NULL,temp;
    long naxes[3]={0,0,0},fpixel[3]={1,1,1};
    
    /* create the filename and open the file */
    sprintf(filename,"st_beam_%03d.fits",stn->st_id);
    remove(filename); /* just overwrite any existing file */
    fits_create_file(&fp, filename, &status);
    if (status!=0) {
        msg("DumpStationbeamToFITS: failed to create file <%s>",0,filename);
        goto EXIT;
    }

    /* if there is only 1 pol, then create a 2D image, 
     * otherwise create a 3D image */
    if (npols > 1) {
        naxis = 3;
        naxes[2] = npols;
    }
    
    /* set up the image axis sizes and allocate temp space */
    naxes[0] = naxes[1] = beam->size;
    data = malloc(sizeof(float)*naxes[0]*naxes[1]);

    /* Write the main data as single precision float. 
     * Currently the beams are double precision complex. */
    fits_create_img( fp, FLOAT_IMG, naxis, naxes, &status);
    if (status!=0) {
        msg("DumpStationbeamToFITS: status %d on fits_create_img",0,status);
        goto EXIT;
    }
    for (pol=0; pol<npols; pol++) {
        fpixel[2] = pol+1;
        /* copy the beam mag into a temp array */
	/* for(i=0; i<naxes[0]*naxes[1]; i++) 
	 *   data[i] = c_mag(beam->beam[pol][i]); */
        for(i=0; i<naxes[0]*naxes[1]; i++) 
	  data[i] = c_real(beam->beam[pol][i]);
        /*write out array for this pol*/
        fits_write_pix(fp, TFLOAT, fpixel,naxes[0]*naxes[1] ,data, &status);
    }
    
    /* write various metadata to the file */
    fits_update_key(fp, TDOUBLE, "FREQ" , &freq, 
		    "Obs freq in MHz", &status);
    fits_update_key(fp, TDOUBLE, "GHA" , &gha, 
		    "Greenwhich HA in radian", &status);
    fits_update_key(fp, TDOUBLE, "REF_RA" , &ra_ref, 
		    "Pointing reference RA in radian", &status);
    fits_update_key(fp, TDOUBLE, "REF_DEC" , &dec_ref, 
		    "Pointing reference DEC in radian", &status);
    fits_update_key(fp, TINT,    "ST_ID" , &(stn->st_id), 
		    "station ID (starts from zero)", &status);
    fits_update_key(fp, TSTRING, "CTYPE1", "RA---SIN",
		    "Simple WCS. Sine projection",&status);
    fits_update_key(fp, TSTRING, "CTYPE2", 
		    "DEC--SIN","Simple WCS",&status);
    temp = beam->xref_pixel+1; /* FITS pixels start at 1 */
    fits_update_key(fp, TFLOAT,  "CRPIX1" , &temp, 
		    "Simple WCS reference pixel - RA", &status);
    temp = beam->yref_pixel+1;
    fits_update_key(fp, TFLOAT,  "CRPIX2" , &temp, 
		    "Simple WCS reference pixel - DEC", &status);
    fits_update_key(fp, TFLOAT,  "CRVAL1" , &crval_ra, 
		    "RA at reference pixel", &status);
    fits_update_key(fp, TFLOAT,  "CRVAL2" , &crval_dec, 
		    "DEC at reference pixel", &status);
    fits_update_key(fp, TFLOAT,  "CDELT1" , &cdelt_ra, 
		    "angular size of pixels - RA", &status);
    fits_update_key(fp, TFLOAT,  "CDELT2" , &cdelt_dec, 
		    "angular size of pixels - DEC", &status);

    /* clean up */
    fits_close_file(fp, &status);
    if (status!=0) {
        msg("DumpStationbeamToFITS: status %d on closefile",1,status);
        goto EXIT;
    }
    msg("wrote %dx%dx%d beam array for station %d to %s",
	-1, (int)naxes[0], (int)naxes[1], (int)naxes[2], stn->st_id, filename);
 EXIT:
    if(data!=NULL) free(data);
    return status;
}

