/************************************************************************/
/*                                                                      */
/* Given the station voltage beams at each end of the baseline, this    */
/* routine figures out what the convolution function in the uv plane    */
/* needs to be in order to properly modify the visibilities to simulate */
/* the effects of the beams                                             */
/*                                                                      */
/*      Inputs:     beam1, beam2    Filled station beam structs         */
/*                                                                      */
/*      Output:     beamconvl       Derived convolution function        */
/*                  return value    0=OK, else bad                      */
/*                                                                      */
/* Created 18 june 2002 by SSD						*/
/*                                                                      */
/************************************************************************/
#include <stdlib.h>
#include "fitsio.h"
#include "station_beam.h"
#include "convolve.h"
#include "comp_func.h"
#include "fourc.h"
#include "utility.h"
#include "ionosphere.h"
#include "compute_beamconvl.h"

#define NDIM 2 	/* Dimension of FFT */

/* private function prototypes */
int DumpComplexArrayToFITS(complex *beam, char *base_filename, 
			   int sizex, int sizey);
int DumpFloatArrayToFITS(float *beam,char *base_filename, 
			 int sizex, int sizey);
                    
static const int m_conv[16] = {0,1,4,5,2,3,6,7,8,9,12,13,10,11,14,15};
static int debug=0;


void compute_beamconvl_set_debug(int level) {
    debug = level;
}


int
compute_beamconvl(/* stbeams, st1, st2, ionoscreen, freq, do_ion, beamconvl) */
		  struct stnbeam *stbeams,
		  int st1_index,
		  int st2_index,
		  struct iono_screen *ionoscreen,
		  double freq,
		  int do_ion,
		  struct conv_fn *beamconvl) {

  int i,j,k,l,ni,nj,isign,nn[2],num_pixels;
  int nindex, index,m_index;
  complex *temp=NULL;
  struct stnbeam *beam2,*beam1;

  beam1 = stbeams + st1_index;
  beam2 = stbeams + st2_index;

  /* Check to make sure both beams have same dimension	*/
  if (beam1->size!=beam2->size) 	{
    msg ("compute_beamconvl(): station beam arrays in baseline pair " \
	 "have different sizes",2);
    exit(1);
  }
  num_pixels = (beam1->size)*(beam1->size);

  /* make working space */
  temp = calloc(num_pixels,sizeof(complex));
  if (temp ==NULL) {
    msg("compute_beamconvl: no malloc",3);
    return -1;
  }

  /* Mueller matrix must look like:
     M= | p1x1*p2x2 | p1x1*p2y2 | p1y1*p2x2 | p1y1*p2y2 |
        | p1x1*q2x2 | p1x1*q2y2 | p1y1*q2x2 | p1y1*q2y2 |
        | q1x1*p2x2 | q1x1*p2y2 | q1y1*p2x2 | q1y1*p2y2 |
        | q1x1*q2x2 | q1x1*q2y2 | q1y1*q2x2 | q1y1*q2y2 |

        ***NOTE***
     There is similar code to this function in pol_response.c 
     These codes MUST be kept in sync if changes
     are made.

  */

  /* loop over each Mueller matrix element */
  for (k=0; k< 4; k++) {  /* beam1 pols p1x1,p1y1,q1x1,q1y1 */
    for (l=0; l<4; l++) { /* beam2 pols p2x2,p2y2,q2x2,q2y2 */

      /* index into the Mueller matrix. The structure of the  *
       * Mueller matrix doesn't lend itself nicely to two     *
       * loops across the beam products, so we use m_conv to  *
       * decode the appropriate index into the Mueller matrix *
       * given the beam pols                                  */
      m_index = m_conv[k*4 +l];

      /* don't do anything if there are no beams to multiply  *
       * (unpolarised case)                                   */
      if (beam1->beam[k]==NULL || beam2->beam[l]==NULL) continue;

      /* make space for result if it is required and if it    *
       * hasn't been allocated already                        */
      if (beamconvl->convl[m_index] ==NULL) {
		beamconvl->convl[m_index] = calloc(num_pixels,sizeof(complex));
        if (beamconvl->convl[m_index] ==NULL) {
          msg("compute_beamconvl: no malloc",3);
          return -1;
        }
        msg("alloced %d complex for Mueller matrix index %d", 0,
	    num_pixels, m_index);
      }

      /* Complex multiply beams together. Complex conjugate second beam. */
      for (i=0; i<beam1->size; i++) {
        for(j=0; j<beam1->size; j++) {
          index = i*beam1->size+j;
          beamconvl->convl[m_index][index] = 
	    c_mult(beam1->beam[k][index],c_conj(beam2->beam[l][index]));
          if (do_ion) {
            double phase_diff;
            /* NOTE: the mathematics says that the phase should be *
	     * negative and we should be subtracting station 1     *
             * from station 2. However, the sources move in the    *
	     * wrong direction, so there is a negative sign        *
	     * problem somewhere else in the software. We add the  *
	     * extra negative here to make it correct, but this    *
	     * might need to be updated if negative signs are      *
	     * fixed elsewhere.                                    */
            phase_diff = -(ionoscreen->phase_delay[st2_index][index] - 
			   ionoscreen->phase_delay[st1_index][index])/freq;
            msg("applying phase %g\tradian to beam grid "\
		"pixel %d,%d for m_index %d", -1, phase_diff, i, j, m_index);
            beamconvl->convl[m_index][index] = 
	      c_mult(beamconvl->convl[m_index][index], rect(1.0, phase_diff));
          }
        }
      }

      /* Print out station beam product when debugging */
      if (debug && beam1->station==0 && beam2->station==1) {
        DumpComplexArrayToFITS(beamconvl->convl[m_index], "beamconv_pre", 
			       beam1->size, beam1->size);
        if (do_ion) 
	  DumpFloatArrayToFITS(ionoscreen->phase_delay[st1_index],
			       "ionoscreen", beam1->size, beam1->size);
      }
      
      /* Now shift center of beam to (0,0) pixel in 	*/
      /* preparation for FFT.				*/
      for (i=0;i<beam1->size;i++) {
        if (i-beam1->yref_pixel < 0) ni = (i-beam1->yref_pixel) + beam1->size;
        else ni = (i-beam1->yref_pixel);
        for (j=0;j<beam1->size;j++) {
          if (j-beam1->xref_pixel < 0) 
	    nj = (j-beam1->xref_pixel) + beam1->size;
          else nj = (j-beam1->xref_pixel);
          
          index = beam1->size*i + j;
          nindex = beam1->size*ni + nj;
          
          temp[index] = beamconvl->convl[m_index][nindex];
        }
      }

      /* Print out shifted station beam product */
      if (debug && beam1->station==0 && beam2->station==1) {
        DumpComplexArrayToFITS(temp, "beamconv_post", beam1->size, beam1->size);
      }

      /* Now FFT beam product	*/

      /* Have to setup size array for FFT routine 	*/
      /* isign=1 gives FFT, isign=-1 gives inverse 	*/
      nn[0]=beam1->size;
      nn[1]=beam2->size;
      isign=-1;
	     	/* Make inverse inplace FFT (isign == -1) on temp    */ 
      fourc(temp,nn,NDIM,isign);

      /* Print out fft'd station beam product */
      if (debug  && beam1->station==0 && beam2->station==1) {
        DumpComplexArrayToFITS(temp, "beamconv_fft_noshift", 
			       beam1->size, beam1->size);
      }

      /* Now shift center of convolving function to center pixel in   */
      /* preparation for convolution.		*/
      for (i=0;i<beam1->size;i++){
        if (i+beam1->yref_pixel < beam1->size) ni = (i+beam1->yref_pixel);
        else ni = (i+beam1->yref_pixel)-beam1->size;
        for (j=0;j<beam1->size;j++) {
          if (j+beam1->xref_pixel < beam1->size) nj = (j+beam1->xref_pixel);
          else nj = (j+beam1->xref_pixel) - beam1->size;
          
          index = beam1->size*i + j;
          nindex = beam1->size*ni + nj;
          
          beamconvl->convl[m_index][index]=temp[nindex];
        }
      }
      
      /* Print out shifted fft'd station beam product */
      if (debug  && beam1->station==0 && beam2->station==1) {
        DumpComplexArrayToFITS(temp, "beamconv_fft_shifted", 
			       beam1->size, beam1->size);
      }

    }
  }
  
  beamconvl->xref_pixel = beam1->size/2;
  beamconvl->yref_pixel = beam1->size/2;
  beamconvl->xsize = beam1->size;
  beamconvl->ysize = beam1->size;
  if (temp !=NULL) free(temp);
  return (0);
}


/*********************************
 Dump a 2D complex array to a FITS file. This is for debugging and testing.
 Could be used for station beams, convolution functions and the like.
 FITS does not support complex intrinsically, so write out the real and 
 imag parts in 2 separate arrays.
 *********************************/
int DumpComplexArrayToFITS(complex *beam,
                    char *base_filename,
                    int sizex, int sizey)
{
    char filename[FILENAME_MAX];
    fitsfile *fp=NULL;
    int status=0,naxis=3,i;
    float *data=NULL;
    long naxes[3]={0,0,2},fpixel[3]={1,1,1};

    /* create the filename and open the file */
    sprintf(filename,"%s.fits",base_filename);
    remove(filename); /* just overwrite any existing file */
    fits_create_file(&fp, filename, &status);
    if (status!=0) {
        msg("DumpComplexArrayToFITS: failed to create file <%s>", 0, filename);
        goto EXIT;
    }

    /* set up the image axis sizes and allocate temp space */
    naxes[0] = sizex;
    naxes[1] = sizey;
    data = malloc(sizeof(float)*naxes[0]*naxes[1]*naxes[2]);

    /* write the main data as single precision float. Currently   *
     * the beams are double precision complex.                    */
    fits_create_img( fp, FLOAT_IMG, naxis, naxes, &status);
    if (status!=0) {
        msg("DumpComplexArrayToFITS: status %d on fits_create_img", 0, status);
        goto EXIT;
    }
    /* copy the beam mag into a temp array        */
    /*    for(i=0;i<naxes[0]*naxes[1]; i++)       */
    /*      data[i] = c_mag(beam->beam[pol][i]);  */
    for(i=0;i<naxes[1]*naxes[0]; i++) {
        data[i] = c_real(beam[i]);
        data[i+naxes[1]*naxes[0]] = c_imag(beam[i]);
    }
    /*write out array for this pol*/
    fits_write_pix(fp, TFLOAT, fpixel, naxes[0]*naxes[1]*naxes[2], data, 
		   &status);

    /* clean up */
    fits_close_file(fp, &status);
    if (status!=0) {
        msg("DumpComplexArrayToFITS: status %d on closefile",1,status);
        goto EXIT;
    }
    msg("wrote %dx%dx%d conv beam array to %s", 1, (int)naxes[0],
	(int)naxes[1], (int)naxes[2], filename);
EXIT:
    if(data!=NULL) free(data);
    return status;
}


/*********************************
 Dump a 2D float array to a FITS file. This is for debugging and testing.
 Could be used for station beams, convolution functions and the like.
*******************************/
int DumpFloatArrayToFITS(float *beam,
                    char *base_filename,
                    int sizex, int sizey)
{
    char filename[FILENAME_MAX];
    fitsfile *fp=NULL;
    int status=0,naxis=2;
    long naxes[2]={0,0},fpixel[2]={1,1};

    /* create the filename and open the file */
    sprintf(filename,"%s.fits",base_filename);
    remove(filename); /* just overwrite any existing file */
    fits_create_file(&fp, filename, &status);
    if (status!=0) {
        msg("DumpFloatArrayToFITS: failed to create file <%s>",0,filename);
        goto EXIT;
    }

    /* set up the image axis sizes and allocate temp space */
    naxes[0] = sizex;
    naxes[1] = sizey;

    /* write the main data as single precision float.      *
     * Currently the beams are double precision complex.   */
    fits_create_img( fp, FLOAT_IMG, naxis, naxes, &status);
    if (status!=0) {
        msg("DumpFloatArrayToFITS: status %d on fits_create_img",0,status);
        goto EXIT;
    }

    /*write out array */
    fits_write_pix(fp, TFLOAT, fpixel,naxes[0]*naxes[1] ,beam, &status);

    /* clean up */
    fits_close_file(fp, &status);
    if (status!=0) {
        msg("DumpFloatArrayToFITS: status %d on closefile",1,status);
        goto EXIT;
    }
    msg("wrote %dx%d float array to %s", 1, (int)naxes[0], (int)naxes[1], 
	filename);
EXIT:
    return status;
}

