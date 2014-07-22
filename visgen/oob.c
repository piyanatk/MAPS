/************************************************************************
* Functions to read in from file and add to visibilities a list of point sources
*
* Created July 2007 by Randall Wayth
*
************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "utility.h"
#include "type_comp.h"
#include "comp_func.h"
#include "pol_response.h"
#include "makebeam.h"
#include "ionosphere.h"
#include "coord_trans.h"
#include "oob.h"

static int debug_oob=0;

int oob_beam_response(struct observing_param *obsparam, struct array_specification *arrayspec, int st_index,
                      double gha, double freq, source_model_t *oobs, int num_src, oob_beam_t *beams );
int oob_calc_geom_phase(struct observing_param *obsparam, struct array_specification *arrayspec,double gha,double freq, int st1_index,int st2_index,source_model_t *oob, complex *phase);

/****************************
****************************/
void oob_set_debug(int level) {
  debug_oob = level;
}


/****************************
  read a text file with positions and fluxes of point sources into an array of source_model_t structs.
  File format: RA (decimal hours), DEC (decimal degs) I,Q,U,V (stokes)
  blank lines and lines beginning with "#" are ignored
****************************/
int read_oob(char *filename, source_model_t **ooblist, int *noob) {
  FILE *fp=NULL;
  char line[BUFSIZ];
  source_model_t *list=NULL;
  int result=0;

  *noob=0;
  
  /* open file */
  if ((fp=fopen(filename,"r"))==NULL) {
    msg("read_oob: cannot open file <%s>",2,filename);
    return -1;
  }


  /* process each row of the file */
  while(fgets(line,BUFSIZ,fp) != NULL) {
    /* ignore blanks and comments */
    rm_whitespace(line);
    if (line[0] == '\0' || line[0] =='\n' || line[0] == '#') continue;

    /* make space for new item */
    list = realloc(list,sizeof(source_model_t)*(*noob+1));
    if (list==NULL) {
      msg("read_oob: realloc failed after %d sources",3,*noob);
      return -1;
    }
    /* line format: ra, dec, I, Q, U ,V */
    result = sscanf(line,"%lf %lf %f %f %f %f",&(list[*noob].ra),&(list[*noob].dec),
                    &(list[*noob].flux[0]),&(list[*noob].flux[1]),
                    &(list[*noob].flux[2]),&(list[*noob].flux[3]));
    if (result != 6) {
      msg("read_oob: failed to scan line: <%s>",3,line);
      return 1;
    }
    msg("oob source: (RA,DEC): %g,%g (I,Q,U,V): %g,%g,%g,%g",0,list[*noob].ra,list[*noob].dec,
        list[*noob].flux[0],list[*noob].flux[1],list[*noob].flux[2],list[*noob].flux[3]);

    (*noob)++;
  }
 
  *ooblist = list;
  fclose(fp);
  return 0;
}


/****************************
 add the flux from OOB point sources into visibilities.

 this function calculates the station beam response for each station for each OOB point source. The results are
 stored along the way in the "beams" array to prevent re-calculaing station responses for each baseline that a
 station forms. The flux from each source is added into "vis".

 gha: the greenwich HA of the correlator phase center (radian)
 freq: observing freq (MHz)
 st_index: integer index for the stations forming a baseline (0 .. (n_stations-1))
 oobs: array of OOB sources previously read in
 num_src: size of oobs array
 beams: array the size of n_antennas which store the station beam responses for each station, for each source
****************************/
int oob_integrate(struct observing_param *obsparam, struct array_specification *arrayspec,double gha, double freq,
                  int st1_index,int st2_index,source_model_t *oobs, int num_src, oob_beam_t *beams,
                  complex vis[4]) {
  int src,pol,result,num_pol_prod,time_ind,freq_ind;
  complex mueller[16],m_s[16],j1[4],j2[4],temp[4],phase,temp_phase;

  num_pol_prod = arrayspec->num_pols*arrayspec->num_pols;
  temp[0] = temp[1] = temp[2] = temp[3] = c_zero();

  /* if either of the beam Jones matrices haven't been calculated yet, do it now */
  if (beams[st1_index].valid ==0) {
    result = oob_beam_response(obsparam,arrayspec,st1_index,gha,freq,oobs,num_src,beams);
    if (debug_oob) {
      int i;

      fprintf(stdout,"Jones matrix for station %d (%s) source 0: ",st1_index,arrayspec->station[st1_index].layout->layout_name);
      for (i=0; i<num_pol_prod; i++) {
        fprintf(stdout,"(%g,%g) ",c_real(beams[st1_index].response[i][0]),c_imag(beams[st1_index].response[i][0]));
      }
      fprintf(stdout,"\n");
    }
  }
  if (beams[st2_index].valid ==0) {
    result = oob_beam_response(obsparam,arrayspec,st2_index,gha,freq,oobs,num_src,beams);
    if (debug_oob) {
      int i;

      fprintf(stdout,"Jones matrix for station %d (%s) source 0: ",st2_index,arrayspec->station[st2_index].layout->layout_name);
      for (i=0; i<num_pol_prod; i++) {
        fprintf(stdout,"(%g,%g) ", c_real(beams[st2_index].response[i][0]),c_imag(beams[st2_index].response[i][0]) );
      }
      fprintf(stdout,"\n");
    }
  }

  /* loop over each source and accumulate the contribution to the visibility */
  for(src=0; src < num_src; src++) {
    int n_time, n_freq;
    /* extract Jones matrix for each station */
    for (pol=0; pol<num_pol_prod; pol++) {
      j1[pol] = beams[st1_index].response[pol][src];
      j2[pol] = beams[st2_index].response[pol][src];
    }

    /* calculate the Mueller matrix */
    pol_make_mueller(num_pol_prod,j1,j2, mueller);

    if (debug_oob) {
      fprintf(stdout,"Mueller matrix for baseline %d-%d source %d, freq: %g:\n",st1_index,st2_index,src,freq);
      print_complex_matrix(num_pol_prod,mueller);
    }
    /* calculate M*S*I where M is the Mueller matrix, S is the Stokes converter and I is the input Stokes. */
    pol_make_MS(num_pol_prod,mueller,stokes_2_receptor, m_s);
    if (debug_oob) {
      fprintf(stdout,"M*S matrix for baseline %d-%d source %d:\n",st1_index,st2_index,src);
      print_complex_matrix(num_pol_prod,m_s);
    }
    pol_make_MI(num_pol_prod,m_s, oobs[src].flux, temp);
    if (debug_oob) {
      fprintf(stdout,"MI: (%g,%gi),(%g,%gi),(%g,%gi),(%g,%gi)\n",c_real(temp[0]),c_imag(temp[0]),c_real(temp[1]),c_imag(temp[1]),
              c_real(temp[2]),c_imag(temp[2]),c_real(temp[3]),c_imag(temp[3]) );
    }
    
    /* calculate the geometric phase */
    /* here we do sub time and freq averaging to account for the bandwidth pattern
       and the fringe rate. Divide the correlator integration time and channel up
       into a number of chunks as defined in the obsparam and average */
    phase = c_zero();
    n_time = (obsparam->nt_cells==0? 1:obsparam->nt_cells);
    n_freq = (obsparam->nf_cells==0? 1:obsparam->nf_cells);
    for(time_ind=0; time_ind < n_time; time_ind++) {
        for(freq_ind=0; freq_ind < n_freq; freq_ind++) {
            oob_calc_geom_phase(obsparam, arrayspec, gha+time_ind*obsparam->integ_time/n_time/3600.0*M_PI/12.0*(24.0/23.9344696),
                        freq+freq_ind*obsparam->corr_chan_bw/n_freq, st1_index, st2_index, oobs+src, &temp_phase);
            phase = c_add(phase,temp_phase);
        }
    }
    /* now get average */
    phase = s_mult(phase,1./(n_time*n_freq));
    if (debug_oob) fprintf(stdout,"Integrated phase: (%g,%g)\n",c_real(phase),c_imag(phase));

    /* apply the ionospheric phase */
    if(beams[st1_index].valid_phase && beams[st2_index].valid_phase) {
        complex iono_phase = rect(1.0,(beams[st2_index].phase[src] - beams[st1_index].phase[src])/freq);
        phase = c_mult(phase,iono_phase);
    }

    /* accumulate visibilities */
    for (pol=0; pol<num_pol_prod; pol++) {
#ifdef _COMPLEX_H
        vis[pol] += temp[pol]*phase;
#else
        /* old style manual complex */
      vis[pol].re += temp[pol].re*phase.re - temp[pol].im*phase.im;
      vis[pol].im += temp[pol].im*phase.re + temp[pol].re*phase.im;
#endif
    }
    if(debug_oob) {
      printf("Vis: (%g,%g) (%g,%g) (%g,%g) (%g,%g)\n",c_real(vis[0]),c_imag(vis[0]),c_real(vis[1]),c_imag(vis[1])
                                                     ,c_real(vis[2]),c_imag(vis[2]),c_real(vis[3]),c_imag(vis[3]));
    }
  }

  return 0;
}


/****************************
  calculate the geometric phase of a source at position "look_dir" vs the phase reference position "ref_dir".
  Use a generic, vector based approach that does not make assumptions about flat sky or array config.
  All units are in conventional radio astronomy X,Y,Z units with X pointing to the Greenwich meridian, Y east
  and Z to north pole. Follows from Thompson, Moran & Swenson, chapter 3.
****************************/
int oob_calc_geom_phase(struct observing_param *obsparam, struct array_specification *arrayspec,double gha,
                        double freq, int st1_index,int st2_index,source_model_t *oob, complex *phase) {

  double ref_vec[3],lk_vec[3]; /* unit vectors in X,Y,Z coord system */
  double bl_X,bl_Y,bl_Z,d_dot_s=0,look_gha,look_dec;

  /* calculate unit vector to phase reference position */
  sperical2DirectionCosines( -gha, obsparam->phase_cent_DEC, ref_vec);

  /* calculate X,Y,Z of baseline (X points to Greenwich meridian) in meters. */
  /* NOTE: definition of baseline must be ant1-ant2 to have negative exponent in exp below*/
  bl_X = arrayspec->station[st1_index].x_coord - arrayspec->station[st2_index].x_coord;
  bl_Y = arrayspec->station[st1_index].y_coord - arrayspec->station[st2_index].y_coord;
  bl_Z = arrayspec->station[st1_index].z_coord - arrayspec->station[st2_index].z_coord;

  /* now unit vec to look direction */
  look_dec = oob->dec*(M_PI/180.0); /* decimal degs to radian */
  look_gha = gha + obsparam->phase_cent_RA - (oob->ra*(M_PI/12.0));
  sperical2DirectionCosines(-look_gha, look_dec, lk_vec);

  /* dot product baseline (in wavelengths) with (look-ref) vec */
  d_dot_s = ((lk_vec[0]-ref_vec[0])*bl_X+(lk_vec[1]-ref_vec[1])*bl_Y+(lk_vec[2]-ref_vec[2])*bl_Z)*freq*(1e6/VLIGHT);

  /* calculate phase, negative sign for Fourier transform should be in exponent by convention. */
  *phase = c_exp(-2.0*M_PI*d_dot_s);

  if (debug_oob) {
    fprintf(stdout,"ref gha/dec: %g,%g. look gha/dec: %g,%g\n",gha,obsparam->phase_cent_DEC,look_gha,look_dec);
    fprintf(stdout,"station 1: x,y,z: %.10g,%.10g,%.10g\n",arrayspec->station[st1_index].x_coord,
	    arrayspec->station[st1_index].y_coord,
	    arrayspec->station[st1_index].z_coord);
    fprintf(stdout,"station 2: x,y,z: %.10g,%.10g,%.10g\n",arrayspec->station[st2_index].x_coord,
	    arrayspec->station[st2_index].y_coord,
	    arrayspec->station[st2_index].z_coord);
	    
    fprintf(stdout,"phaseref vec: (%g,%g,%g)\n",ref_vec[0],ref_vec[1],ref_vec[2]);
    fprintf(stdout,"source   vec: (%g,%g,%g)\n",lk_vec[0],lk_vec[1],lk_vec[2]);
    fprintf(stdout,"skypos   vec: (%g,%g,%g)\n",lk_vec[0]-ref_vec[0],lk_vec[1]-ref_vec[1],lk_vec[2]-ref_vec[2]);
    fprintf(stdout,"baseline vec: (%g,%g,%g)\n",bl_X,bl_Y,bl_Z);
    fprintf(stdout,"product: %g (wavelenghts), phase: (%g,%g)\n",d_dot_s,c_real(*phase),c_imag(*phase));
  }

  return 0;
}



/****************************
 calculate the station beam response for a single station (st_index) for all the directions of the OOB sources.
****************************/
int oob_beam_response(struct observing_param *obsparam, struct array_specification *arrayspec, int st_index,
                      double gha, double freq,source_model_t *oobs, int num_src, oob_beam_t *beams ) {

  skypos ref_dir;
  double hang,ra,dec;
  int pol_index,src_index,receptor,skypol,result;
  complex beam[4];
  
  beam[0] = beam[1] = beam[2] = beam[3] = c_zero();

  /* set up the reference dir for the pointing center of the station. Note that may be different from the phase
     center of the array (correlator). Each station will have a (slightly) different location and since the stations
     follow the curvature of the Earth, have a slightly different ha */
  /* The GHA is for the correlator phase center. We need to know the HA of the pointing center.
     The HA(point) = LST - RA(point). We know the RA of the pointing center and correlator phase center.
     LST = GHA(correlator) + lon + RA(correlator), so we use this to calculate the HA of the pointing center below */
  MakeSkyPos(gha + arrayspec->station[st_index].longitude + obsparam->phase_cent_RA - obsparam->point_cent_RA,
            obsparam->point_cent_DEC, arrayspec->station[st_index].latitude, &ref_dir);

  for (src_index=0; src_index<num_src; src_index++) {
    ra = oobs[src_index].ra*(M_PI/12.0);
    dec= oobs[src_index].dec*(M_PI/180.0);
    hang = gha + arrayspec->station[st_index].longitude + (obsparam->phase_cent_RA - ra);

    result = station_response(arrayspec->station+st_index, &ref_dir, hang, dec, freq, beam);
    if (result != 0) return result;

    /* loop through receptor and sky pol combinations. The order of the loops here is important, following
       makebeam.c. The receptor loop must be outer. */
    for(receptor=0; receptor < arrayspec->num_pols; receptor++) {
      for(skypol=0; skypol < arrayspec->num_pols; skypol++) {
        pol_index = receptor*arrayspec->num_pols + skypol;

        beams[st_index].response[pol_index][src_index] = beam[pol_index];
      }
    }
    if(debug_oob) {
      fprintf(stdout,"src RA/DEC: %g,%g voltage: (%g,%g),(%g,%g),(%g,%g),(%g,%g)\n",ra,dec,
              c_real(beam[0]),c_imag(beam[0]),c_real(beam[1]),c_imag(beam[1]),c_real(beam[2]),c_imag(beam[2]),c_real(beam[3]),c_imag(beam[3]) );
    }
  }
  beams[st_index].valid = 1;
  return 0;
}


/************************
************************/
void init_oob(const int n_oob, const int n_stations, const int n_pol_prodcuts, oob_beam_t **oob_beams) {

    oob_beam_t *beams=NULL;
    int i,j;

    /* alloc space for the OOB station beams, one for each station */
    beams = (*oob_beams) = calloc(n_stations,sizeof(oob_beam_t));
    if (beams==NULL) {
        msg("no malloc for oob beams",3);
        exit(1);
    }
    
    /* now within each oob station beam, allocate space for each source */
    for (j=0; j< n_stations; j++) {
        for (i=0; i<n_pol_prodcuts;i++){
          beams[j].response[i] = calloc(n_oob,sizeof(complex));
          if (beams[j].response[i] ==NULL) {
            msg("no malloc for oob beams",3);
            exit(1);
          }
        }
        beams[j].phase = malloc(n_oob*sizeof(double));
        if (beams[j].phase ==NULL) {
            msg("no malloc for oob phases",3);
            exit(1);
        }
    }
}

