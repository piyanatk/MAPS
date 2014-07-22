#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include "uvfits.h"

/* make some fake data for a 6 antenna E-W array */

#define NVIS 10      /* number of time instants */
#define NANT 6
#define NPOL 4
#define NFREQ 13
#define NBASE (NANT*(NANT-1)/2)
#define TIMESTEP 10

void calcUVW(float HA, float DEC, float X, float Y, float Z, double *u, double *v, double *w);

int main(int argc, char **argv){

  int i,j,k,l,mon,day,year,siz_visdump=NBASE*NPOL*NFREQ,ant1,ant2;
  double *u[NVIS],*v[NVIS],*w[NVIS],t[NVIS];
  double *u_ptr, *v_ptr, *w_ptr;
  float *vis[NVIS],*wt[NVIS];
  float *vis_ptr, *wt_ptr, ant_y_offsets[NANT]={64,84,100,110,113,392};
  uvdata data;
  ant_table ant[NANT];
  source_table source;
  array_data array;
  time_t thetime;
  struct tm *gm_time;
  double jdtime;

  /* test calendar date to Julian date and back */
  Cal_to_JD(2007,02,18,&jdtime);
  JD_to_Cal(jdtime,&year,&mon,&day);
  printf("calcheck: 2007/2/18 in jd is: %g. Back again: %d/%d/%d \n",jdtime,year,mon,day);

  time(&thetime);
  gm_time = gmtime(&thetime);
  Cal_to_JD(gm_time->tm_year+1900,gm_time->tm_mon+1,gm_time->tm_mday,&jdtime);
  printf("now is: %d/%d/%d JD: %f\n",gm_time->tm_year+1900,gm_time->tm_mon+1,gm_time->tm_mday,jdtime);
  JD_to_Cal(jdtime,&year,&mon,&day);
  printf("calcheck: now is: %d/%d/%d \n",year,mon,day);

  printf("size of vis dump: %d\n",siz_visdump);

  /* make data space */
  data.n_baselines = calloc(NVIS,sizeof(int));
  data.baseline = calloc(NVIS,sizeof(float));
  for (i=0; i<NVIS; i++){
    vis[i] = calloc(siz_visdump,sizeof(float)*2); /*nbase*npol*nfreq complex visibilities per time dump */
    wt[i] = calloc(siz_visdump,sizeof(float));    /*nbase*npol*nfreq weights per time dump */
    u[i] = calloc(NBASE,sizeof(double));          /* ther are NBASE baselines per time time */
    v[i] = calloc(NBASE,sizeof(double));
    w[i] = calloc(NBASE,sizeof(double));
    data.n_baselines[i] = NBASE;
    data.baseline[i] = calloc(NBASE,sizeof(float));
  }

  /* fill in some values for our fake data. Basically populate all the data of the uvdata structure. */
  array.n_ant=NANT;
  data.n_pol=NPOL;
  data.n_vis=NVIS;
  data.n_freq = NFREQ;
  data.cent_freq = 140.0e6;
  data.freq_delta = 32e3;
  data.antennas = ant;
  data.source = &(source);
  data.u = u;
  data.v = v;
  data.w = w;
  data.date =t;
  data.array = &array;
  array.xyz_pos[0] = -2613097.0;
  array.xyz_pos[1] = 5083361.0;
  array.xyz_pos[2] = -2822005.0;
  data.visdata = vis;
  data.weightdata = wt;
  sprintf(source.name,"TEST SOURCE");
  source.ra=22.6;
  source.dec=3.08;
  sprintf(array.name,"MWA-SIM");

  /* create antenna indexes for baselines. starts at 1, not 0 */
  for(l=0; l<NVIS; l++) {
    i=0;
    for (j=1; j<NANT; j++) {
      for (k=j+1; k<=NANT; k++) {
        EncodeBaseline(j,k,data.baseline[l]+i);
	    i++;
      }
    }
  }

  fprintf(stderr,"filling antenna info.\n");
  for (i=0; i<NANT; i++) {
    sprintf(ant[i].name,"ANT%d",i+1);
    ant[i].xyz_pos[0] = 0.0;  /* E-W array has only Y offsets */
    ant[i].xyz_pos[1] = ant_y_offsets[i]*15;
    ant[i].xyz_pos[2] = 0.0;
    ant[i].xyz_deriv[0] = 0;
    ant[i].xyz_deriv[1] = 0;
    ant[i].xyz_deriv[2] = 0;
    sprintf(ant[i].pol_typeA,"X");
    sprintf(ant[i].pol_typeB,"Y");
    ant[i].pol_angleA = 0.0;
    ant[i].pol_angleB = 90.0;
    ant[i].pol_calA = 0.0;
    ant[i].pol_calB = 0.0;
    ant[i].mount_type = 0;
  }

  fprintf(stderr,"filling visibilities.\n");
  for (i=0; i<NVIS; i++) {
    /* temp pointers for this time dump */
    u_ptr = u[i];
    v_ptr = v[i];
    w_ptr = w[i];
    vis_ptr = vis[i];
    wt_ptr = wt[i];
    t[i] = TIMESTEP/86400.0*i+jdtime+gm_time->tm_hour/24.0+gm_time->tm_min/(24.0*60.0)+gm_time->tm_sec/(86400.0);
    for (l=0; l<data.n_baselines[i]; l++){
      DecodeBaseline(data.baseline[i][l],&ant1, &ant2);
      calcUVW(TIMESTEP/86400.0*i*2.0*M_PI,-90.0,
	      ant[ant1-1].xyz_pos[0]-ant[ant2-1].xyz_pos[0],
	      ant[ant1-1].xyz_pos[1]-ant[ant2-1].xyz_pos[1],
	      ant[ant1-1].xyz_pos[2]-ant[ant2-1].xyz_pos[2],
	      u_ptr+l,v_ptr+l,w_ptr+l);
      u_ptr[l] /= 3e8;
      v_ptr[l] /= 3e8;
      w_ptr[l] /= 3e8;
      /*      printf("%g,%g,%g\n",u_ptr[l],v_ptr[l],w_ptr[l]); */
      for (j=0; j<data.n_freq; j++) {
        for(k=0; k<data.n_pol; k++) {
        vis_ptr[(l*data.n_freq*data.n_pol+j*NPOL+k)*2  ] = 0.70+l*1e-3+j*1e-5;
        vis_ptr[(l*data.n_freq*data.n_pol+j*NPOL+k)*2+1] = -0.70-l*1e-3-j*1e-5;
        wt_ptr[l*data.n_freq*data.n_pol+j*NPOL+k] = 100.0;
        }
      }
    }
  }

  fprintf(stderr,"writing fits file.\n");
  writeUVFITS("testfile.fits",&data);
/*
  printUVData(&data,stdout);
*/

  return 0;
}

void calcUVW(float HA, float DEC, float X, float Y, float Z, double *u, double *v, double *w) {
  float delta;

  delta = DEC*M_PI/180.0;

  *u =  sin(HA)*X + cos(HA)*Y;
  *v = -sin(delta)*cos(HA)*X + sin(delta)*sin(HA)*Y + cos(delta)*Z;
  *w =  cos(delta)*cos(HA)*X - cos(delta)*sin(HA)*Y + sin(delta)*Z;
}
