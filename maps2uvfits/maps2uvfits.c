/* Program to convert the binary output format of MAPS visgen into a
   UV FITS format. Written by Randall Wayth. Sep, 2006.
 */

/* a note about sign conventions etc:
  The de-facto standard for UVFITS files follows the AIPS convention:
  - a baseline is defined as location(ant1)-location(ant2) where the antenna indices are such that ant2 > ant1
  - using this definition of a baseline, a visibility (for a point source) is V(u,v) = I*exp(-2*pi*(ul+vm))
  - this means that if a baseline is defined as (ant2-ant1),
    then the exponent in the exponential must be +2pi, not -2pi
  - The AIPS standard is consistent with the definition of a visibility from a coherence average point of view
    (e.g section 14.1 of Thompson, Moran and Swenson or Chapter 1 of "Synthesis imaging in radio astronomy II".)
    

  - In MAPS, the baseline is defined as (ant1-ant2), and visibilities have a -2pi exponent.
  - Miriad stores visibilities internally in a way which is consistent with this standard, but which has the
    opposite definition of baseline (i.e ant2-ant1), so when UVFITS files are read into miriad, the baseline
    is reversed (u,v,w)<-(-u,-v,-w) AND the vis is conjugated.
    
    There are also various memos/articles about this lying in the deep bowels of the internet:
    http://www.atnf.csiro.au/observers/memos/d97383~1.pdf
    http://adsabs.harvard.edu/abs/1981A%26AS...44..371G
*/


#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <ctype.h>
#include "visibility.h"
#include "comp_func.h"
#include "uvfits.h"

#define MAX_ANT 10000                // INCREASED 2010 Apr 16
#define MAX_LINE 1000                // INCREASED 2010 Apr 16
#define MWA_LAT -26.43	// Array latitude. degrees North
#define MWA_LON 0.0   // Array longitude. degrees East
#define MWA_HGT 400	// Array altitude. meters above sea level
#define EARTH_RAD 6378100.0  // meters

/* private function prototypes */
void printheader(FILE *fp, struct main_header *head);
void printusage(char *progname);
int readScan(FILE *fp, uvdata *data);
int createPolIndex(char *polprods, int *index);
int readArray(char *filename, uvdata *data, double lat);

int debug=0;
int pol_index[4];
FILE *fpd=NULL;


/************************
************************/
int main(int argc, char *argv[]) {
  FILE *fpin;
  char *stationfilename=NULL;
  int n_read,scan=0,res=0,nstokes=0,i;
  struct main_header main_header;
  uvdata data;
  array_data *arraydat=NULL;
  source_table *source=NULL;
  ant_table *antennas=NULL;
  double jdtime_base,lat=MWA_LAT*(M_PI/180.0),lon=MWA_LON*(M_PI/180.0),height=MWA_HGT;

  fpd=stdout;

  /* print a 'usage' message if necessary */
  if (argc < 6) printusage(argv[0]);
  if (strcmp(argv[argc-1],"-d")==0) debug=1;

  if (argc == 7 && !debug){
    /* check for a station file filename */
    stationfilename=argv[6];
  }
  if (argc == 8 && debug){
    /* check for a station file filename */
    stationfilename=argv[6];
  }

  lat = atof(argv[3]);
  lon = atof(argv[4]);
  height = atof(argv[5]);
  if (fabs(lat) > 90 || fabs(lon) > 180 || fabs(height) > 1e4) {
    fprintf(stderr,"Illegal lat or lon or height\n");
    printusage(argv[0]);
  }
  lat *= (M_PI/180.0);
  lon *= (M_PI/180.0);

  /* initialise some values for the UV data array. Because the MAPS output
   * does not include how many baselines there are in the header, we can't know in
   * advance how much mem to allocate. We will make some here and use realloc()
   * where necessary */
  antennas = calloc(MAX_ANT,sizeof(ant_table));
  arraydat = calloc(1,sizeof(array_data));
  source = calloc(1,sizeof(source_table));
  if(antennas==NULL || arraydat==NULL || source ==NULL) {
    fprintf(stderr,"no malloc\n");
    exit(1);
  }
  data.date = calloc(1,sizeof(double));
  data.array = arraydat;
  data.source= source;
  data.antennas=antennas;
  data.n_pol=0;
  data.n_baselines=NULL;
  data.n_freq=0;
  data.n_vis=0;
  data.cent_freq=0.0;
  data.freq_delta = 0.0;
  data.u=NULL;
  data.v=NULL;
  data.w=NULL;
  data.baseline=NULL;
  data.weightdata=NULL;
  data.visdata=NULL;
  data.pol_type=1;    /* default is Stokes pol products */
  data.fq = NULL;
  arraydat->n_ant=0;
  strcpy(arraydat->name,"MWA-SIM");

  if ((fpin=fopen(argv[1],"r"))==NULL) {
    fprintf(stderr,"cannot open input file %s\n",argv[1]);
    exit(1);
  }

  /* for debugging: print out the sizes of various structs */
  if (debug) {
    fprintf(fpd,"sizeof struct visfreq: %d\n",(int)sizeof(struct visfreq));
    fprintf(fpd,"sizeof struct visgroup: %d\n",(int)sizeof(struct visgroup));
    fprintf(fpd,"sizeof complex: %d\n",(int)sizeof(complex));
  }

  /* assign XYZ positions of the array for the site. */
  /* FIXME this is not accurate. UVFITS writer now contains correct geodetic
     conversion. Need to update MAPS to use this also before fixing here.*/
  arraydat->xyz_pos[0] = (EARTH_RAD+height)*cos(lat)*cos(lon);
  arraydat->xyz_pos[1] = (EARTH_RAD+height)*cos(lat)*sin(lon);
  arraydat->xyz_pos[2] = (EARTH_RAD+height)*sin(lat);
  /* assign XYZ positions for antennas */
  if (stationfilename!=NULL) {
    readArray(stationfilename, &data,lat);
  }
  else {
    fprintf(stderr,"ERROR: no array file.\n");
    exit(1);
  }

  /* populate antenna info */
  if (debug) fprintf(fpd,"there are %d antennas\n",arraydat->n_ant);
  for (i=0; i<arraydat->n_ant; i++){
    sprintf(antennas[i].name,"ANT%03d",i+1);
    sprintf(antennas[i].pol_typeA,"X");
    sprintf(antennas[i].pol_typeB,"Y");
    antennas[i].pol_angleA = 0.0;
    antennas[i].pol_angleB = 90.0;
    antennas[i].pol_calA = 0.0;
    antennas[i].pol_calB = 0.0;
    antennas[i].mount_type = 0;
  }

  /* read the main header */
  n_read = fread(&main_header,sizeof(struct main_header),1,fpin);
  if (n_read != 1) {
    fprintf(stderr,"Failed to read main header. exiting.\n");
    exit(1);
  }
  if (debug) printheader(fpd,&main_header);

  /* how many pol products... */
  createPolIndex(main_header.stokes, pol_index);
  for(i=0; i<4; i++) if (main_header.stokes[2*i] !='\0') nstokes++;
  if(nstokes<1 || nstokes > 4) {
    fprintf(stderr,"bad number of stokes: %d\n",nstokes);
    exit(1);
  }
  /* set the polarisation product type. linear, circular or stokes */
  if (toupper(main_header.stokes[0]) == 'X' || toupper(main_header.stokes[0]) == 'Y') data.pol_type=-5;
  if (toupper(main_header.stokes[0]) == 'R' || toupper(main_header.stokes[0]) == 'L') data.pol_type=-1;
  if (debug) fprintf(fpd,"Found %d pol products. pol_type is: %d\n",nstokes,data.pol_type);

  /* fill in some data */
  /* everyone has to use their own special date format... */
  Cal_to_JD(main_header.reftime.year,1,1,&jdtime_base);
  jdtime_base += main_header.reftime.day + main_header.reftime.hour/24.0 + main_header.reftime.minute/1440.0 +
                 main_header.reftime.second/86400.0;
  data.date[0] = jdtime_base;
  if (debug) fprintf(fpd,"JD time is %.2lf\n",jdtime_base);

  data.n_pol = nstokes;
  strncpy(source->name,main_header.field_name,SIZE_SOURCE_NAME-1);

  /* extract RA, DEC from header. Note negative zero bug. */
  source->ra = (main_header.field_center.ra_hrs + main_header.field_center.ra_mins/60.0 +
                main_header.field_center.ra_secs/3600.0);

  source->dec = fabs(main_header.field_center.dec_degs) + fabs(main_header.field_center.dec_mins)/60.0 +
                fabs(main_header.field_center.dec_secs/3600.0);
  if (main_header.field_center.dec_degs<0||main_header.field_center.dec_mins<0||main_header.field_center.dec_secs<0) {
    source->dec = - source->dec;
  }

  //  printf("ra: %g, dec: %g\n",source->ra,source->dec);
  /* for each scan */
  for (scan=0; scan < main_header.ntime_block; scan++) {
    if(debug) fprintf(fpd,"Reading scan %d of %d\n",scan+1,main_header.ntime_block);
    res = readScan(fpin,&data);
    if(res!=0) {
      fprintf(stderr,"Problems in readScan(). exiting\n");
      exit(1);
    }
  }

  if (debug){
    fprintf(fpd,"writing the UVFITS file\n");
    fprintf(fpd,"there are %d time sets of visibilities with %d baselines, %d IFs, %d channels\n",
			data.n_vis,data.n_baselines[0], data.n_IF, data.n_freq);
  }

  /* we now have all the data, write it out */
  writeUVFITS(argv[2],&data);

  // Testing: write it out to a test file
  //FILE *fp=fopen("original.uvfits", "w");
  //printUVData(&data, fp);
  
  if(debug) fprintf(fpd,"finished writing UVFITS file\n");
  freeUVFITSdata(&data);
  return 0;
}


/***************************
 ***************************/
int readScan(FILE *fp, uvdata *uvdata) {
  static int found_timeblock=0;

  char temp[20];
  int n_read,i,j,visindex=0,int_no=0,f;
  int n_index; // for storing array index from baseline, freq info
  int more_visblocks=1,have_visfreq=1;
  int baseline=0,n_integrations=0,max_baselines=0;
  struct freq_span *freqdata=NULL;
  struct time_block_header timehead;
  struct visblock visblock;
  struct visfreq *visfreq;
  struct visgroup *visgroup;
  float cent_freq=0.0, bw=0.0;

  visfreq = visblock.visfreq;
  visgroup = visfreq->visgroup;
  visindex = uvdata->n_vis;
  max_baselines = (uvdata->array->n_ant*(uvdata->array->n_ant + 1)/2); // can include autocorrelations...
  /* read the timeblock and header */
  if (found_timeblock) {
    if(debug) fprintf(fpd,"found a timeblock header previously\n");
  }
  else{
    n_read = fread(temp,sizeof(char),8,fp);
    if (n_read != 8 || strcmp(temp,"timblck")) {
      fprintf(stderr,"expected 'timblck', but it was not found. found: %s\n",temp);
      return 1;
    }
  }

  /* time block and vis block structure definitions have clunky 1-element array
   * in them which is extended in the writer to an arbitrary length array. So
   * don't read that section here, leave it for the real array once we know the
   * length */
  fread(&timehead,sizeof(struct time_block_header)-sizeof(struct freq_span),1,fp);
  n_integrations = rint(timehead.duration/timehead.intg_time);
  if (debug) {
    fprintf(fpd,"There are %d IFs\n",timehead.num_IFs);
    fprintf(fpd,"start time: %g, duration: %g, integration time: %g, total integrations: %d\n",
            timehead.start_time, timehead.duration, timehead.intg_time,n_integrations);
  }
  //if(timehead.num_IFs > 1) {
  //  fprintf(stderr,"ERROR: this code does not yet handle multiple IFs in the one file.\n");
  //  return 1;
  //}

  /* we know the number of IFs now, so make space for that array */
  freqdata = calloc(timehead.num_IFs,sizeof(struct freq_span));
  if(freqdata ==NULL) {
    fprintf(stderr,"readScan: no malloc for freq_span\n");
    exit(1);
  }

  /* read in the IF freq and bandwidth info */
  n_read = fread(freqdata,sizeof(struct freq_span),timehead.num_IFs,fp);
  if (n_read !=timehead.num_IFs) {
    fprintf(stderr,"readScan: expected to read %d freq_span records, only read %d\n",timehead.num_IFs,n_read);
    return 1;
  }
  if (debug) for(i=0; i<timehead.num_IFs;i++) fprintf(fpd,"IF %d: Freq: %g, BW: %g. (MHz)\n",i,freqdata[i].frequency,freqdata[i].bandwidth);
	
  /* read the visblock header for the number of freq channels in the IF */
  n_read = fread(temp,sizeof(char),8,fp);
  if (n_read != 8 || strcmp(temp,"visblck")) {
    fprintf(stderr,"expected 'visblck', but it was not found. read %d bytes. found: <%s>. EOF: %d\n",n_read,temp,feof(fp));
    return 1;
  }
  fread(&visblock,sizeof(struct visblock)-sizeof(struct visfreq),1,fp);
  
  if (debug) {
    fprintf(fpd,"There are %d IFs in this visblock\n",visblock.nfreq);
  }
  if (visblock.nfreq != timehead.num_IFs) {
    fprintf(stderr,"ERROR: num IFs in timeblock header (%d) is different to num in visblock header (%d)\n",
            timehead.num_IFs,visblock.nfreq);
    return 1;            
  }

  /* now read the visfreq header since that tells us how many actual channels there are */
  n_read = fread(visfreq,1,sizeof(struct visfreq)-sizeof(struct visgroup),fp);
  if (n_read != (sizeof(struct visfreq)-sizeof(struct visgroup))) {
    fprintf(stderr,"failed reading read %d bytes for visfreq header. Only read %d bytes.\n",
                (int)(sizeof(struct visfreq)-sizeof(struct visgroup)),n_read);
    return 1;
  }
  if(visfreq->chan != 0) {
      fprintf(stderr,"readScan: read index for freq channel %d, expected 0\n",visfreq->chan);
      return 1;
  }
  if (debug) fprintf(fpd,"read vis_freq. IF: %d of %d, n corr chans: %d\n",visfreq->chan+1,visblock.nfreq, visfreq->ncchan);
  uvdata->n_freq = visfreq->ncchan;
    
  // Set the IF in uvdata to match the number of IFs in the data, and allocate memory
  uvdata->n_IF = timehead.num_IFs;
  uvdata->fq = calloc(uvdata->n_IF, sizeof(fq_table));
  
  // Fill the data based on the freqdata list
  for(i=0; i<uvdata->n_IF; i++) {
    uvdata->fq[i].freq = freqdata[i].frequency*1e6;  // frequency of this IF
    uvdata->fq[i].bandwidth = freqdata[i].bandwidth*1e6; // bandwidth of this IF
    uvdata->fq[i].chbw = freqdata[i].bandwidth*1e6/uvdata->n_freq; // width of each band (based on number of bands)
    uvdata->fq[i].sideband = 1; // Not sure which to set this as??
  }

  /* make space for the actual visibilities and weights */
  uvdata->n_vis += n_integrations;
  uvdata->date=realloc(uvdata->date,uvdata->n_vis*sizeof(double));
  uvdata->visdata=realloc(uvdata->visdata,uvdata->n_vis*sizeof(double *));
  uvdata->weightdata=realloc(uvdata->weightdata,uvdata->n_vis*sizeof(double *));
  uvdata->u=realloc(uvdata->u,uvdata->n_vis*sizeof(double *));
  uvdata->v=realloc(uvdata->v,uvdata->n_vis*sizeof(double *));
  uvdata->w=realloc(uvdata->w,uvdata->n_vis*sizeof(double *));
  uvdata->n_baselines = realloc(uvdata->n_baselines,uvdata->n_vis*sizeof(int));
  uvdata->baseline = realloc(uvdata->baseline,uvdata->n_vis*sizeof(float *));
	
  // Allocate arrays of data for each time step (baselines * n_IF * n_freq * n_pol)
  // (new dimension is for each IF, was previously only 1 IF)
  if (debug) fprintf(fpd,"callocing %d arrays of %d floats\n",n_integrations,max_baselines*uvdata->n_IF*uvdata->n_freq*uvdata->n_pol*2);
  for (i=0; i< n_integrations; i++) {
    uvdata->n_baselines[visindex+i] = 0;
    uvdata->baseline[visindex+i] = NULL;
    uvdata->visdata[visindex+i] = calloc(max_baselines*uvdata->n_IF*uvdata->n_freq*uvdata->n_pol*2,sizeof(float));
    uvdata->weightdata[visindex+i] = calloc(max_baselines*uvdata->n_IF*uvdata->n_freq*uvdata->n_pol  ,sizeof(float));
    uvdata->u[visindex+i] = calloc(max_baselines,sizeof(double));
    uvdata->v[visindex+i] = calloc(max_baselines,sizeof(double));
    uvdata->w[visindex+i] = calloc(max_baselines,sizeof(double));
	  
	  
	  
	  
    if(uvdata->visdata[visindex+i]==NULL || uvdata->weightdata[visindex+i]==NULL || uvdata->visdata[visindex+i]==NULL
        || uvdata->visdata[visindex+i]==NULL || uvdata->visdata[visindex+i]==NULL) {
        fprintf(stderr,"readScan: no malloc for BIG arrays. Try reducing MAX_ANT.\n");
        exit(1);
    }
    
    /* set the date */
    /* MAPS counts time from the start of a window, but it should
     be the middle of a window, so correct for this time offset */
    if (visindex+i==0) {
        //this isn't needed any more (time done properly in visgen). Not sure about below the else.
        //uvdata->date[0] += 0.5*timehead.intg_time/86400.0;
    } else {
		// Removed 0.5 from the line below to get time offsets right (relative to first and that's offset already)
        uvdata->date[visindex+i] = uvdata->date[0]+(timehead.start_time+i*timehead.intg_time)/86400.0;
    }
    if (debug) {
        fprintf(fpd,"Vis index %d, date offset is: %g\n",visindex+i,uvdata->date[visindex+i]-uvdata->date[0]);
    }

  }

  /* in MAPS, the freq is the bottom of the frequency channel, but it
     is supposed to be the middle of the channel, so correct for this */
  // check that the freq hasn't changed
  bw = freqdata[0].bandwidth/uvdata->n_freq;
  cent_freq = (freqdata[0].frequency+bw/2)*1e6;
  //cent_freq = (freqdata[0].frequency+freqdata[0].bandwidth/2)*1e6;
  if (uvdata->cent_freq==0.0) {
    uvdata->cent_freq  = cent_freq;  /* freq and bandwidth in Hz */
    uvdata->freq_delta = bw*1e6;
  }
  else {
    if (uvdata->cent_freq != cent_freq) {	
        fprintf(stderr,"ERROR: freq changed. Not supported. Old %g, new: %g\n",uvdata->cent_freq,cent_freq);
        return -1;
    }
  }

  /* sub-scan integrations are written to this file with index "intg_no" since the MAPS MPI is split in this way.
     when there is more than one sub-integration, we want to keep them separate and write them as if they were
     individual correlator dumps for a time period. */

  /* there is now one visblock for each baseline (unknown number). Keep looping
     until we hit something that isn't a visblock header */
  while (more_visblocks) {

    if(visblock.intg_no != int_no) {
        int_no = visblock.intg_no;
        baseline=0;
        if (debug) fprintf(fpd,"found new integration block %d\n",visblock.intg_no);
    }

    if (debug) {
      fprintf(fpd,"integ num: %d, baseline: %d-%d (%d)\n",visblock.intg_no,visblock.station1,visblock.station2,baseline);
      fprintf(fpd,"uvw: %g, %g, %g\n",(float)visblock.u,(float)visblock.v,(float)visblock.w);
      fprintf(fpd,"num freqs %d\n",visblock.nfreq);
    }
    
    /* populate the baseline info. Antenna numbers start at 0 in MAPS, 1 in UVFITS.  */
    uvdata->baseline[visindex+int_no] = realloc(uvdata->baseline[visindex+int_no],(baseline+1)*sizeof(float));
    uvdata->baseline[visindex+int_no][baseline] = 0;
    EncodeBaseline(visblock.station1+1, visblock.station2+1, uvdata->baseline[visindex+int_no]+baseline);

    uvdata->u[visindex+int_no][baseline] = visblock.u/VLIGHT;
    uvdata->v[visindex+int_no][baseline] = visblock.v/VLIGHT;
    uvdata->w[visindex+int_no][baseline] = visblock.w/VLIGHT;

	  
    // Need to read all of the IFs (1 visfreq for each IF)
	for(f=0; f<uvdata->n_IF; f++) {
		/* get the visfreq struct. The first time around, we have already read it above */
		if(!have_visfreq) {
			n_read = fread(visfreq,1,sizeof(struct visfreq)-sizeof(struct visgroup),fp);
			
			// Will not always be zero if there are multiple IFs
			//if(visfreq->chan != 0) {
			//	fprintf(stderr,"readScan: read index for freq channel %d, expected 0\n",visfreq->chan);
			//	return 1;
			//}
			
			if (debug) fprintf(fpd,"read visfreq (%d bytes). IF: %d of %d, n corr chans: %d\n",n_read, visfreq->chan+1,visblock.nfreq, visfreq->ncchan);
		}
		else {
			have_visfreq=0;
		}
				  
		// sanity check to make sure the number of correlator channels hasn't suddenly changed
		if(visfreq->ncchan != uvdata->n_freq) {
			fprintf(stderr,"ERROR: number of correlator channels has changed from %d to %d\n",uvdata->n_freq,visfreq->ncchan);
			return 1;
		}


		// now read all the visgroups, one for each correlator channel
		for (i=0; i < visfreq->ncchan; i++) {
		  
		  //      n_read = fread(visgroup,1,sizeof(struct visgroup)+(uvdata->n_pol-1)*sizeof(complex),fp);
		  n_read = fread(visgroup,1,sizeof(struct visgroup),fp);
		  if (debug) {
			fprintf(fpd,"read visgroup (%d bytes). cchan: %d, flag: %d, weight: %g\n",n_read,
					visgroup->cchan,visgroup->flag,(float)visgroup->weight);
			for (j=0; j<uvdata->n_pol; j++) {
			  fprintf(fpd,"pol: %d, re: %g, im: %g\n",j,c_real(visgroup->vis[j]),c_imag(visgroup->vis[j]) );
			}
		  }

		  /* populate the visibility arrays */
		  /* NOTE: we conjugate the visibility here to be consistent with the format of UVFITS files */
		  for (j=0; j<uvdata->n_pol;j++) {
        n_index = baseline*uvdata->n_IF*uvdata->n_freq*uvdata->n_pol + f*uvdata->n_freq*uvdata->n_pol + i*uvdata->n_pol + j;
        
        uvdata->visdata[visindex+int_no][n_index*2   ] = c_real(visgroup->vis[pol_index[j]]);
			  uvdata->visdata[visindex+int_no][n_index*2 +1] = -c_imag(visgroup->vis[pol_index[j]]);
			  uvdata->weightdata[visindex+int_no][n_index] = visgroup->weight;
        
			  //uvdata->visdata[visindex+int_no][(baseline*uvdata->n_pol*visfreq->ncchan + i*uvdata->n_pol+j)*2   ] = c_real(visgroup->vis[pol_index[j]]);
			  //uvdata->visdata[visindex+int_no][(baseline*uvdata->n_pol*visfreq->ncchan + i*uvdata->n_pol+j)*2 +1] = -c_imag(visgroup->vis[pol_index[j]]);
			  //uvdata->weightdata[visindex+int_no][(baseline*uvdata->n_pol*visfreq->ncchan + i*uvdata->n_pol)+j] = visgroup->weight;
		  } // n_pol
		} // ncchan
	} // n_IF

    /* look for another visblock header */
    temp[0] = '\0';
    found_timeblock=0;
    n_read = fread(temp,sizeof(char),8,fp);
    if (n_read != 8 || strcmp(temp,"visblck")) {
      if (debug) fprintf(fpd,"read %d bytes. Done with visblcks. Got '%s' instead\n",n_read,temp);
      if (strcmp(temp,"timblck")==0) found_timeblock=1;
      more_visblocks=0;
    } else{
      n_read = fread(&visblock,sizeof(struct visblock)-sizeof(struct visfreq),1,fp);
      if(debug) fprintf(fpd,"Found next visblck. Read %d structure\n",n_read);
    }
    /* count baselines for each scan (time interval)  */
    (uvdata->n_baselines[visindex+int_no])++;
    baseline++;
  }

  if(freqdata !=NULL) free(freqdata);
  return 0;
}


/***************************
 ***************************/
void printusage(char *progname) {
  fprintf(stderr,"Usage:\n");
  fprintf(stderr,"%s <mapsfile> <uvfile> <latitude degrees> <longitude degrees East> <altitude meters> [stationfile] [-d (debug)]\n",progname);
  exit(1);
}


/****************************
*****************************/
void printheader(FILE *fp, struct main_header *head) {
  char temp[12];

  fprintf(fp,"name:\t%s\n",head->simname);
  fprintf(fp,"date:\t%d-%03d %d:%d:%f\n",(int)head->reftime.year,(int)head->reftime.day,
	  (int)head->reftime.hour,(int)head->reftime.minute,(float)head->reftime.second);
  fprintf(fp,"Field:\t%d:%d:%g (hours) %d:%d:%g (deg)\n",
          (int)head->field_center.ra_hrs,(int)head->field_center.ra_mins,head->field_center.ra_secs,
          (int)head->field_center.dec_degs,(int)head->field_center.dec_mins,head->field_center.dec_secs);
  fprintf(fp,"Field:\t%s\n",head->field_name);
  memset(temp,'\0',12);
  strncpy(temp,head->stokes,SIZ_PRODNAME);
  fprintf(fp,"Pol prods:\t<%s>\n",temp );
  fprintf(fp,"nscan:\t%d\n",head->ntime_block);
}


/***************************************
 examine the ordering of polarisation products from maps and create an
 index to order them the way miriad likes: i.e. XX, YY, XY, YX
 typically they will be XX,XY,YX,YY from MAPS
 **************************************/
int createPolIndex(char *polprods, int *index) {
  int p1='\0',p2='\0',i;

  /* find the unique letters representing the pol products. i.e. X,Y or R,L or II */
  for (i=0; i<SIZ_PRODNAME; i++) {
    if (p1=='\0' && polprods[i] != '\0') {
      p1 = polprods[i];
      continue;
    }
    if (p2=='\0' && polprods[i] != '\0' && p1 != polprods[i]) {
      p2 = polprods[i];
      continue;
    }
  }
  if (debug) fprintf(fpd,"Found pol keys '%c' and '%c'\n",p1,p2);

  /* find the index of products */
  for (i=0; i<4; i++) {
    if (polprods[i*2]==p1 && polprods[i*2+1]==p1) index[0] = i;
    if (polprods[i*2]==p2 && polprods[i*2+1]==p2) index[1] = i;
    if (polprods[i*2]==p1 && polprods[i*2+1]==p2) index[2] = i;
    if (polprods[i*2]==p2 && polprods[i*2+1]==p1) index[3] = i;
  }
  if (debug) {
    for (i=0; i<4; i++) fprintf(fpd,"polindex: %d ",index[i]);
    fprintf(fpd,"\n");
  }
  return 0;
}


/******************************
 read the station locations from a MAPS array
 file and populate the antenna positions in the
 data structure.
*******************************/
int readArray(char *filename, uvdata *data, double lat) {
  FILE *fp=NULL;
  char line[MAX_LINE];
  char items[8][MAX_LINE];
  int index=0,nscan=0,vlbi_coords=0;
  float east=0,north=0,height=0;

  if( (fp=fopen(filename,"r"))==NULL) {
    fprintf(stderr,"ERROR: readArray: failed to open array file <%s>\n",filename);
    return 1;
  }

  /* scan through lines. convert MAPS east,north,height units to XYZ units */
  /* array specs with more than 3 items in a line are special VLBI units and need
     no conversion */
  while((fgets(line,MAX_LINE-1,fp)) !=NULL) {
    if(line[0]=='\n' || line[0]=='#' || line[0]=='\0') continue; // skip blank/comment lines
    nscan = sscanf(line,"%s %s %s %s %s %s %s %s",items[0],items[1],items[2],items[3],items[4],items[5],items[6],items[7]);
    if (nscan < 5) {
        /* old style coords in East, North, Height */
        east = atof(items[0]);
        north = atof(items[1]);
        ENH2XYZ_local(east,north,height,lat, data->antennas[index].xyz_pos, data->antennas[index].xyz_pos+1,
	            data->antennas[index].xyz_pos+2);
	}
	else {
	    /* new style VLBI: name,X,Y,Z,station,el_low,el_high,sefd */
	    data->antennas[index].xyz_pos[0] = atof(items[1]);  //X
	    data->antennas[index].xyz_pos[1] = atof(items[2]);  //Y
	    data->antennas[index].xyz_pos[2] = atof(items[3]);  //Z
	    strncpy(data->antennas[index].name,items[0],SIZE_ANT_NAME); // copy the name across also
	    vlbi_coords=1;
	}
    if (debug) {
      fprintf(fpd,"ant %s. Pos (ENH) (%g,%g,%g).\tPos (XYZ): (%g,%g,%g)\n",data->antennas[index].name,
	      east,north,height,
	      data->antennas[index].xyz_pos[0],data->antennas[index].xyz_pos[1],data->antennas[index].xyz_pos[2]);
    }
    index++;
  }
  fclose(fp);
  
  /* if we're using VLBI coords, set the array center to center of Earth */
  if (vlbi_coords) {
    data->array->xyz_pos[0] = 0.0;
    data->array->xyz_pos[1] = 0.0;
    data->array->xyz_pos[2] = 0.0;
  }
  if (debug) {
    fprintf(fpd,"Array coords XYZ: (%g,%g,%g), lat: %g\n",
	    data->array->xyz_pos[0],data->array->xyz_pos[1],data->array->xyz_pos[2],lat);
  }
  data->array->n_ant = index;
  return 0;
}


