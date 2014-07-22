/*
 * check_baseline_unsorted_bblt.c
 * In this version of check_baseline.c the oversized baselines are added 
 * to struct bl_trespass line by line without any sorting.
 * I am working now on a new version of check_baseline.c, where 
 * struct bl_trespass will be filled in using insertion between the 
 * lines so that after the run the table is sorted in descending order
 * if u and v. The insertion will not stop when the table is full;
 * The lowest lines will just vanish after the insertion.
 */
/************************************************************************/
/* check_baseline()                                                    */
/* This program is a "mini visgen" in the sence it generally has the    */
/* same structure and logics, but it does not fulfil the integration to */
/* generate the visibilities, and it writes no files.                   */
/* It finds the conditions (maximum frequency etc.) at which the        */
/* baseline can reach maximum lengths, loops over time and baseline    */
/* and checks if the uv-patches are contained within the specified      */
/* uv-plane dimensions. Successful run of check_baseline() guarantees  */
/* absense of error conditions due to forming patches beyond the uv-    */
/* plane boundary and due to the patch size exceeding the maximum       */
/* during the subsequent visgen run. */
/*                                                                      */
/* in order to generate individual simulated visibilities, which are    */
/* written to an output file, ready for reading by AIPS++.  The program */
/* requires as input a description of the simulated array, a uv grid    */
/* representing the fourier transformed simulated sky, a description of */
/* the simulated observing parameters, and optionally, files containing */
/* a model ionosphere, and a model RFI environment.                     */
/*                                                                      */
/*      Inputs:     Files and descriptions, as above                    */
/*                                                                      */
/*      Output:     Simulated visibility dataset                        */
/*                                                                      */
/* Created Nov 19th 2010 by L. benkevitch                               */
/*                                                                      */
/************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <math.h>
#include <string.h>
#include <errno.h>
#include <ctype.h>
#include <sys/time.h>
#include <time.h>
#include <float.h>

/* header files of the original visgen */
#include "check_baseline.h"
#include "visgen.h"
#include "visibility.h"
#include "observing_param.h"
#include "array_specification.h"
#include "station_beam.h"
#include "rfi_model.h"
#include "convolve.h"
#include "station_beam.h"
#include "uvgrid.h"
#include "account.h"
#include "comp_func.h"
#include "type_radec.h"
#include "sizes.h"
#include "add_noise.h"
#include "makebeam.h"
#include "oob.h"
#include "ionosphere.h"
#include "pol_response.h"
#include "get_patch.h"
#include "coord_trans.h"
#include "compute_beamconvl.h"
#include "visgen_prototypes.h"
/* global variables */
extern int do_acc;

#define RADASEC 206264.80624709635515647 /* number of arcsec in a radian  */

//extern int msglev
//extern char timblock;
//extern char vblock;
//extern char progname[]; 

int compute_ionoscreen_oob(const struct stationspec *stn,
			   const int                 st_index,
			   const double              gha,
			   const double              fov_ra,
			   const int                 n_oob,
			   const source_model_t     *ooblist,
			   struct ion_model         *isphere,
			   const int                 center_iono,
			   oob_beam_t               *oob_beams);


int check_baseline(struct vg_param vgparam)
{
  off_t refpoint, write_pos;
  
  /* part of the original visgen */
  int  i,scan=0, intg=0, st1, st2, ch, ncch=0, cch;
  int  bufsize=0, noob=0, nelev=0, ret=0; 
  /* int timeflag; */ 
  int  intg_start=0, intg_num=0, elev_flag;
  double start, time, bw, cbw, obschan_freq, freq;
  double u, v,w, u_len, v_len, w_len;
  double scanstart, gha, fov;  // gha in radian. 
  double udim, vdim,reltime=0, time_ut=0.0,time_scan0=0;
  double ainit, na;
  /* double umean, vmean;    */
  /* FILE *outfp;            */
  struct observing_param obsparam;
  struct array_specification arrayspec;
  struct ion_model isphere;
  struct rfi_model rfimodel;
  struct stationspec *st1spec, *st2spec;
  struct stnbeam *beams=NULL;
  struct stnbeam *beam1=NULL;  /* Added by L. Benkevitch, 01-Dec-2010 */
  struct beamgrid bgrid;
  oob_beam_t     *oob_beams=NULL;
  source_model_t *ooblist=NULL;
  struct conv_fn beamconvl;
  struct iono_screen ionoscreen;
  struct uvgrid_parms uvgrid;
  struct uvpatch input_uv,conv_result;
  complex visibility[4];
  struct visfreq *visfreq=NULL;
  /* struct visblock visblock; */
  struct scaninfo *sc=NULL;
  int myid, numprocs;  /* Fake MPI declarations */
  extern char progname[]; 
  extern char *progstr[]; 
  struct bl_trespass badbl[N_TRESPASSED_BASELINES]; 
  long ib, ibbl = 0;    /* index into badbl */
  FILE *fh_bltab;   /* Where badbl is saved in case *
		     * of baseline exceptions.        */ 
  int baseline_error_flag = 0;  
  int baseline_table_overflow = 0;  
  struct stationspec *stn;
  double dx, dy;    /* Array station coordinates */
  double st_radius; /* A station distance from array centre */
  double umin_len, vmin_len, umax_len, vmax_len;
  double umin, vmin, umax, vmax;
  double freq_umin = 0, freq_vmin = 0, freq_umax = 0, freq_vmax = 0;
  double bl, blmax = 0;
  const double pi = 3.1414926, radeg = 180./pi, dtor = pi/180.;
  int log2N, log2N_1;
  double uvmax, fov_1;
  extern int msglev;

  /* Comparison function of station radii, *
   * used for station qsorting             */

  int station_cmp(const struct stationspec *st1, 
		  const struct stationspec *st2) {
    double dx1, dy1, st1_radius2, dx2, dy2, st2_radius2;
    dx1 = st1->east_offset;
    dy1 = st1->north_offset;
    st1_radius2 = pow(dx1,2) + pow(dy1,2); /* Distance^2 from array center */
    dx2 = st2->east_offset;
    dy2 = st2->north_offset;
    st2_radius2 = pow(dx2,2) + pow(dy2,2); /* Distance^2 from array center */
    /* return (int)(st2_radius2 - st1_radius2); */
    if      (st2_radius2 > st1_radius2) return  1;
    else if (st2_radius2 < st1_radius2) return -1;
    else                                return  0;
  }

  strcpy(progname, progstr[2]); /* progname = "CHECK_BASELINE" */

  printf("Entering check_baseline() =========================== \n");
  /* initialise stuff */
  ainit = 0.0;
  na = 0.0;
  for(i=0; i<NPOL*NPOL; i++) {
    beamconvl.convl[i] = NULL;
  }
  for(i=0; i<NPOL; i++) {
    input_uv.patch[i] = NULL;
    conv_result.patch[i] = NULL;
    uvgrid.fh[i]=NULL;
  }
  bgrid.grid = NULL;
  bgrid.size = 0;
  ionoscreen.phase_delay = NULL;
  ionoscreen.valid = NULL;
  ionoscreen.size_x = ionoscreen.size_y = 0;
  uvgrid.n_pol_inputs=1;
  uvgrid.ncells =0;
  uvgrid.padfactor=0;
  uvgrid.cellsize=0;
  arrayspec.array_name[0] = '\0';
  arrayspec.nst = arrayspec.num_pols=0;
  numprocs = 1;
  myid = 0;

  for (i=0; i<N_TRESPASSED_BASELINES; i++) {
    badbl[i].st1 = 0;      /* Station 1 */
    badbl[i].st2 = 0;      /* Station 2 */
    badbl[i].u_len = 0.0;  /* U in meters */
    badbl[i].v_len = 0.0;  /* V in meters */
    badbl[i].freq = 0.0;   /* (Minimum) frequency */
    badbl[i].ch = 0;       /* Frequency channel (IF) */
    badbl[i].cch = 0;      /* Correlator freq channels within an IF. */
    badbl[i].usize = 0;    /* Patch size along U  */
    badbl[i].vsize = 0;    /* Patch size along V */ 
    badbl[i].ercode = 0;   /* Error code as given in get_patch() */
  }  

  /* Process array site specification */
  if (get_site (vgparam.site, &arrayspec) != 0) {
    msg ("Failure processing site '%s'", 2, vgparam.site);
    exit (1);
  }

  /* Read array properties */
  if (read_array_spec (vgparam.arrayfilename, vgparam.layoutdirname, 
		       &arrayspec) != 0) {
    msg ("Failure reading array specification from file %s", 3, 
	 vgparam.arrayfilename);
    exit (1);
  }

  /*
   * Print out array coordinates
   * Added by L. Benkevitch, 06-Dec-2010
   */
  /* printf("\n==== Array Coordinates Unsorted (meters) =====\n"); */
  /* printf("St_ID        East        North       Radius \n"); */
  /* for (i = 0; i < arrayspec.nst; i++) { */
  /*   stn = &arrayspec.station[i]; */
  /*   dx = stn->east_offset; */
  /*   dy = stn->north_offset; */
  /*   st_radius = sqrt(pow(dx,2) + pow(dy,2)); /\* Dist. from center *\/ */
  /*   printf("%5d %12g %12g %12g\n", stn->st_id, dx, dy, st_radius); */
  /* } */

  /*
   * Sort stations in arrayspec in order of descendance of 
   * their distance from the array center
   * Added by L. Benkevitch, 06-Dec-2010
   */

  qsort (arrayspec.station, arrayspec.nst, sizeof(struct stationspec), 
	 station_cmp);

  /* Read in observation specification */
  msg("Observation file name: %s", 3, vgparam.obsfilename);
  if (read_obs_spec (vgparam.obsfilename, 
		     &obsparam,arrayspec.location.lon_radian) != 0) {
    msg ("Failure reading observation specification from file %s", 3, 
	 vgparam.obsfilename);
    exit (1);
  }
    
  /* We have to copy the number of input pols to the uvgrid structure. */
  /* This is found during the parse_cmdline call but the uvgrid structure */
  /* does not yet know. */

  uvgrid.n_pol_inputs = vgparam.n_pol_inputs;

  /* there really should only be one structure to perform this function */


  /* Open the main visibility grid file */
  /* FOV in radians */
  fov = obsparam.fov_RA / 206264.8062471;  // radians; i.e., (.)/RADASEC
  if (vgparam.do_vis && open_gridfile(vgparam.gridfile, fov, &uvgrid) != 0) {
    msg ("Failure opening main visibility grid file(s)", 3);
    exit (1);
  }

  /* discover how many polarisation products are required */
  if (arrayspec.num_pols < 1 || arrayspec.num_pols>2) {
    msg("ERROR: bad number of required output polarisations: %d",3,
	arrayspec.num_pols);
    exit(1);
  }
  vgparam.n_pol_products = arrayspec.num_pols*arrayspec.num_pols;
  msg("There are %d pols, %d pol products",1,arrayspec.num_pols,
      vgparam.n_pol_products);

  /* initialise the polarisation converter matrix */
  if (vgparam.n_pol_products == 1) {
    init_stokes_I_only();
  }
  else {
    init_stokes_2_xy();
  }

  /* Write out the station file */
  if (write_stnfile (&arrayspec, vgparam.simname) != 0) {
    msg ("Failure writing station output file", 3);
    exit (1);
  }


  /* Allocate memory for beams */
  beams = (struct stnbeam *)calloc(arrayspec.nst, sizeof (struct stnbeam));
  if (beams == NULL) {
    msg ("Failure allocating memory for station beams", 3);
    exit (1);
  }

  /* Read in out-of-beam sources */
  if (vgparam.do_oob) {
    if (read_oob(vgparam.oobfilename, &ooblist, &noob) != 0) {
      msg ("Failure reading oob sources from file %s", 3, 
	   vgparam.oobfilename);
      exit (1);
    }
    msg("there are %d OOB sources",1,noob);
    init_oob(noob, arrayspec.nst, arrayspec.num_pols*arrayspec.num_pols, 
	     &oob_beams);
  }
    
  /* Read in the RFI model */
  if (vgparam.do_rfimodel)
    if (read_rfimodel (vgparam.rfimodel_filename, &rfimodel) != 0)
      {
	msg ("Failure reading RFI model from file %s", 
	     3, vgparam.rfimodel_filename);
	exit (1);
      }
  if (do_acc) account ("Read input information");
    
  ///* Only node 0 writes main header	*/
  ///* note that node 0 refers to the first node in the list of machines */
  //if (myid == 0) write_main_header (&obsparam, &vgparam, &arrayspec,outfp);

  /* get current write position in output file    */
  /* needed for streaming outputs from the different nodes */
  refpoint = sizeof(struct main_header);

  /* Set up station beam sampling grid */
  if (vgparam.do_vis && 
      fill_beamgrid(&bgrid, uvgrid.cellsize, &obsparam) != 0) {
    msg ("ERROR: failure in fill_beamgrid", 3);
    exit (1);
  }
  /* Initialize the ionosphere. This must come after init of station beam 
   * sampling grid. The size/spatial location of the ionosphere grid is 
   * the same as the size/location of the beam sampling grid. */
  if (vgparam.do_ionomodel) {
    strncpy (isphere.iono_fname, vgparam.ionomodel_filename, IONO_MAX_FNAME);
    strncpy (isphere.settings_fname, vgparam.ionosettings_filename, 
	     IONO_MAX_FNAME);
    init_ionoscreen(arrayspec.nst,bgrid.size,bgrid.size,&ionoscreen);
    if (open_ionosphere(&isphere) != 0) {
      msg ("Failure reading ionospheric model from file %s",3, 
	   isphere.iono_fname);
      exit (1);
    }
  }
    
  /* Initialization of random seeds for noise */
  if (vgparam.do_noise) {
    struct timeval now;
    gettimeofday(&now,NULL);
    seed_noise ((int)now.tv_usec, myid);
  }


  /*******************************/
  /* ***** START ACTUAL WORK *****/
  /*******************************/
    
  nelev = 0;
  /*
   * Prepare min and max valies for search
   */
  umin_len = vmin_len =  DBL_MAX_10_EXP;
  umax_len = vmax_len = -DBL_MAX_10_EXP;
  blmax = 0.0;
  ; 
 
  /* Loop over scans */
  for (scan=0; scan<obsparam.nscan; scan++) {
    if (ibbl == N_TRESPASSED_BASELINES) break; //==========================>>>
    /* only the node 0 will print out the uvgrid info */
    if (vgparam.do_vis && (myid == 0)) {
      msg("Number of cells is : %d", 2, uvgrid.ncells);
      msg("Padfactor is : %.5f", 2, uvgrid.padfactor);
      msg("Cellsize  is : %.5f (wavelengths)", 2, uvgrid.cellsize);
    }

    ///* write timeblock marker, Only node 0 will write this */
    //if (myid == 0) fwrite (timblock, 1, 8, outfp);
    ///* Write timeblock header itself, Only node 0 will write this */
    //if (myid == 0) write_timeblock_header(&obsparam, scan, &vgparam, outfp);
    
    /* Convenience pointer */
    sc = obsparam.scan + scan;
    /* Update file position pointer to account  */
    /* for timblock and timeblock header writing    */
    refpoint += 8 + size_of_tbh(sc->nfreq); 
       
    /* file_offset sets output file pointer to proper position,	*/
    /* gets number of slices for this processor, starting	*/
    /* slice for this processor, and sets the refpoint to	*/
    /* the end of the scan	*/
    file_offset (numprocs, myid, &refpoint, &obsparam, scan, arrayspec.nst, 
    		 &intg_start, &intg_num, vgparam.n_pol_products,
    		 &write_pos);
    msg("num_slices=%d, start_slice=%d, start pointer=%lld, " \
	"final pointer=%lld", 2, intg_num, intg_start, 
	(long long) write_pos, (long long) refpoint);

    /* If no slices for this processor, 
     * skip out of scan loop */
    if (intg_num == 0) continue;

    /* Loop over integration times */
    /* at this point, the simulation is split into 
     * multiple chunks of contiguous 
     * time slices, and farmed out to 
     * different nodes in the cluster. note that 
     * intg_start (starting time slice) and 
     * intg_num (# of time slices per node) 
     * are node-dependent. in non-parallel mode 
     * (i.e., non-mpi mode), the entire 
     * simulation is done by the master node itself */

    scanstart = time_to_double(sc->start);

    /* timeflag = 0; */

    for (intg=intg_start; intg<(intg_start+intg_num); intg++) {
      if (ibbl == N_TRESPASSED_BASELINES) break; //==========================>>>
      start = scanstart + (intg * obsparam.integ_time);

      /* Only compute GHA once 
       * per integration */
      time = start + obsparam.integ_time/2.0;
      if (scan == 0)              // save time of 0th scan for
	time_scan0 = time;      // relative ionosphere times

      if (obsparam.scan[scan].gha_used) {
	gha = obsparam.scan[scan].gha_start + 
	  intg*obsparam.integ_time/3600.0*M_PI/12.0*(24.0/23.9344696);
      }
      else if (compute_gha(&obsparam, time, &gha) != 0) {
	msg ("Failure in compute_gha, skip this integration", 2);
	continue;
      }
      msg ("scan: %d, intg: %d GHA = %g (radian), lst = %g (radian)", 1,
	   scan,intg, gha,
	   gha+obsparam.phase_cent_RA+arrayspec.location.lon_radian);

      // if desired, center ionosphere such that
      // 0th scan time is at 2nd tabular point
      // of ionosphere model (for interpolation).
      if (vgparam.center_ion)
	reltime = time - time_scan0 + isphere.origin[T] + isphere.delta[T];
      else
	reltime = time;

      /* If ionosphere flag set then */
      /* Read 4D ionosphere file; interpolate
       * to 3D array for this time */
      if (vgparam.do_ionomodel) {
	    
	// if force time then use time from iono_settings else use reltime
	if (isphere.do_force == 1) {
	  time_ut = isphere.force_ut;
	  msg ("Time for ionosphere model forced to &g seconds UT", 2, 
	       time_ut);
	} else {
	      
	  // Calculate current time in UT given fact that scanstart 
	  // is seconds from jan 1 2000
	  time_ut = start-86400.0*floor(start/86400.0);
	  msg ("Time for ionosphere model &g seconds UT derived " \
	       "from start time", 2, time_ut);
	}
	    
	if (interp_t (&isphere, time_ut)) {
	  msg ("Failure interpolating time slice, skipping AP", 2);
	  continue;
	}
      }
	  
	  
      /* Initialize station beam flags	*/
      for (st1=0; st1<arrayspec.nst; st1++) beams[st1].valid_flag =0;
      /* if there are OOB sources, reset the station response since
	 the time has changed */
      if (vgparam.do_oob) {
	for(i=0; i<arrayspec.nst; i++) oob_beams[i].valid=0;
      }
      /* Loop over baselines */
      for (st1=0; st1<arrayspec.nst-1; st1++) {
	if (ibbl == N_TRESPASSED_BASELINES) break; //======================>>>
	for (st2=st1+1; st2<arrayspec.nst; st2++) {
	  if (ibbl == N_TRESPASSED_BASELINES) break; //====================>>>

	  msg ("baseline %d-%d", -2, st1,st2);
	  if (st2%500 == 0) 
	    msg ("baseline %d-%d", 0, st1,st2);
	  /* Convenience pointers */
	  st1spec = arrayspec.station + st1;
	  st2spec = arrayspec.stationation + st2;
	      
	  /* Compute u, v, w in meters */
	  elev_flag = compute_uvw(st1spec, st2spec, gha, &obsparam, 
				  &u_len, &v_len, &w_len);
		
	  /* Elevation limits */
	  /* for the sake of VLBI sims, we need to tolerate 
	   * some antennas being below their elevation
	   * limit, so the elevation flag for the baseline 
	   * doesn't make sense */
	  if (elev_flag) {
	    msg("WARNING: elevation out of range. Setting weights " \
		"for baseline %d-%d to zero.",1,st1,st2);
	  }
	      
	  /* Make ionosphere phase screen */
	  if (!elev_flag && vgparam.do_ionomodel && vgparam.do_vis) {
	    if (ionoscreen.valid[st1]==0) {
	      ret = compute_ionoscreen(st1spec, st1, &bgrid, gha,
	  			       obsparam.phase_cent_RA,
	  			       uvgrid.padfactor,
	  			       &isphere, vgparam.center_ion,
	  			       &ionoscreen);
	      if (ret != 0) {
	  	msg("compute_ionoscreen failed. exiting.",3);
	  	exit (1);
	      }
	    }
	    if (ionoscreen.valid[st2]==0) {
	      ret = compute_ionoscreen(st2spec, st2, &bgrid, gha,
	  			       obsparam.phase_cent_RA,
	  			       uvgrid.padfactor,
	  			       &isphere, vgparam.center_ion,
	  			       &ionoscreen);
	      if (ret != 0) {
	  	msg("compute_ionoscreen failed. exiting.",3);
	  	exit (1);
	      }
	    }
	  }
	  if (!elev_flag && vgparam.do_ionomodel && vgparam.do_oob) {
	    int result;
	    if (oob_beams[st1].valid_phase==0) {
	      result = compute_ionoscreen_oob(st1spec,st1,gha,
	  				      obsparam.phase_cent_RA,
	  				      noob,ooblist,&isphere,
	  				      vgparam.center_ion,
	  				      oob_beams);
	      if (result != 0) {
	  	msg("compute_ionoscreen_oob failed with code %d. exiting.",
	  	    3,result);
	  	exit (1);
	      }
	    }
	    if (oob_beams[st2].valid_phase==0) {
	      result = compute_ionoscreen_oob(st2spec,st2,gha,
	  				      obsparam.phase_cent_RA,
	  				      noob,ooblist,&isphere,
	  				      vgparam.center_ion,
	  				      oob_beams);
	      if (result != 0) {
	  	msg("compute_ionoscreen_oob failed with code %d. exiting.",
	  	    3,result);
	  	exit (1);
	      }
	    }
	  }

	  /* Loop over observing frequency IFs */
	  for (ch=0; ch<sc->nfreq; ch++) {
	    if (ibbl == N_TRESPASSED_BASELINES) break; //==================>>>
	    obschan_freq = sc->freq[ch].frequency;
	    bw = sc->freq[ch].bandwidth;
	    cbw = obsparam.corr_chan_bw;
	    ncch = ceil (bw / cbw);
            
	    /* Allocate and initialize visfreq struct */
	    bufsize = size_of_vf(ncch, vgparam.n_pol_products);
	    // fprintf(stderr,"bufsize is: %d. size for maps2uvfits: %d\n",
	    //   bufsize,sizeof(struct visgroup)+
	    //   (vgparam.n_pol_products-1)*sizeof(complex));
	    if (visfreq==NULL) visfreq = (struct visfreq *) calloc(1,bufsize);
	    visfreq->chan   = ch;
	    visfreq->ncchan = ncch;

	    /* Loop over correlator freq channels within an IF. */
	    for (cch=0; cch<ncch; cch++) {
	      if (ibbl == N_TRESPASSED_BASELINES) break; //================>>>
	      /* Center frequency of this channel */
	      freq = obschan_freq + (cch + 0.5) * cbw;

	      /* if do.vis is FALSE, then skip vis computation 
	       * from pixellised input sky */
	      if (!elev_flag && vgparam.do_vis) {
            
		/* Compute u,v in lambdas */
		u = u_len * freq * 1.0e6 / VLIGHT;
		v = v_len * freq * 1.0e6 / VLIGHT;
		w = w_len * freq * 1.0e6 / VLIGHT;

		//printf("u,v,w (visgen) = %g, %g, %g\n", u,v,w);
		/* Compute station beams for this time, freq, baseline */
		/* stn, beam, bgrid, freq, gha, dec_ref, ra_ref) */
		if (!beams[st1].valid_flag) {
		  beams[st1].station = st1;
		  ret = makebeam(st1spec, beams+st1, &bgrid, freq, gha, 
				 &obsparam);
		  if (ret != 0) {
		    fprintf(stderr,"ERROR: makebeam failed for " \
			    "station %d, freq %f, gha %f, point cent: %f,%f\n",
			    st1,freq,gha,obsparam.point_cent_DEC, 
			    obsparam.point_cent_RA);
		    exit(ret);
		  }
		}
		if (!beams[st2].valid_flag) {
		  beams[st2].station = st2;
		  ret = makebeam(st2spec, beams+st2, &bgrid, freq, gha, 
				 &obsparam);
		  if (ret != 0) {
		    fprintf(stderr,"ERROR: makebeam failed for " \
			    "station %d, freq %f, gha %f, point cent: %f,%f\n",
			    st2,freq,gha,obsparam.point_cent_DEC, 
			    obsparam.point_cent_RA);
		    exit(ret);
		  }
		}
		//if (do_acc) account ("Make station beams");
		msg("Making station beams",-2);
		/*
		 * Search for min and max u_len and v_len
		 * and the longest baseline blmax.
		 * This is required when all the baselines are fit
		 * within the uv-plane, all patches are small enough
		 * and hence the bad baseline table is not being 
		 * filled in.
		 * Added by L. Benkevitch 08-Dec-2010.                     
		 */
		if (umax_len < u_len) {
		  umax_len = u_len;
		  freq_umax = freq;
		}
		if (umin_len > u_len) {
		  umin_len = u_len;
		  freq_umin = freq;
		}
		if (vmax_len < v_len) {
		  vmax_len = v_len;
		  freq_vmax = freq;
		}
		if (vmin_len > v_len) {
		  vmin_len = v_len;
		  freq_vmin = freq;
		}
		bl = sqrt(pow(u,2) + pow(v,2));  /* baseline length */
		if (blmax < bl) blmax = bl;
		
  /* printf("umin_len, umax_len = %9g %9g (meters)\n", umin_len, umax_len); */
  /* printf("vmin_len, vmax_len = %9g %9g (meters)\n", vmin_len, vmax_len); */
		/*                                                         *
		 * Since we do not do actual simulation here, and only     *
		 * beamconvl.xsize and beamconvl.ysize are used in         *
		 * patch_size() below, the time-consuming convolution      *
		 * computations in compute_beamconvl() are not needed.     *
		 * That is why the call to it is commented out, and the    *
		 * x and y sizes are given the same values as in           *
		 * compute_beamconvl().                                    *
		 * Added by L. Benkevitch 01-Dec-2010.                     *
		 *                                                         */
		/* Make station beam convolving fn */
		/* compute_beamconvl(beams, st1, st2, &ionoscreen, freq,   */ 
		/*		     vgparam.do_ionomodel, &beamconvl);    */  
		beam1 = beams + st1;
		// beam2 = beams + st2;
		beamconvl.xsize = beam1->size;
		beamconvl.ysize = beam1->size;
		//if (do_acc) account ("Compute station beam conv fn");
                
		/* Figure out patch size */
		msg("Getting patch size",-2);
		patch_size(st1spec, st2spec, start, freq, &obsparam, scan, 
			   intg, uvgrid.cellsize, &beamconvl, &udim, &vdim);
		/* msg("Getting patch",-2); // Why waste precious time?.. */
		/* Snip out relevant part of input UV grids */
		
		/* Dummy routine instead of get_patch() to avoid actual *
		 * work with file and disk access                       *
		 * Replaced by L. Benkevitch, 01-Dec-2010               */

		ret = get_patch_dummy (&uvgrid, u, v, freq, udim, vdim, 
				       &input_uv);
		/*  return:  0=OK                             */
		/*           1=udim or vdim is not positive   */
		/*           2=u out of range                 */
		/*           3=v out of range                 */
		/*           4=excessive patch size requested */
		/*           5=bad read                       */
		/*                                            */
		/* msg("Got patch",-2); // Why waste precious time?..   */

		/*                                                      *
		 * We do not print an error message here. Instead, the  *
		 * error-associated data are stored on badbl table      *
		 */
		if (ret) { /* Non-zero, hence, ret is error code */
		  /* 
		   * Check if the pair (st1,st2) is already in the table *
		   */
		  int bl_in_table = 0;
		  for (ib = 0; ib < ibbl; ib++) {
		    if ((badbl[ib].st1 == st1) && (badbl[ib].st2 == st2)) {
		      bl_in_table = 1;
		      break; //=====>>>
		    }
		  }

		  if (bl_in_table == 0) {
		    badbl[ibbl].st1 = st1;     /* Station 1 */
		    badbl[ibbl].st2 = st2;     /* Station 2 */
		    badbl[ibbl].u_len = u_len; /* U in meters */
		    badbl[ibbl].v_len = v_len; /* V in meters */
		    badbl[ibbl].freq = freq;   /* (Minimum) frequency */
		    badbl[ibbl].gha = gha; /* Greenwich Hour Angle, radians */
		    badbl[ibbl].ch = ch;       /* IF channel */
		    badbl[ibbl].cch = cch;    /* Correlator freq chann. */
		    badbl[ibbl].usize = input_uv.usize; /* Patch U-size */
		    badbl[ibbl].vsize = input_uv.vsize; /* Patch V-size */
		    badbl[ibbl].ercode = ret; /* get_patch() error code */
		    ibbl++;
		    baseline_error_flag = 1;
		  }
		  if (ibbl == (long)N_TRESPASSED_BASELINES) {
		    msg("WARNING: Table of oversized baselines and/or " \
			"patches is full (%ld entries); exiting main loop "\
			"for analysis.", 2, ibbl);
		    baseline_table_overflow = 1;
		  }
		}
		if (ibbl == N_TRESPASSED_BASELINES) {
		  msg("Breaking innermost loop\n", 0);
		  break; //=============>>>
		}
		

	      } /* At end of do_vis if statement */
	      else {
		visibility[0] = visibility[1] = visibility[2] = 
		  visibility[3] = c_zero();
	      }


	      /* flag the data if one or more stations have 
	       * the phase center at bad elevation */
	      if (elev_flag) visfreq->visgroup[cch].weight = 0.0;
 
	      msg("Vis: BL: %d-%d. chan: %d, %g,%g %g,%g %g,%g, %g,%g",
		  -1, st1, st2, ch,
		  c_real(visibility[0]), c_imag(visibility[0]),
		  c_real(visibility[1]), c_imag(visibility[1]),
		  c_real(visibility[2]), c_imag(visibility[2]),
		  c_real(visibility[3]), c_imag(visibility[3]));
	      /* if the frequency has changed, then we should really 
	       * recompute the station beams. However, if this is 
	       * for a single freq then we can reuse the result 
	       * from before and make a significant saving
	       by not calling makebeam() */
	      if (ch > 2) {
		for (i=0; i<arrayspec.nst; i++) {
		  beams[i].valid_flag = 0;
		}
	      }
	    }  /* End correlator freq loop */
	  }    /* End observing freq channel loop */
	}      /* End baseline loop (st2) */
      }        /* End baseline loop (st1) */
    }          /* End integration time loop */
  }            /* End scan loop */
  
  if (baseline_error_flag && (msglev <= -1)) { 
    msg("\n", 3);
    msg("Bad baseline table ============================================\n", 3);
    msg("ERR  st1  st2       u_len         v_len        freq" \
	"        ch  cch      usize  vsize    radius", 3);
  
    for (ib = 0; ib < ibbl; ib++) {
      //    msg("%2d %4d %4d  %11.4g %11.4g %11.4g %4d %4d %6d %6d", 3, 
      msg("%2d %4d %4d  %11.4e   %11.4e   %11.4e   %4d  %4d    %6d   %6d %g",
	  3, 
	  badbl[ib].ercode,
	  badbl[ib].st1,
	  badbl[ib].st2,
	  badbl[ib].u_len,
	  badbl[ib].v_len,
	  badbl[ib].freq,
	  badbl[ib].ch,
	  badbl[ib].cch,
	  badbl[ib].usize,
	  badbl[ib].vsize,                  //);
	  sqrt(pow(badbl[ib].u_len,2) + pow(badbl[ib].v_len,2)));
    } 
  }
  if (baseline_table_overflow) 
    msg("WARNING: Table of oversized baselines and/or "			\
	"patches was filled up (%ld entries) before all of them "	\
	"had been brought up.", 2, ibbl);
  else
    msg("%d bad baselines found", 2, ibbl);

  /*
   * Print out sorted array coordinates
   */
  if (msglev <= -1) {
    printf("\n==== Array Coordinates Sorted (meters) =====\n");
    printf("St_ID        East        North       Radius \n");
    for (i = 0; i < arrayspec.nst; i++) {
      stn = &arrayspec.station[i];
      dx = stn->east_offset;
      dy = stn->north_offset;
      st_radius = sqrt(pow(dx,2) + pow(dy,2)); /* Distance from array center */
      printf("%5d %12g %12g %12g\n", stn->st_id, dx, dy, st_radius);
    }
  }
  /*
   * Find min and max u and v over the bad baselines
   */
  if (baseline_error_flag) { /* If all baselines are good, table is empty */
    umin_len = umax_len = badbl[0].u_len;
    vmin_len = vmax_len = badbl[0].v_len; 
    for (ib = 0; ib < ibbl; ib++) {
      if (umax_len < badbl[ib].u_len) {
	umax_len = badbl[ib].u_len;
	freq_umax = badbl[ib].freq;
      }
      if (umin_len > badbl[ib].u_len) {
	umin_len = badbl[ib].u_len;
	freq_umin = badbl[ib].freq;
      }
      if (vmax_len < badbl[ib].v_len) {
	vmax_len = badbl[ib].v_len;
	freq_vmax = badbl[ib].freq;
      }
      if (vmin_len > badbl[ib].v_len) {
	vmin_len = badbl[ib].v_len;
	freq_vmin = badbl[ib].freq;
      }
    }  
  }
  msg("umin_len, umax_len = %9g %9g (meters)", 2, umin_len, umax_len);
  msg("vmin_len, vmax_len = %9g %9g (meters)", 2, vmin_len, vmax_len);
  umin = umin_len * freq_umin * 1.0e6 / VLIGHT;
  umax = umax_len * freq_umax * 1.0e6 / VLIGHT;
  vmin = vmin_len * freq_vmin * 1.0e6 / VLIGHT;
  vmax = vmax_len * freq_vmax * 1.0e6 / VLIGHT;
  uvmax = fmax(fmax(abs(umin),abs(umax)), fmax(abs(vmin),abs(vmax))); 
  /* uvmax = fabs(umin); */
  /* if (uvmax < fabs(umax)) uvmax = fabs(umax); */
  /* if (uvmax < fabs(vmin)) uvmax = fabs(vmin); */
  /* if (uvmax < fabs(vmax)) uvmax = fabs(vmax); */

  /* (open_gridfile.c): uvgrid->cellsize = 1.0 / (fov * uvgrid->padfactor); */

  msg("umin, umax = %12g %12g (wl, at %g MHz)", 2, umin, umax, freq_umax);
  msg("vmin, vmax = %12g %12g (wl, at %g MHz)", 2, vmin, vmax, freq_vmax);
  msg("Absolute maximum dimensions along u or v: uvmax = %lg (wavelenghts)", 
      2, uvmax);
  msg("The longest baseline is %lg (wavelengths)", 2, blmax);
  msg("The array resolution is %lg (radians) or %lg (degrees) "\
      "at frequency %g MHz", 2, 1./blmax, radeg/blmax, freq);

  /*
   * Recommendations:
   * either decrease the uv-grid cellsize or shrink the field-of-view.
   */
  if (baseline_error_flag) {
    log2N = (int)log2((double)uvgrid.ncells);
    /* ncells_1 = 2.*uvmax/uvgrid.cellsize; */
    log2N_1 = (int)ceil(log2(2.*uvmax/uvgrid.cellsize));
    fov_1 = radeg*(0.5*uvgrid.ncells*uvgrid.padfactor/uvmax); /* in degrees */
    fov_1 = fov_1*3600; /* in arcseconds */
    msg(" ", 3);
    msg("Recommendations:", 3);
    msg("In order to continue simulation with this array,", 3);
    msg("EITHER increase the patameter 'Visibility grid x-log2N' in file", 3);
    msg("'%s%s' from current %ld to %ld, ", 3, 
	vgparam.simname, "_Description.txt", log2N, log2N_1);
    msg("and then run LOsim and MAPS_im2uv programs, ", 3);
    msg("OR decrease the parameters FOV_size in file", 3);
    msg("'%s' from current %lg down to %lg (arcseconds) or less.", 3, 
	vgparam.obsfilename, obsparam.fov_RA, fov_1);
    msg("'%s' Added for check.", 3, 
	vgparam.obsfilename, obsparam.fov_RA, fov_1);
  }

  printf("Exit check_baseline() ==============================\n");
  strcpy(progname, progstr[0]); /* progname = "VISGEN" */

  return baseline_error_flag;

}

