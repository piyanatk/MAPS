/************************************************************************/
/*                                                                      */
/* This program implements the loops over time, baseline and frequency  */
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
/* Created Jan 7th 2002 by CJL                                          */
/*                                                                      */
/* Testing for MPI version, Oct 17th 2002, NDRB                         */
/* 09 Dec 2002                                                          */
/* Outputs from different nodes are streamed into a single file         */
/* 27 Dec 2002:                                                         */
/* version 2.2: Streaming the output into a single file,                */
/*                            but only for a _single_ time slice/node   */
/* version 2.3: Modified/tested for multiple time slices per node       */
/* 21 Jan 2003: Cleaned up for better readability                       */
/* 30 Jan 2003: Modified to set pointer offsets (*outfp) _only_ once    */
/*									*/
/* Latest (synchronized) version: 18 March 2003		NDRB@reber	*/
/*									*/
/* 24 Nov 2003: added more comments                                     */
/* 12 Apr 2010: conditionalized MPI code compilation                    */
/*  6 January 2011: added ASCII printout option (Jeremy Steeger)        */
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

/* additional header file for MPI versionget_patch.c */
#ifdef USEMPI
#include "mpi.h"
#endif /* USEMPI */

/* header files of the original visgen */
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

/* private function prototypes */
void init_module_debugging(const char *module_debug_flags);

/* global variables */
int do_acc = FALSE;

#define RADASEC 206264.80624709635515647 /* number of arcsec in a radian  */

/* #ifdef USEMPI */
/* char progname[] = "mpivisgen";  /\* extern variables for messages *\/ */
/* #else */
/* char progname[] = "visgen";     /\* extern variables for messages *\/ */
/* #endif /\* USEMPI *\/ */

char progname[] = "*******************************************************";

char *progstr[] = {
  "VISGEN",
  "MPIVISGEN",
  "CHECK_BASELINE"
};

int msglev = 2;
char version_no[] = "0.0";             /* Update this with new releases */

/* DONT EVER CHANGE the length of these or you will need to 
 * find all the places where 8 is
 hard coded into the various files. Good luck to you if you do. */
char timblock[8] = "timblck";
char vblock[8]   = "visblck";

#ifdef USEMPI
/* declarations for MPI version */
MPI_Status status;
#endif /* USEMPI */

/* this function is defined here because if we put it in ionosphere.c 
 * where it belongs, then it creates circular dependencies on header 
 * files which is too much fuss to fix...
 */
int compute_ionoscreen_oob(const struct stationspec *stn,
			   const int                 st_index,
			   const double              gha,
			   const double              fov_ra,
			   const int                 n_oob,
			   const source_model_t     *ooblist,
			   struct ion_model         *isphere,
			   const int                 center_iono,
			   oob_beam_t               *oob_beams)
{
  int src,result;
  double ha,az,el,xs[3],tec;

  msg("computing ionoscreen for oob sources for station %d\n",0,st_index);

  for(src=0; src<n_oob; src++) {
    ha = gha + stn->longitude + (fov_ra - ooblist[src].ra*(M_PI/12.0));
    /* calc alt/az. check source is above horizon */
    simple_E2H(ha,ooblist[src].dec*(M_PI/180.0),stn->latitude,&az,&el);
        
    msg("RA/DEC of source: %g,%g (hr,deg). Elevation of source: %g, " \
	"elevation limit of station: %g",0,
	ooblist[src].ra,ooblist[src].dec,el,stn->low_el_limit);
    if(el <= stn->low_el_limit) {
      continue;
    }
    xs[0] = xs[1] = xs[2] = 0;
    if (center_iono) {
      // offset relative to ion. center
      xs[0] = 1.0e-3 * stn->east_offset  + isphere->origin[Y] + 
	(isphere->delta[Y] * (isphere->dim[Y] - 1)) / 2;
      xs[1] = 1.0e-3 * stn->north_offset + isphere->origin[X] + 
	(isphere->delta[X] * (isphere->dim[X] - 1)) / 2;
      xs[2] = 0.0;
    }
    // find phase delay for this station via path integral 
    // through ionosphere model.
    result = path(isphere, xs, az, el, &tec);
    if (result) return result;

    // phase at 1 MHz, noting that the path integral was in km. 
    //(see Thompson et al. eqn. 13.138)
    oob_beams[st_index].phase[src] = -8.438e-10 * tec;
    msg("stn %d, oob source %d. az/el: %g,%g. tec: %g, phase: %g\n",0,
	st_index,src,az,el,tec,oob_beams[st_index].phase[src]);
  }
  oob_beams[st_index].valid_phase=1;
  return 0;
}

/*=========================================================*/
/*                        main()                           */          
/*=========================================================*/

int main (/* argc, argv) */
	  int argc,
	  char **argv)
{

  /* More MPI declarations	*/   int myid, numprocs;
#ifdef USEMPI
  int namelen;
  char processor_name[MPI_MAX_PROCESSOR_NAME];
#endif /* USEMPI */

  off_t refpoint, write_pos;

  /* part of the original visgen */
  int  i, j, scan=0, intg=0,st1, st2, ch, ncch=0, cch;
  int  bufsize=0, noob=0, nelev=0, ret=0, timeflag; 
  int  intg_start, intg_num=0, elev_flag;
  double start, time, bw, cbw, obschan_freq, freq;
  double u, v,w, u_len, v_len, w_len;
  double scanstart, gha, fov;  // gha in radian. 
  double udim, vdim,reltime=0, time_ut=0.0,time_scan0=0;
  double ainit, na, umean, vmean;
  FILE *outfp = NULL, *asciifp = NULL, *visfp = NULL;
  struct vg_param vgparam;
  struct observing_param obsparam;
  struct array_specification arrayspec;
  struct ion_model isphere;
  struct rfi_model rfimodel;
  struct stationspec *st1spec, *st2spec;
  struct stnbeam *beams=NULL;
  struct beamgrid bgrid;
  oob_beam_t     *oob_beams=NULL;
  source_model_t *ooblist=NULL;
  struct conv_fn beamconvl;
  struct iono_screen ionoscreen;
  struct uvgrid_parms uvgrid;
  struct uvpatch input_uv,conv_result;
  complex visibility[4];
  struct visfreq *visfreq=NULL;
  struct visblock visblock;
  struct scaninfo *sc=NULL;
  int ichkbl = 0;  /* Results of the check_baseline() call */   
  /* int ib; */
  int vlbi = FALSE; /* TRUE, if the array is distributed globally (-v key) */
  /********jsteeger variables*/
  int itime, minute_print=0, hour_print=0, day_print=0, year_print=0;
  double time_min, visibility_amp, visibility_phase, v_weight; 
  double time_print=0, second_print=0, printout[8];
  float vnoise;
  /********end of jsteeger variables*/

  visfp = fopen("vis.txt", "w");

#ifdef USEMPI
  strcpy((char *)progname, progstr[1]); /* extern variable for messages */
#else
  strcpy((char *)progname, progstr[0]); /* extern variable for messages */
#endif /* USEMPI */
  /* //progname = progstr[0];    /\* extern variables for messages *\/ */
  /* msg("++++++++++++++ Start _++++++++++++++\n",2); */
  /* printf("VISGEN: progname = %s, ptr = %p\n", progname, progname); */
  /* printf("progstr[0,1,2] = ptr = %p %p %p\n",  */
  /* 	   progstr[0], progstr[1], progstr[2]); */
  /* printf("progstr[0,1,2] = ptr = %s %s %s\n",  */
  /* 	   progstr[0], progstr[1], progstr[2]); */

  /* initialise stuff */
  ainit = 0.0;
  na = 0.0;
  for(i=0; i<NPOL*NPOL; i++) {
    beamconvl.convl[i] = NULL;
  }
  for(i=0; i<NPOL; i++) {
    input_uv.patch[i] = NULL;
    conv_result.patch[i] = NULL;
    uvgrid.fh[i]=NULL;  /* Brightness image file handles (for upto 4 pols) */
    vgparam.gridfile[i]=NULL;
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
  vgparam.n_pol_products = 1;
  vgparam.global_sefd = 0;
    
  arrayspec.array_name[0] = '\0';
  arrayspec.nst = arrayspec.num_pols=0;
    
  environment();

#ifdef USEMPI
  MPI_Init(&argc,&argv);		/* initialize MPI */
					/* how many nodes are included 
					 * in this cluster? */
  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
  /* get rank and name of the processor */
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);
  MPI_Get_processor_name(processor_name,&namelen);
#else
  numprocs = 1;
  myid = 0;
#endif /* USEMPI */

  strncpy(vgparam.layoutdirname,getenv("SIM"),sizeof(vgparam.layoutdirname));
  if ( vgparam.layoutdirname == NULL) {
    msg("SIM environment variable is not set. exiting\n",2);
    exit(1);
  }
  strncat(vgparam.layoutdirname,"/stn_layout",sizeof(vgparam.layoutdirname));
    
  /* Start accounting */
  account ("!BEGIN");

  /* Get command line arguments */
  if (parse_cmdline (argc, argv, &vgparam, &outfp, &asciifp) != 0) {
    msg ("Fatal error interpreting command line arguments", 3);
    exit(1);
  }

  /* Execute a visgen dry-run to check we will remain
   * within the UV-plane limits during the actual computations 
   * Added by L. Benkevitch*/

  ichkbl = check_baseline(vgparam); /*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
  if (ichkbl) {
    msg("Did not pass the check -- exiting.",3);
    exit (1);
  }

  /* enable module-level debugging if there is only one process */
  if (vgparam.module_debug_flags != NULL) {
    if (numprocs ==1) init_module_debugging(vgparam.module_debug_flags);
    else msg("Cannot enable module level debugging in MPI mode.",3);
  }

  /* Process array site specification */
  if (get_site(vgparam.site, &arrayspec) != 0) {
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
   * Setup the VLBI flag
   */
  if (vgparam.vlbi == 1) vlbi = TRUE;         /* cmd option -v */
  else if (vgparam.vlbi == 0) vlbi = arrayspec.vlbi; /* no cmd option */
  else if (vgparam.vlbi == -1)  vlbi = FALSE; /* cmd option -l */
  else {
    msg("ERROR: Wrong vgparam.vlbi value; may only be -1;0;1",3);
    exit(1);
  }
  if (vlbi) msg("The array is processed as VLBI",3);
  else      msg("The array is processed as non-VLBI",3);

  /* Read in observation specification */
  if (read_obs_spec (vgparam.obsfilename, 
		     &obsparam, &vgparam, arrayspec.location.lon_radian) != 0) {
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
    
  /* Only node 0 writes main header	*/
  /* note that node 0 refers to the first node in the list of machines */
  if (myid == 0) write_main_header (&obsparam, &vgparam, &arrayspec,outfp);

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
    
  /* Loop over scans */
  for (scan=0; scan<obsparam.nscan; scan++) {
    /* only the node 0 will print out the uvgrid info */
    if (vgparam.do_vis && (myid == 0)) {
      msg("Number of cells is : %d",2,uvgrid.ncells);
      msg("Padfactor is : %.5f",2,uvgrid.padfactor);
      msg("Cellsize  is : %.5f",2,uvgrid.cellsize);
    }

    /* write timeblock marker, Only node 0 will write this */
    if (myid == 0) fwrite (timblock, 1, 8, outfp);
    /* Write timeblock header itself, Only node 0 will write this */
    if (myid == 0) write_timeblock_header(&obsparam, scan, &vgparam, outfp);
    
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

    /* Added by L. Benkevitch because file_offset() does not      *
     * make output file positioning any more; it only calculates  *
     * the position (in write_pos).                               */ 
    fseeko (outfp, write_pos, SEEK_SET);

    msg("num_slices=%d, start_slice=%d, start pointer=%lld, "	\
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

    /*****NOTE:  Opening file here for jsteeger file printout.*/
    /* asciifp = fopen("visgenUVdata.txt", "w"); */
    if (asciifp) {
      fprintf(asciifp, "UV data taken from visgen.\n");
      fprintf(asciifp, "Scan Start (UT)            U(klam)      V(klam)");
      fprintf(asciifp, "      W(klam)  Baseline   Channel");
      fprintf(asciifp, "         Visibility (amp, phase)");
      fprintf(asciifp, "   Weight      Sigma\n");
    }
    /*****End of jsteeger changes.*/


    scanstart = time_to_double(sc->start);

    timeflag = 0;
    /* Loop over integration time slices */
    for (intg=intg_start; intg<(intg_start+intg_num); intg++) {
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
      for (st1=0; st1<arrayspec.nst-1; st1++)
	for (st2=st1+1; st2<arrayspec.nst; st2++) {
	  msg ("baseline %d-%d", -2, st1,st2);
	  if (st2%500 == 0) 
	    msg ("baseline %d-%d", 0, st1,st2);
	  /* Convenience pointers */
	  st1spec = arrayspec.station + st1;
	  st2spec = arrayspec.station + st2;
	      
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

	  /* Initialize visblock struct */
	  visblock.intg_no = intg;
	  visblock.station1 = st1;
	  visblock.station2 = st2;
	  visblock.nfreq = sc->nfreq;
	  visblock.u = u_len; 
	  visblock.v = v_len; 
	  visblock.w = w_len;
	  /* Notify reader that visblock 
	   * is next */
	  fwrite (vblock, 1, 8, outfp);
	  /* Write out visblock header */
	  fwrite (&visblock, 
		  sizeof(struct visblock)-sizeof(struct visfreq), 
		  1, outfp);

	  /* Loop over observing frequency IFs */
	  for (ch=0; ch<sc->nfreq; ch++) {
            obschan_freq = sc->freq[ch].frequency;
            bw = sc->freq[ch].bandwidth;
            cbw = obsparam.corr_chan_bw;
            ncch = ceil (bw / cbw);
            
            /* Allocate and initialize visfreq struct */
            bufsize = size_of_vf(ncch,vgparam.n_pol_products);
            // fprintf(stderr,"bufsize is: %d. size for maps2uvfits: %d\n",
	    //   bufsize,sizeof(struct visgroup)+
	    //   (vgparam.n_pol_products-1)*sizeof(complex));
            if (visfreq==NULL) visfreq = (struct visfreq *) calloc(1,bufsize);
            visfreq->chan   = ch;
            visfreq->ncchan = ncch;

            /* Loop over correlator freq channels within an IF. */
            for (cch=0; cch<ncch; cch++) {

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
                if (do_acc) account ("Make station beams");
                msg("Making station beams",-2);
                /* Make station beam convolving fn */
                compute_beamconvl(beams, st1, st2, &ionoscreen, freq, 
				  vgparam.do_ionomodel, &beamconvl);
                if (do_acc) account ("Compute station beam conv fn");
                
                /* Figure out patch size */
                msg("Getting patch size",-2);
                patch_size(st1spec, st2spec, start, freq, &obsparam, scan, intg,
			   uvgrid.cellsize, &beamconvl, &udim, &vdim);
		/* printf("Patch sizes in cells: udim = %g, vdim = %g cells\n",
		 *        udim/uvgrid.cellsize, vdim/uvgrid.cellsize); */
                msg("Getting patch",-2);
                /* Snip out relevant part of input UV grids */
                ret = get_patch (&uvgrid, u, v, freq, udim, vdim, &input_uv);
		if (vgparam.do_vis && (myid == 0)) {
		  msg("Patch sizes: udim = %g,  vdim = %g (wavelengths)", -2, 
		      udim, vdim);
		}
                msg("Got patch",-2);
                if (ret != 0) {
		  fprintf(stderr,"ERROR: get_patch() failed " \
			  "with code %d. u,v: %g,%g. udim,vdim: %g,%g\n",
			  ret,u,v,udim,vdim);
		  exit (1);
                }
                if (do_acc) account ("get patch");
                
                /* copy across relevant UV grid size params to 
		 * convolution result structure. This is a bit clunky
                 * but it is necessary to have the output data struct 
		 * be different to in the input for the full
                 * polarised case. There can be 4 output products 
		 * even with only one input Stokes product */
                conv_result.frequency = input_uv.frequency;
                conv_result.umin = input_uv.umin;
                conv_result.vmin = input_uv.vmin;
                conv_result.usize = input_uv.usize;
                conv_result.vsize = input_uv.vsize;
                conv_result.cellsize = input_uv.cellsize;

                /* Perform mini-convolution of patch 'input_uv' and 
		 * baseline beam 'beamconvl' to save result in 'conv_result'.
		 * Currently, beamconvl.xsize and ysize are
		 * defined as NUM_SAMPS=32 in fill_beamgrid.c 
		 * The size of conv_result is the same as that of input_uv */
                ret = do_convl(&input_uv, &conv_result, &beamconvl, 
			       vgparam.n_pol_products);
                if (ret != 0) {
		  fprintf(stderr,"ERROR: do_convl() failed\n");
		  exit (1);
                }
                if (do_acc) account ("convolve");
                
                /* Integrate over time/frequency cell */
                ret = integrate(&conv_result, st1spec, st2spec, time, 
				freq, scan, intg, &obsparam, visibility,
                                &umean, &vmean, vgparam.n_pol_products);
                if (ret != 0) {
		  fprintf(stderr,"ERROR: integrate() failed\n");
		  exit (1);
                }
                if (do_acc) account ("integrate");

              } /* At end of do_vis if statement */
              else {
                visibility[0] = visibility[1] = visibility[2] = 
		  visibility[3] = c_zero();
              }

              /* Add system noise to the visibility */
              /* note that seeds are node-dependent */
              if (vgparam.do_noise) {
                msg("Adding noise",-2);
                add_noise(st1spec, st2spec, beams+st1, beams+st2, 
			  obsparam.corr_chan_bw*1e6, obsparam.integ_time,
			  &uvgrid, vgparam.global_sefd, 1.0, 
			  visibility, vgparam.n_pol_products, &vnoise);
              }

              /* Add "out of beam" point sources to the integration */
              if (vgparam.do_oob) {
		oob_integrate(&obsparam,&arrayspec,gha,freq,st1,st2,
			      ooblist,noob,oob_beams,visibility);
              }

              /* Append to visibility freq block for this baseline/time */
              visfreq->visgroup[cch].cchan = cch;
              visfreq->visgroup[cch].flag = 0;
              visfreq->visgroup[cch].weight = 
		cbw * obsparam.integ_time* pow(st1spec->layout->nant * 
					       st2spec->layout->nant, 0.5);

              /* flag the data if one or more stations have 
	       * the phase center at bad elevation */
              if (elev_flag) visfreq->visgroup[cch].weight = 0.0;
 
	      /***********jsteeger addition*/ 
	      v_weight = visfreq->visgroup[cch].weight;
	      /***********end jsteeger addition*/

              for (i=0; i<vgparam.n_pol_products; i++) {
                visfreq->visgroup[cch].vis[i] = visibility[i];
              }
              msg("Vis: BL: %d-%d. chan: %d, %g,%g %g,%g %g,%g, %g,%g",
		  -1,st1,st2,ch,
                  c_real(visibility[0]),c_imag(visibility[0]),
		  c_real(visibility[1]),c_imag(visibility[1]),
                  c_real(visibility[2]),c_imag(visibility[2]),
		  c_real(visibility[3]),c_imag(visibility[3]));
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
            }               /* End correlator freq loop for (cch = 0; ... */

            /* Write out visfreq */
            fwrite (visfreq, bufsize, 1, outfp);
            free (visfreq);
            visfreq=NULL;
            if (do_acc) account ("Write datafile");

	    /*******This section added by jsteeger.  
	     * It prints out a file with all the UV data.*/
	    /* if (asciifp) { /\* i.e. if the file is open *\/ */
	    /*   visibility_amp = 0; */
	    /*   visibility_phase = 0; */
	    /* for (j=0; j<=3; j++) { */
	    /* 	if (c_real(visibility[j]) != 0 &&  */
	    /*            c_imag(visibility[j]) != 0) { */
	    /* 	  visibility_amp = sqrt(c_real(visibility[j])* */
	    /* 				c_real(visibility[j]) +  */
	    /* 				c_imag(visibility[j])* */
	    /* 				c_imag(visibility[j])); */
	    /* 	  visibility_phase = -1*atan2(c_imag(visibility[j]),  */
	    /* 				      c_real(visibility[j])); */
	    /* 	} */
	    /* } */

	    /* L. Benkevitch: */
	    if (asciifp) { /* i.e. if the file is open */
	      complex vis = visibility[0]; /* I Stokes component only! */
	      visibility_amp = sqrt(c_real(vis)*c_real(vis) + 
				    c_imag(vis)*c_imag(vis));
	      visibility_phase = atan2(c_imag(vis), c_real(vis));
	      
	      /*Clear printout array.*/
	      for (j=0; j<9; j++) {
		printout[j] = 0;
	      }

	      /*Find current time.*/
	      time_print = start;
	      second_print = fmod (time_print, 60.0);   /* Floating seconds */
	      time_min = (time_print - second_print) / 60.0; /* t in minutes */
	      itime = rint (time_min);   /* Convert to integer for the rest */
	      minute_print = itime%60;
	      itime /= 60;               /* Time now in hours */
	      hour_print = itime%24;
	      itime /= 24;                /* Time now in days */
	      /* for (qs=0; i<99; i++) Jeremy, what is "qs"? */ 
	      for (i=0; i<99; i++) 
		{                                 /* Clumsy but effective */
		  if (i%4 == 0 && itime <= 366) break;
		  else if (itime <= 365) break;
		  if (i%4 == 0) itime -= 366;
		  else itime -= 365;
		}
	      year_print = i + 2000;
	      day_print = itime + 1;              /* Days are 1-relative */

	      itime = 0;
	      time_min = 0;
			
	      /*Data order for printout array:
		printout[0] = UU (in meters)
		printout[1] = VV (in meters)
		printout[2] = WW (in meters)
		printout[3] = baseline1
		printout[4] = baseline2
		printout[5] = channel
		printout[6] = visibility amplitude
		printout[7] = visibility phase
	      */

	      printout[0] = (u_len * freq * 1.0e6 / VLIGHT)/1000.;
	      printout[1] = (v_len * freq * 1.0e6 / VLIGHT)/1000.;
	      printout[2] = (w_len * freq * 1.0e6 / VLIGHT)/1000.;
	      printout[3] = st1+1;
	      printout[4] = st2+1;
	      printout[5] = ch;
	      printout[6] = visibility_amp;
	      printout[7] = visibility_phase*180/3.14159265;
	      if (v_weight != 0) { /* Do not write lines with 0 weight */
		if (second_print<10) {
		  fprintf(asciifp, "%d:%03d:%02d:%02d:0%-5.2f ",
			  year_print,day_print,hour_print,
			  minute_print,second_print);
		} /* if (second_print<10) { */
		else {
		  fprintf(asciifp, "%d:%03d:%02d:%02d:%-6.2f ",
			  year_print,day_print,hour_print,
			  minute_print,second_print);
		} /* end of else */
		fprintf(visfp, "%12.6f\t%12.6f\n", c_real(vis), c_imag(vis));
		fprintf(asciifp, "% 12.2f % 12.2f % 12.2f  %02d-%02d      %02d"\
			"        (% 13f, % 12f)% 9g %10.5f\n",		\
			printout[0], printout[1], printout[2],
			(int)printout[3], (int)printout[4], (int)printout[5],
			printout[6], printout[7],
			v_weight, vnoise);
	      } /* if (v_weight != 0) { */
	    } /* if (asciifp)  */
	    /******* End of jsteeger changes.*/

	  }                   /* End observing freq channel loop */
          
          
	}                       /* End baseline loop */
      timeflag++;
    }                           /* End integration time loop */
    /* since the time has changed, all the ionospheric phase screens 
     * should be re-computed */
    if (vgparam.do_ionomodel) {
      for (st1=0; st1<ionoscreen.num_stations; st1++) 
	ionoscreen.valid[st1] =0;
      if (oob_beams!=NULL) {
	for (st1=0; st1<noob; st1++) oob_beams[st1].valid_phase =0;
      }
    }
  }    /* End scan loop for (scan = ...*/
    
  /* only node 0 prints profiling -- for debugging... */
#ifdef USEMPI
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);
#endif /* USEMPI */

  if (myid == 0) {
    if (do_acc) {
      printf("myid is %d \n",myid); fflush(stdout);
      account ("!REPORT");
    }
    account ("!SYSTEM");
  }

  /*****Closing file previously opened for printout.*/
  if (asciifp) fclose(asciifp);
  if (visfp) fclose(visfp);
 
  /*****End of jsteeger changes.*/

  msg("All done. exiting",2);

#ifdef USEMPI
  /* MPI calls */
  /* wait until everyone finishes 
   * their tasks */
  MPI_Barrier(MPI_COMM_WORLD);
  /* finalize the mpi mode */
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);
  if (myid == 0) fclose (outfp);

  MPI_Finalize();
#else
  fclose (outfp);
#endif /* USEMPI */

  exit(ret);
}
    

void init_module_debugging(const char *module_debug_flags) {
  int c;
  while ( (c = toupper(*(module_debug_flags++))) != '\0') {
    if (c=='O') {
      oob_set_debug(1);
      msg("Enabling OOB module debugging",2);
    }
    if (c=='M') {
      makebeam_set_debug(1);
      msg("Enabling makebeam module debugging",2);
    }
    if (c=='I') {
      ionosphere_set_debug(1);
      msg("Enabling ionosphere module debugging",2);
    }
    if (c=='B') {
      compute_beamconvl_set_debug(1);
      msg("Enabling compute_beamconvl module debugging",2);
    }
    if (c=='N') {
      addnoise_set_debug(1);
      msg("Enabling addnoise module debugging",2);
    }
  }
}
