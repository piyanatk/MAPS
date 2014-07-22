/************************************************************************/
/*                                                                      */
/* This handles everything to do with the command line.  It identifies  */
/* acts upon various UNIX-style flags, and reads in control information */
/* from an external file.  Option flag arguments override the contents  */
/* of the external file, where applicable.                              */
/*                                                                      */
/*      Inputs:         argc, argv      Command line information        */
/*                                                                      */
/*      Output:         vgparam         Basic make_vis control info     */
/*                                                                      */
/* Created Jan 8th 2002 by CJL                                          */
/*                                                                      */
/* Modified 28 Jan 2003 by RB (for mpi compatability)                   */
/*         to pass around the environment VISDIR to the cluster nodes   */
/* Modified 13 Jan 2011 by J. Steeger and L. Benkevitch.                */
/*         Parameter 'FILE **asciifp' added for ascii visibility output */
/*                                                                      */
/************************************************************************/

#include <string.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>

#include "visgen.h"
#include "utility.h"

#ifdef USEMPI
#include "mpi.h"
#endif /* USEMPI */

#define FALSE 0
#define TRUE 1
#define DEFAULT_SEFD    500.0

void usage(char *filename);

/* new declarations for the mpi version */
int *buf, i, rank, namelen, nints, len, flag;
char *filename, *tmp;
int numprocs, myid, namelen, partner;
int opt_v_set = FALSE, opt_l_set = FALSE;
#ifdef USEMPI
char processor_name[MPI_MAX_PROCESSOR_NAME];
#endif /* USEMPI */
// MPI_File fh;  // -- NOT USED!
// MPI_Status status;  // -- NOT USED!

/* part of the original version; added mvisdir for mpi version */
int
parse_cmdline (/* argc, argv, vgparam, outfp) */
	       int argc,
	       char **argv,
	       struct vg_param *vgparam,
	       FILE **outfp, FILE **asciifp) 
{
  char c, outfilename[VISGEN_SIZE_NAME], asciioutfilename[VISGEN_SIZE_NAME];
  int array_flag, grid_flag, iono_flag, rfi_flag;
  int oob_flag, obs_flag, meta_flag;
  extern char *optarg;
  extern int do_acc, msglev;
  int npol_readin = 0;
#ifdef USEMPI
  extern char visdir[];
  char mvisdir[VISGEN_SIZE_NAME];
#endif /* USEMPI */
  char *mfilename, polfilename[VISGEN_SIZE_NAME],*textdir=NULL;
  FILE *mfile;
  enum polarisations { POLI = 0, POLQ, POLU, POLV };

  /* quick check command line args now and print usage message */
  if (argc < 2) usage(argv[0]);

  /*  MPI routines */
#ifdef USEMPI
  /*  how many nodes are included in this cluster */
  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
  /*  get rank and name of the processor */
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);
  MPI_Get_processor_name(processor_name,&namelen);
#else
  numprocs = 1;
  myid = 0;
#endif /* USEMPI */

  /* Initialize */
  vgparam->do_oob = FALSE;
  vgparam->do_ionomodel = FALSE;
  vgparam->do_rfimodel = FALSE;
  vgparam->do_noise = TRUE;
  vgparam->do_vis = TRUE;
  vgparam->obsfilename[0] = '\0';
  vgparam->simname[0] = '\0';
  vgparam->arrayfilename[0] = '\0';
  vgparam->oobfilename[0] = '\0';
  vgparam->ionomodel_filename[0] = '\0';
  vgparam->rfimodel_filename[0] = '\0';
  vgparam->site[0] = '\0';
  vgparam->center_ion = TRUE;
  vgparam->global_sefd = DEFAULT_SEFD;
  vgparam->n_pol_inputs = 1;
  vgparam->module_debug_flags = NULL;
  vgparam->vlbi = 0; 
  meta_flag = array_flag = grid_flag = iono_flag = oob_flag = 
    rfi_flag = obs_flag = FALSE;
  /* Make default simulation name */
  /*     default_name (vgparam->simname); */

  /* make default file for ionospheric settings based on "text" subdir */
  textdir=getenv("TEXTDIR");
  sprintf(vgparam->ionosettings_filename,"%s/iono_settings.txt",
	  textdir==NULL? ".": textdir);

  /* Interpret command line flags */
  while ((c = getopt (argc, argv, 
		      "am:n:vlp:ors:A:G:NM:ZV:I:i:O:R:S:D:")) != -1) 
  {
    switch (c) 
    {
      /* Profiling */
      case 'a':
        do_acc = TRUE;
        break;
        /* Verbosity control */
      case 'm':
        if (sscanf (optarg, "%d", &msglev) != 1)
        {
          msg ("Invalid -m flag argument '%s'", 2, optarg);
          msg ("Message level remains at %d", 2, msglev);
        }
        if (msglev > 3) msglev = 3;
        if (msglev < -3) msglev = -3;
        break;
        /* Name of this simulation */
      case 'n':
        strncpy (vgparam->simname, optarg, VISGEN_SIZE_NAME);
        break;
        /* Site identifier */
      case 's':
        strncpy (vgparam->site, optarg, VISGEN_SIZE_NAME);
        break;
        /* Process array as VLBI (global array) */
      case 'v':
	if (opt_l_set) {
	  msg("Conflicting options -v and -l. Exiting.", 3);
	  exit(1);
	}
	opt_v_set = TRUE;	
        vgparam->vlbi = 1;
        break;
        /* Process array as non-VLBI (localized array) */
      case 'l':
	if (opt_v_set) {
	  msg("Conflicting options -v and -l. Exiting.", 3);
	  exit(1);
	}
	opt_l_set = TRUE;
        vgparam->vlbi = -1;
        break;
	/* ASCII output file name */
      case 'p':
	if (strlen(optarg) >= VISGEN_SIZE_NAME) {
	  msg ("The ASCII output file name, '%s', is too long.\n", 2, 
	       asciioutfilename);
	  exit(1);
	}
	strcpy(asciioutfilename, optarg);
	*asciifp = fopen (asciioutfilename, "w");
	if (*asciifp == NULL) {
            msg ("Cannot open ascii output file %s\n", 2, asciioutfilename);
            exit(1);
	}
	break;
         /* Array description file */
      case 'A':
        array_flag = TRUE;
        strncpy (vgparam->arrayfilename, optarg, VISGEN_SIZE_NAME);
        break;
        /* UV grid file */
      case 'G':
        grid_flag = TRUE;
        vgparam->gridfile[npol_readin] = strdup(optarg);
        npol_readin++;
        vgparam->n_pol_inputs = npol_readin;
        break;
        /* File containing list of gridfiles */
      case 'M':
        {
          grid_flag = TRUE;
          mfilename = strdup(optarg);
          /* may as well read it in here but feel free to move it if it looks 
             out of place */
          mfile = fopen (mfilename, "r");
          if (mfile == NULL) {
            msg ("Invalid metafile (-M flag) %s, %s\n", 2, mfilename, 
		 strerror(errno));
            return(1);
          }
          else {
            while (!feof(mfile)) {
              if (fscanf(mfile,"%c %s \n",&c,polfilename) == 2) {
                /* re-using c */
                msg ("parsing %c and %s\n",2,c, polfilename);
                switch (c) {
                  case '#':
                    while (getc(mfile) != '\n') {
                        /* comment line */
                    }
                    break;
                  case 'I':
                    vgparam->gridfile[POLI] = strdup(polfilename);
                    msg ("readin %s\n",2,vgparam->gridfile[POLI]);
                    npol_readin++;
                    break;
                  case 'Q':
                    vgparam->gridfile[POLQ] = strdup(polfilename);
                    msg ("readin %s\n",2,vgparam->gridfile[POLQ]);
                    npol_readin++;
                    break;
                  case 'U':
                    vgparam->gridfile[POLU] = strdup(polfilename);
                    msg ("readin %s\n",2,vgparam->gridfile[POLU]);
                    npol_readin++;
                    break;
                  case 'V':
                    vgparam->gridfile[POLV] = strdup(polfilename);
                    msg ("readin %s\n",2,vgparam->gridfile[POLV]);
                    npol_readin++;
                    break;
                  default:
                    break;
                }
              }
            }
            vgparam->n_pol_inputs = npol_readin;
            fclose(mfile);
          }
        }
        break;
        /* Switch off noise calculations */
      case 'N':
        vgparam->do_noise = FALSE;
        break;
        /* Turn off visibility computation        */
      case 'Z':
        vgparam->do_vis = FALSE;
        break;
        /* Obs_spec file */
      case 'V':
        obs_flag = TRUE;
        strncpy (vgparam->obsfilename, optarg, VISGEN_SIZE_NAME);            
        break;

        /* Out of beam source list file */
      case 'O':
        oob_flag = TRUE;
            vgparam->do_oob = TRUE;
        strncpy (vgparam->oobfilename, optarg, VISGEN_SIZE_NAME);
        break;
        /* Ionospheric model file */
      case 'I':
        iono_flag = TRUE;
        vgparam->do_ionomodel = TRUE;
        strncpy (vgparam->ionomodel_filename, optarg, VISGEN_SIZE_NAME);
        break;
      case 'i':
        strncpy (vgparam->ionosettings_filename, optarg, VISGEN_SIZE_NAME);
        break;
        /* RFI model file */
      case 'R':
        rfi_flag = TRUE;
        vgparam->do_rfimodel = TRUE;
        strncpy (vgparam->rfimodel_filename, optarg, VISGEN_SIZE_NAME);
        break;
        /* set system-wide station SEFD */
      case 'S':
            vgparam->global_sefd = atof(optarg);
            break;
        /* enable module-level debugging */
      case 'D':
            vgparam->module_debug_flags = optarg;
            break;
      case '?':
      default:
        msg ("Unrecognized command line flag '%c'", 2, c);
        return (1);
    }
  }

  /* Check for complete input info */
  if (vgparam->n_pol_inputs < 1 || vgparam->n_pol_inputs > 4) 
  {
    msg ("ERROR: Unacceptable number of input polarisations: %d.", 2, 
	 vgparam->n_pol_inputs);
    return 1;
  }
  if (strlen (vgparam->simname) == 0)
  {
    msg ("Missing -n flag", 2);
    return (1);
  }
  if (strlen (vgparam->site) == 0) 
  {
    msg ("You must specify a site (control file or -s cmd line flag)", 2);
    return (1);
  }
  if (vgparam->obsfilename[0] == '\0') {
    fprintf(stderr,"ERROR: no observation specification file was given\n");
    exit(1);
  }
  if (vgparam->simname[0] == '\0') {
    fprintf(stderr,"ERROR: no output file name specified\n");
    exit(1);
  }

#ifdef USEMPI
  /* pass around the environment variable VISDIR to the other nodes        */

  // fprintf(stdout, " VISDIR = %s\n", visdir); fflush(stdout);
  //
  if (myid != 0)
  {
    MPI_Status stat;
    MPI_Recv(mvisdir, sizeof(mvisdir), MPI_BYTE,
        MPI_ANY_SOURCE, 2, MPI_COMM_WORLD, &stat);
    // // fprintf(stdout, "myid %d %s VISDIR = %s\n",
    //            myid, processor_name, mvisdir);
  }
  else
  {
    for (partner=1; partner< numprocs; partner++)
    {
      // strcpy(visdir,"/data2/simulator/visdata"); 
      strcpy(mvisdir, visdir);
      MPI_Send(mvisdir, strlen(mvisdir)+1, MPI_BYTE,
          partner, 2, MPI_COMM_WORLD);
    }
  }

  if (myid == 0) { strcpy(mvisdir, visdir); }
#endif /* USEMPI */

  /* Open output file */
  sprintf (outfilename, "%s.vis", vgparam->simname);
  if ((*outfp = fopen (outfilename, "w")) == NULL)  {
    msg ("Unable to open output file '%s'", 2, outfilename);
    return (1);
  }
  if (myid == 0) msg ("Opened output file '%s'", 2, outfilename);

  return (0);
}

void usage(char *filename) {
  fprintf(stdout,"Usage: %s [options]\n",filename);
  fprintf(stdout,"\t-n name\t name of the simulation\n");
  fprintf(stdout,"\t-s name\t name of the site\n");
  fprintf(stdout,"\t-A name\t name of file with the array definition " \
	  "(station types and locations)\n");
  fprintf(stdout,"\t-G name\t name of a UV grid data output file\n");
  fprintf(stdout,"\t\t If multiple input polarisations are required\n");
  fprintf(stdout,"\t\t specify multiple -G options.\n");
  fprintf(stdout,"\t\t Assumed order is I,Q,U,V\n\n");
  fprintf(stdout,"\t-M name\t name of a metafile listing gridfiles\n");
  fprintf(stdout,"\t\t Format: pol identifier [I,Q,U or V] and filename " \
	  "per line\n");
  fprintf(stdout,"\t\t [over-rides -G options]\n\n");

  fprintf(stdout,"\t-V name\t name of the observation specification file\n");
  fprintf(stdout,"\t-N\t do not add noise\n");
  fprintf(stdout,"\t-S SEFD\t set global station system equivalent flux " \
	  "density. Default: %g\n",DEFAULT_SEFD);  
  fprintf(stdout,"\t-Z\t do not compute visibilities from UV grid data " \
	  "(does not apply to oob)\n");
  fprintf(stdout,"\t-a\t turn on accounting\n");
  fprintf(stdout,"\t-I name\t name of the ionospheric model file. " \
	  "Turns on ionospheric modelling\n");
  fprintf(stdout,"\t-i name\t name of the ionospheric settings file. " \
	  "Default: $TEXTDIR/iono_settings.txt \n");
  fprintf(stdout,"\t-O name\t include 'out-of-beam' (oob) sources " \
	  "listed in file 'name'\n");
  fprintf(stdout,"\t-R name\t name of the RFI model file. " \
	  "Turns on RFI modelling (currently not implemented)\n");
  fprintf(stdout,"\t-m num\t turn debugging output to level num. " \
	  "(-2 lots thru 2 little)\n");
  fprintf(stdout,"\t-l\t treat the array as localized (non-VLBI). Only the "\
	  "baselines formed by the most remote stations are " \
	  "checked before the actual simulation.\n");
  fprintf(stdout,"\t-v\t treat the array as global (VLBI). All the " \
	  "baselines are checked before the actual simulation.\n");
  fprintf(stdout,"\t-p name\t name of a UV grid data ASCII output file\n");
  fprintf(stdout,"\t-D [IOMBN] turn on module-level dubugging " \
	  "(copious output) for one or many modules:\n");
  fprintf(stdout,"\t      \t I ionosphere\n");
  fprintf(stdout,"\t      \t O oob\n");
  fprintf(stdout,"\t      \t M makebeam\n");
  fprintf(stdout,"\t      \t B compute_beamconvl\n");
  fprintf(stdout,"\t      \t N add_noise\n");
  
  exit(0);
}
