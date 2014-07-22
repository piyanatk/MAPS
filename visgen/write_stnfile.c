/************************************************************************/
/*                                                                      */
/* As soon as the array specification information is known, we can      */
/* write out a compact description for use by AIPS++.  This routine     */
/* does that job.                                                       */
/*                                                                      */
/*      Inputs:     arrayspec       Complete array description          */
/*                  simname         Root name of this simulation run    */
/*                                                                      */
/*      Output:     Filled output file with .stn extension              */
/*                  return value    0=OK, else bad                      */
/*                                                                      */
/* Created 12 Feb 2002 by CJL                                           */
/*                                                                      */
/* Modified 30 Jan 2003 by RB 						*/
/* 	for compatability with the MPI version of visgen		*/ 
/*									*/
/************************************************************************/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "array_specification.h"
#include "station.h"
#include "utility.h"

#ifdef USEMPI
#include "mpi.h"
#endif /* USEMPI */

/* mpi related declarations */
int numprocs, myid, namelen, partner;
#ifdef USEMPI
char processor_name[MPI_MAX_PROCESSOR_NAME];
#endif /* USEMPI */

/* expanded with mpi related declarations */
int
write_stnfile (/* arrayspec, simname) */
struct array_specification *arrayspec, char *simname){
  int n, size, ant, st,maxant=0;
  char filename[256];
  FILE *fp;
  struct station *stn;
  struct antenna *anten;

#ifdef USEMPI
  extern char visdir[];
  char mvisdir[100];
#endif /* USEMPI */

#ifdef USEMPI

/*  MPI routines */

/*  how many nodes are included in this cluster */
    MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
/*  get rank and name of the processor */
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);
    MPI_Get_processor_name(processor_name,&namelen);
#else
    numprocs = 1;
    myid = 0;
#endif /* USEMPI */

#ifdef USEMPI
/* pass around the environment variable VISDIR to the other nodes       */

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
      // strcpy(visdir,"/data2/simulator/visdata");
      if (myid == 0) { strcpy(mvisdir, visdir); }
#endif /* USEMPI */

                                        // Open the output file
    /*
    sprintf (filename, "%s/%s.stn", visdir, simname);
    */
    sprintf (filename, "%s.stn", simname);

    /* RB 13/10/02 temporary fix */
    /*
    strcpy (filename, "blanknoise.stn");
    */

    if ((fp = fopen (filename, "w")) == NULL){
      msg ("Unable to open output station file '%s'", 2, filename);
      return (1);
    }

    // Allocate enough memory for any station
    for (st=0; st < arrayspec->nst; st++) {
      if (arrayspec->station[st].layout->nant > maxant) maxant = arrayspec->station[st].layout->nant;
    }
    size = 2 * sizeof (int) + 3 * sizeof (double) 
                                + maxant * sizeof (struct ant_loc);
    if ((stn = (struct station *)malloc (size)) == NULL)
        {
        msg ("Memory allocation failure in write_stnfile", 2);
        return (1);
        }
                                        // Loop over stations
    for (st=0; st<arrayspec->nst; st++)
        {
        stn->station_no = st;
        stn->nant = arrayspec->station[st].layout->nant;
        stn->stn_x = arrayspec->station[st].x_coord;
        stn->stn_y = arrayspec->station[st].y_coord;
        stn->stn_z = arrayspec->station[st].z_coord;
        for (ant=0; ant<stn->nant; ant++)
            {
            anten = arrayspec->station[st].layout->ants + ant;
            stn->ants[ant].north = anten->north;
            stn->ants[ant].east = anten->east;
            stn->ants[ant].height = anten->height;
            }
                                        // Write out this station
        size = 2 * sizeof (int) + 3 * sizeof (double) 
                            + stn->nant * sizeof (struct ant_loc);
        n = fwrite (stn, size, 1, fp);
        if (n != 1)
            {
            msg ("Error writing out station %d", 2, st);
            free (stn);
            fclose (fp);
            return (1);
            }
        }

    free (stn);
    fclose (fp);

    msg ("Successfully wrote %d stations to file '%s'", 1, st, filename);

    return (0);
    }

