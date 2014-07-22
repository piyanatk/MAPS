
// CHECK fseeko() parameters!


/************************************************************************/
/*                                                                      */
/* file_offset.c                                                        */
/*                                                                      */
/* Given number of processors, idenity of node, and parameters of       */
/* observation and array, this routine figures out where in the         */
/* output file to write node output.  An FSEEK command is issued        */
/* which positions the output file pointer to the correct write         */
/* position.  This routine also computes the next position at which     */
/* node=0 will begin writing and updates refpoint for the next scan.    */
/*                                                                      */
/*      Inputs:     numproc     number of processors in cluster         */
/*                  myid        identity of node calling this routine   */
/*                  refpoint    point at which node=0 starts writing    */
/*                              (measured in bytes)                     */    
/*                  obsparam    Observing parameter structure           */
/*                  scan_no     scan number                             */
/*                  outfp       output file pointer                     */
/*                  nstat       number of stations in configuration     */
/*		    intg_start 	starting slice for this processor       */
/*		    intg_num	number of slices for this processor     */	
/*                                                                      */
/*      Output:     NONE                                                */
/*                                                                      */
/* Created 30 Jan 03 by SSD/RB                                          */
/************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>

#ifdef USEMPI
#include "mpi.h"
#endif /* USEMPI */

#include "utility.h"
#include "observing_param.h"
#include "sizes.h"
#include "visibility.h"

int numprocs, myid, namelen, partner;
#ifdef USEMPI
char processor_name[MPI_MAX_PROCESSOR_NAME];
#endif /* USEMPI */

static int debug=0;

int
file_offset (/* numproc, myid, refpoint, obsparam, scan_no, nstat, outfp, 
	      * intg_start, intg_num, write_pos) */
	     int numproc,
	     int myid,
	     off_t *refpoint,
	     struct observing_param *obsparam,
	     int scan_no,
	     int nstat,
	     FILE *outfp,
	     int *intg_start,
	     int *intg_num,
	     int n_pol_products,
	     off_t *write_pos)
{
    int n,ch,ncch,nintg_scan,nintg_proc, nbase;
    int start, stop;
    long  bufsize;
    struct scaninfo *sc;
    double bw, cbw;
    div_t num_slice;

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

#ifdef USEMPI
    if (debug) printf("FILE offset for myid %d node %s ",myid, processor_name);
#endif /* USEMPI */
    
            /* set up convenience pointer   */
    sc = obsparam->scan + scan_no;

    cbw = obsparam->corr_chan_bw; // printf("cbw=%e\n",cbw);
    nbase = nstat*(nstat - 1)/2;
    bufsize = 0;

            /* Loop over observing frequency channels   */
            /* to determine total size of visfreq data per baseline.    */
    for (ch=0; ch<sc->nfreq; ch++) {
      bw = sc->freq[ch].bandwidth; 
      if(debug) printf("bw=%e, bw/cbw=%e, ceil(bw/cbw)=%e\n",
		       bw, bw/cbw, ceil(bw/cbw));
      ncch = (int) ceil (bw/cbw); 
      if(debug) printf("ncch=%d\n",ncch);
      bufsize+=size_of_vf (ncch, n_pol_products);
    }
    if(debug) printf ("bufsize = %ld\n",bufsize);

    
            /* Total number of integration times in scan    */
    nintg_scan = ceil(sc->duration / obsparam->integ_time);
    num_slice = div (nintg_scan, numproc);
    if(debug) printf ("results of div: quot=%d   rem=%d\n", 
		      num_slice.quot, num_slice.rem);

            /* Loop over processors to get offset for myid  */
            /* and total offset for this scan   */
    start = 0;
    if(debug) printf ("number of processors = %d\n", numproc);
    for (n=0;n<numproc;n++)
        {
	    /* Figure out how many time slices for this processor   */
        nintg_proc = num_slice.quot;
        if (n < num_slice.rem) nintg_proc += 1;

            /* If no time slices for this processor, skip to end        */
        if (nintg_proc == 0) {
	  if (myid == n) {*intg_num = 0; *write_pos = *refpoint;}
	  if(debug) printf ("myid %d n = %d,  position = %lld,  " \
			    "slices = %d\n  ", myid, n, (long long)*refpoint, 
			    nintg_proc);
	  continue;
	    }

            /* set write pointer to proper point in file for this processor */
	    /* and export start slice and number of slices 	*/
        if (myid == n) {
	  // Shep: why did you comment this out?? 
	  fseeko (outfp, *refpoint, SEEK_SET);
	  *write_pos = *refpoint;
	  *intg_start= start;
	  *intg_num  = nintg_proc;
	}

	    /* Figure out stop integration and update start integration	*/
        stop = start + nintg_proc;

        if(debug){
	  printf ("n = %d,  position = %lld,  slices = %d,  start = %d, " \
		  "stop = %d,  ", n, (long long)*refpoint, nintg_proc, 
		  start, stop);
        }
        start = stop;

            /* update refpoint  */
            /* 8 is the length of vblock string, 48 is length of 
	     * visblock struct    */
        *refpoint += nintg_proc * nbase * (bufsize + 8 + 
		       sizeof(struct visblock)-sizeof(struct visfreq));
        if(debug) printf("myid %d Updated refpoint = %lld\n",
			 myid, (long long)*refpoint);
        }

    return (0);
    }

