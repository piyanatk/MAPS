/*************************************************************************
 *                                                                       *
 * convex_hull_2d.c                                                      *
 *                                                                       *
 * This function finds the convex hull, i.e. a minimal polygon for a set *
 * of N points on a plane. The points are provided as Cartesian pairs    *
 * (x,y) in the 2D array r[N,2]. The function returns a list of the      *
 * border indices.                                                       *
 *                                                                       *
 * Algorithm design: L. Benkevitch.                                      *
 * Created 13 December 2010 by L. Benkevitch.                            *
 *                                                                       *
 *************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "list_proc.h"
#include "convex_hull_2d.h"


struct list_int *convex_hull_2d(/* N, r */
				int N, double r[N][2]) 
{
  char ix[N];
  int imin, imax, i, j, extr, nbd, nbd1, jmaxr;
  struct list_int *pbd, *bd, *pbd1, *bd1;
  double ymin, ymax, ry, lseg, lpt, psinph, dmaxr;
  double seg[2], pt[2];

  /*
   * Each element of the array ix[N] contains current properties
   * of the points with coordinates r[N,2]:
   * -1: candidate for the interior;
   *  2: border point
   *  0: interior point
   *  1: exterior point
   * Obviously, for each ix element two bits (i.e. 4 states) is enough.
   */

  /*
   * Seed all the points with (-1)s 
   */
  for (i = 0; i < N; i++) ix[i] = -1;

  /*
   * Find 2 initial points with minimal and maximal y coordinates.
   *  They will be the initial border.
   */
  ymin = ymax = r[0][1];
  imin = imax = 0;
  for (i = 0; i < N; i++) {
    ry = r[i][1];
    if (ymin > ry) { ymin = ry; imin = i; }
    if (ymax < ry) { ymax = ry; imax = i; }
  }
  ix[imin] = ix[imax] = 2; /* Mark top and bottom points as border elements */

  /* printf("ymin = %g,  ymax = %g\n", ymin, ymax); */
  /* printf("imin = %d,  imax = %d\n", imin, imax); */

  /*
   * Form initial border vertex list from the two extremal points.
   * The border now consists of two (coincident) segments
   */
  /* bd = (struct list_int *) malloc(sizeof(struct list_int)); */
  /* bd->next = (struct list_int *) malloc(sizeof(struct list_int)); */
  /* bd->next->next = bd; /\* Close the border list into the ring *\/ */
  /* bd->num = imin; bd->next->num = imax; */
  bd = list_create_atom(imin); 
  list_insert(bd, imax); /* Insert imax AFTER the *bd atom */
  /* list_print("bd = ", bd); */
  nbd = 2;      /* Current border list length */
  nbd1 = nbd;   /* New border list length */

  /*
   * Main loop
   */
  extr = TRUE;
  /* int ncyc = 0 */
  while (extr) {
    pbd = bd; /* Set pointer at beginning of old border bd */
    bd1 = list_copy(bd);
    pbd1 = bd1; /* Set pointer at beginning of new border bd1 */
    for (i = 0; i < nbd; i++) {
      seg[0] = r[pbd->next->num][0] - r[pbd->num][0]; /* Segment vector, x */
      seg[1] = r[pbd->next->num][1] - r[pbd->num][1]; /* Segment vector, y */
      lseg = sqrt(pow(seg[0],2) + pow(seg[1],2));   /* Its length  */
      dmaxr = -1.0;  /* Right distance can only be positive */
      for (j = 0; j < N; j++) {
	if (ix[j] == 0) continue; /* Avoid interior (ix=0) */ 
	if (ix[j] == 2) continue; /* Avoid border (ix=2)   */ 
	pt[0] = r[j][0] - r[pbd->num][0]; /* Vecx from i-th vert. to j-th pt. */
	pt[1] = r[j][1] - r[pbd->num][1]; /* Vecy from i-th vert. to j-th pt. */
	lpt = sqrt(pow(pt[0],2) + pow(pt[1],2));  /* Its length */
	/* Left(-) or right(+) position of j-th point wrt to i-th segment */
	psinph = (seg[1]*pt[0] - seg[0]*pt[1])/lseg;  /* lpt*sin(phi) */
	if (psinph >= 0.0) {
	  ix[j] = 1; /* Mark r[j,:] as exterior point */
	  if (psinph > dmaxr) {  /* Ext. point farthest from seg. */
	    dmaxr = psinph;
	    jmaxr = j;
	  }
	}
      }    /* for j */
      /* printf("i= %d, jmaxr= %d, dmaxr= %g\n", i, jmaxr, dmaxr); */
      if (dmaxr >= 0) { /* Farthest exterior point found */
	list_insert(pbd1,jmaxr); /* Insrt new vert. bw i-th and (i+1)-th */
	nbd1++; /* Count one more element in new border */ 
	pbd1 = pbd1->next;
	ix[jmaxr] = 2;  /* Mark as border point */
      }   /* if (dmaxr >= 0) */
      pbd = pbd->next;
      pbd1 = pbd1->next;
    }   /* for i */
    extr = FALSE;
    for (j = 0; j < N; j++) {
      if (ix[j] == -1) ix[j] = 0;
      if (ix[j] ==  1) {
	extr = TRUE;
	ix[j] = -1; /* The points outside the current polygon */
      }
    }  /* for j */
    
    /* list_print("bd = ", bd); */
    /* list_print("bd1 = ", bd1); */
    /* printf("\n"); */

    bd = list_copy(bd1);
    list_del(bd1);
    nbd = nbd1;
  }   /* while (extr) */

  return bd;
}






