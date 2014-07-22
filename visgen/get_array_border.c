#include <stdio.h>
#include <stdlib.h>
#include "array_specification.h"
#include "list_proc.h"
#include "convex_hull_2d.h"
#include "get_array_border.h"

/*
 * Find the stations that are most remote from each other 
 * to only check the longest baselines.
 * Such stations form the "convex hull"  or border of the 
 * set of stations, scattered over the plane (we assume the 
 * earth is flat :)). The convex hull of a point set on
 * a plane comprises the verices of a minimal polygon
 * such that all the points other than the border points
 * are situated in the interior of the polygon. 
 * This may greatly reduce the computations for checking 
 * the baselines.
 * This procedure only makes sense for non-VLBI arrays.
 * The Cartesian coordinates of stations are:
 *   arrayspec.station[i].east_offset  -- x[i];
 *   arrayspec.station[i].north_offset -- y[i].
 *
 * The function returns the number of vertices in the border.
 * The indices of the border verices are rerurned in 
 * array_border array.
 *
 */
 
int get_array_border(/* *arrayspec, array_border[] */
		     struct array_specification *arrayspec, 
		     int array_border[]) 
{
  int i, n_border = 0;
  struct list_int *border, *ptr, *next;;
  int nst = arrayspec->nst;
  double coord[nst][2];
  FILE *fh;
  
  /*
   * Copy array (x,y) coordinates from arrayspec->(east_offset,north_offset)
   * to coord[nst][2]. 
   */

  fh = fopen("points.txt", "w");
  for (i = 0; i < nst; i++) {
    coord[i][0] = arrayspec->station[i].east_offset;  /* x station coordinate */
    coord[i][1] = arrayspec->station[i].north_offset; /* y station coordinate */
    fprintf(fh, "%12.4e %12.4e\n", coord[i][0], coord[i][1]);
  }
  fclose(fh);

  /*
   * Find the convex hull of the array stations -- the minimal poligon 
   * vertices. The vertex indices are returned in the list 'border', which
   * is closed into ring (the last element points at the first).
   */

  border = convex_hull_2d(nst, coord);

  /* list_print("border indices: ", border); */

  fh = fopen("bd_idx.txt", "w");
  ptr = border;
  do {
    i = ptr->num;
    ptr = ptr->next;
    fprintf(fh, "%d ", i);
  } while (ptr != border);
  fprintf(fh, "\n");
  fclose(fh);


  /*
   * Fill in array_border[:] with the indices from the list border
   */
  ptr = border;
  next = NULL;
  n_border = 0;
  while(next != border) {
    if (i == MAX_ARRAY_BORDER) {
      msg("ERROR in get_array_border(): the number of border elements " \
	  "exceeds  MAX_ARRAY_BORDER = %d.", 3, MAX_ARRAY_BORDER);
      msg("Increase MAX_ARRAY_BORDER constant in file array_specifications.h",
	  3);
      exit(1);
    }
    array_border[n_border++] = ptr->num;
    next = ptr->next;
    ptr = next;
  }

  list_del(border);

  return n_border;
}


//MAX_ARRAY_BORDER
