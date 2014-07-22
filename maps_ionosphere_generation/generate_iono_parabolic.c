#include <stdio.h>
#include <math.h>
#include "iono-gen.h"

/* A function to generate an ionosphere with a parabolic z distribution
 *
 * Divya 29 March 05
 */

static const double rad2km = 6378.19;

void generate_iono_parabolic(int t_dim, int x_dim, int y_dim, int z_dim, double start_time, double start_lon, double start_lat, double start_alt, double delta_time, double delta_lon, double delta_lat, double delta_alt, double *iono_grid)
{
  int i, j, k, l, index;
  double a = -1*5.09055e+8, b = -1*5.09055e+8, c = 10e+12;
  double c_x = x_dim/2, c_y = y_dim/2;

  for (i = 0; i < t_dim; i++)
    for (j = 0; j < x_dim; j++)
      for (k = 0; k < y_dim; k++)
        for (l = 0; l < z_dim; l++)
	{
          index = i*x_dim*y_dim*z_dim + j*y_dim*z_dim + k*z_dim + l;

	  /* Generating an ionosphere with a parabolic z profile */

          iono_grid[index] = a*(rad2km*delta_lon*(c_x-j))*(rad2km*delta_lon*(c_x-j)) + b*(rad2km*delta_lat*(c_y-k))*(rad2km*delta_lat*(c_y-k)) + c;

	  /*if ((i == 0) &&(j == 0) && (k == 0))
            fprintf(stdout, "%d %d %d %d %.4e %.7e\n", i, j, k, l, start_alt+l*delta_alt, iono_grid[index]);*/
	}
}
