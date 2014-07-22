#include <stdio.h>
#include <math.h>
#include "iono-gen.h"

/* A function to generate an ionosphere
 * with a ramp along a specified direction
 *
 * Divya 02 March 05
 */

#define rad2km 6378.0

double transform_distance(double trans_coord, double z)
{
  return trans_coord * (z+rad2km)/rad2km;
}

void generate_iono_ramp(int t_dim, int x_dim, int y_dim, int z_dim, double start_time, double start_lon,
                        double start_lat, double start_alt, double delta_time, double delta_lon, double delta_lat,
                        double delta_alt, double *iono_grid)
{
  int i, j, k, l, index;
  double m = 1e+10;    /* Slope of the ramp (e m^-3 km^-1), 1e+16 e m^-2= 1 TEC */
  double c = 1e+2;   /* A constant DC level on which the ramp is riding */
  double x0;          /* x0, y0 - A reference location where the ramp */
  double y0;          /* has a value equal to c */
  double theta = 0.0; /* The angle along which the ramp is generated (deg) */
  double deg2rad, cos_theta, sin_theta;

  
  deg2rad = atan(1.0)/45.0;
  cos_theta = cos(theta*deg2rad);
  sin_theta = sin(theta*deg2rad);
  y0 = start_lat;
  x0 = start_lon;
  x0 = y0 = 0.0;

  fprintf(stdout, "cos_theta = %g, sin_theta = %g\n", cos_theta, sin_theta);
  fprintf(stdout, "x0 = %g, y0 = %g\n", x0, y0);
  fprintf(stdout, "delta_lon = %g, delta_lat = %g\n", delta_lon, delta_lat);
  for (i = 0; i < t_dim; i++)
    for (j = 0; j < x_dim; j++)
      for (k = 0; k < y_dim; k++)
        for (l = 0; l < z_dim; l++)
          {
            index = i*x_dim*y_dim*z_dim + j*y_dim*z_dim + k*z_dim + l;

            /* Generating a ramp along arbitrary directions in x,y plane
             * (the latitude, longitude plane) */
            /*iono_grid[index] = m*(rad2km*(cos_theta*(transform_distance(j*delta_lon, start_alt + l*delta_alt)) + sin_theta*(transform_distance(k*delta_lat, start_alt + l*delta_alt)))) + c;*/
            iono_grid[index] = m*(rad2km*(cos_theta*(j*delta_lon) + sin_theta*(k*delta_lat))) + c;
            /*	  fprintf(stdout, "%d %.5e\n", index, iono_grid[index]); */
            
            /*if ((i == 0) && (k == 0)) 
              fprintf(stdout, "%d %d %d %d %.7e\n", i, j, k, l, iono_grid[index]);*/
          }
}


