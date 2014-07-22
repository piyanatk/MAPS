#include <stdio.h>
#include <math.h>
#include "iono-gen.h"

/* A function to generate an ionosphere.
 *
 * Divya 16 Sept. 04
 */
extern int g_debug;

void generate_sine(int t_dim, int x_dim, int y_dim, int z_dim, double *iono_grid, double delta_lat, double start_lat)
{
  int i, j, k, l, index;
  double lambda=2.0;
  double scale = 1e+15; /* 1e+15 = 0.1 TEC */

  if (g_debug) {
    printf("Generating sine wave with lamda: %g, scale: %g, delta_lat: %g\n",lambda,scale,delta_lat);
  }
  
  for (i = 0; i < t_dim; i++)
    for (j = 0; j < x_dim; j++)
      for (k = 0; k < y_dim; k++)
        for (l = 0; l < z_dim; l++) {
          index = i*x_dim*y_dim*z_dim + j*y_dim*z_dim + k*z_dim + l;
          /* Generating a ramp along the latitude axis 
             iono_grid[index] = (j)*scale; */
          
          /* A sine wave along the latitude axis 
           * amplitude = scale; wavelength = */
          iono_grid[index] = scale*sin(2.0*M_PI*(j*delta_lat+start_lat)/lambda);
          /*
            if ((i == 0) &&(j == 50) && (k == 0)) fprintf(stdout, "%.3f\n", iono_grid[index]);
          */
        }
}

