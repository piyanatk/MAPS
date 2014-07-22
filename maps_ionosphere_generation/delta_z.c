#include <stdio.h>
#include <math.h>

int main(void)
{
  double z, ne, delta_z, delta_h, h_i, nu_p, nu, r_0;
  double sz, cz, a, b, c, d, deg2rad, z_incr;
  double z_rad, rad2deg;

  rad2deg = 180.0/(4*atan(1));
  deg2rad = 1/rad2deg;
  h_i = 350000.0; 	/* 350 km */
  ne  = 1e+12;          /* 1e+12 el m^-3 */
  delta_h = 500000.0; 	/* 500 km */
  nu_p = 8.9787*sqrt(ne); 
  nu = 100.0E+6; 	/* 100 MHz */
  r_0 = 6378166.0;	/* Radius of Earth in m */
  z_incr = 0.1;

  b = (nu_p/nu);
  b *= b;
  b *= (1 + h_i/r_0);
  d = 2.0*h_i/r_0;
  z = 0.0;
  while (z < 45.0)
  { 
    z_rad = z*deg2rad;
    sz = sin(z_rad);
    cz = cos(z_rad);
    a = 2.0*delta_h*sz/(3.0*r_0);
    c = cz*cz + d;
    c = pow(c, -1.5);
    delta_z = a*b*c;
    fprintf(stdout, "%.3f %.5e\n", z, delta_z*rad2deg);
    z += z_incr;
  }

  return 0;
}
