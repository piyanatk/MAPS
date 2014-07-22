#ifndef IONO_GEN_H
#define IONO_GEN_H

#define R_EARTH 6378.0 /* Radius of Earth, km */
void generate_sine(int, int, int, int, double *, double delta_lat, double start_lat);
void generate_iono_ramp(int, int, int, int, double, double, double, double, double, double, double, double, double *);
void generate_iono_parabolic(int, int, int, int, double, double, double, double, double, double, double, double, double *);
void generate_iono_turbulent(int t_dim, int x_dim, int y_dim, int z_dim, double c, double o_scale, double i_scale, const char* profile, double delta_lat, double delta_lon, 
                             double start_alt, double delta_alt, double delta_time, double* iono_grid);
void generate_iono_turbulent2D(int t_dim, int x_dim, int y_dim, int z_dim, double c, double o_scale, double i_scale,
                             double delta_lat, double delta_lon,double start_alt, double delta_alt, double delta_time, double* iono_grid);
float ran0(int *);

#endif
