#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include <string.h>
#include "netcdf.h"
#include "iono-gen.h"

/* Generate an ionospheric refraction model and write it to a netcdf data file
 * for use in visgen.
 * History:
 * Divya          Sept, 04   Created
   Sean Ting      June, 06.  Added different models
   Randall Wayth. March, 07. Integrate with MAPS subversion environment.
                             remove hard-codeded stuff. Add command line args and
			     parameter checking, usage message.
 */

#define OUTER_SCALE 5.0
#define INNER_SCALE 25.0

typedef struct _ionogen_options {
    int sim_type;
    int t_dim,x_dim,y_dim,z_dim;
    double inner_scale,outer_scale;
    double width_km, height_km;
    double start_time, start_lat, start_lon, start_alt;
    double delta_time, delta_lat, delta_lon, delta_alt;
    char   *outfilename;
} ionogen_options_t;


/* function prototypes */
void define_iono_grid(ionogen_options_t *options, double **grid);
void usage(int argc, char *argv[],ionogen_options_t *options);
void parse_commandline(int argc, char *argv[], char *optstring, ionogen_options_t *options);
int write_NC_header(const char *outfile, int *outfile_id, int t_dim, int x_dim, int y_dim, int z_dim,
                    int *time_id, int *x_id, int *y_id, int *z_id, int *grid_id,
                    double start_time, double delta_time, double start_lat, double delta_lat,
                    double start_lon, double delta_lon, double start_alt, double delta_alt);
/* global variables */
int g_debug=0;

int main(int argc, char *argv[])
{
  ionogen_options_t options;
  int    status, outfile_id;
  int    time_id, x_id, y_id, z_id, grid_id;
  double *iono_grid = NULL;
  size_t start[4] = {0, 0, 0, 0}, count[4];
  char   optstring[]="do:t:x:z:w:h:I:O:";

    /* set defaults */
  options.width_km  = 100.0;
  options.height_km = 100.0;
  options.sim_type = 'R';
  options.t_dim = 4;
  options.x_dim = options.y_dim = 128;
  options.z_dim = 4;
  options.outfilename = NULL;
  options.start_time = 0.0;
  options.delta_time = 25.0; /* sec */
  options.start_lon  = 2.05;
  options.start_lat  = -.461;
  options.inner_scale = INNER_SCALE;
  options.outer_scale = OUTER_SCALE;
  /* for 400km thick ionosphere over z_dim-1 pixels */
  options.start_alt  = 100.0;  /* starting height (km) */
  options.height_km  = 400.0; /* thickness of ionosphere (km) */

  /* print a usage message if there are no args */
  if (argc < 2) usage(argc,argv,&options);
  parse_commandline(argc, argv, optstring,&options);

  if(options.x_dim < 2 || options.y_dim < 2 || options.z_dim < 2 || options.x_dim%4 ||
            options.y_dim%4 || options.z_dim%4) {
    fprintf(stderr,"Bad dimensions: x,y,z,t (%d,%d,%d,%d). Dimensions must have at least 2 pixels and be power of 2.\n",
            options.x_dim, options.y_dim, options.z_dim,options.t_dim);
    exit(1);
  }

  /* initialise */
  define_iono_grid(&options,&iono_grid);

  /* Write the header information for the ionosphere description */  
  status = write_NC_header(options.outfilename, &outfile_id, options.t_dim, options.x_dim, options.y_dim, options.z_dim,
            &time_id, &x_id, &y_id, &z_id, &grid_id, options.start_time, options.delta_time, options.start_lon, options.delta_lon,
            options.start_lat, options.delta_lat, options.start_alt, options.delta_alt);
  if (status != NC_NOERR) {
    fprintf(stderr, "### Error define_NC, ret = %d\n", status);
    return -1;
  }
  if (g_debug) fprintf(stdout, "# Planted header information for NC data file %s (ID %d)\n", options.outfilename, outfile_id);


/* Generate the ionosphere */
  switch (options.sim_type) {
  case 'S' :
    generate_sine(options.t_dim, options.x_dim, options.y_dim, options.z_dim, iono_grid, options.start_lat,
                options.delta_lat);
    break;
  case 'T':
    generate_iono_turbulent(options.t_dim, options.x_dim, options.y_dim, options.z_dim, .0045, OUTER_SCALE, INNER_SCALE,
         "./Profiles/HeightProfile2.txt", options.delta_lat, options.delta_lon, options.start_alt, options.delta_alt,
            options.delta_time, iono_grid);
    break;
  case '2':
    generate_iono_turbulent2D(options.t_dim, options.x_dim, options.y_dim, options.z_dim, .0045, OUTER_SCALE, INNER_SCALE,
         options.delta_lat, options.delta_lon, options.start_alt, options.delta_alt, options.delta_time, iono_grid);
    break;
  case 'R':
    generate_iono_ramp(options.t_dim, options.x_dim, options.y_dim, options.z_dim, options.start_time, options.start_lon,
                 options.start_lat, options.start_alt, options.delta_time, options.delta_lon, options.delta_lat,
                 options.delta_alt, iono_grid);
    break;
  case 'P':
    generate_iono_parabolic(options.t_dim, options.x_dim, options.y_dim, options.z_dim, options.start_time,
                options.start_lon, options.start_lat, options.start_alt, options.delta_time, options.delta_lon,
                 options.delta_lat, options.delta_alt, iono_grid);
    break;
  default:
    fprintf(stderr,"ERROR: unknown ionosphere type %c\n",options.sim_type);
    usage(argc,argv,&options);
  }

/* Populating the NC dataset */
  count[0] = options.t_dim; count[1] = options.x_dim; count[2] = options.y_dim; count[3] = options.z_dim;
  status = nc_put_vara_double(outfile_id, grid_id, start, count, iono_grid);
  if (status != NC_NOERR) 
  {
    fprintf(stderr, "### Error in nc_put_vara_double, ret = %d\n", status);
    return -1;
  }
  if (g_debug) fprintf(stdout, "# Populated the NC data structure.\n");

/* Close the NC file to flush the buffer */
  status = nc_close(outfile_id);
  if (status != NC_NOERR) 
  {
    fprintf(stderr, "### Error in nc_close, ret = %d\n", status);
    return -1;
  }
  if (g_debug) fprintf(stdout, "# NC data file closed.\n");

  return 0;
}


void define_iono_grid(ionogen_options_t *options, double **iono_grid)
{	
/* Filling in numbers for array dimensions */
  /* un hard-coded. RBW Jan, 2008 */

/* Filling in numbers to define the origin and increments along each axis */
/* Roger's code treats longitude and latitude identically as increments
 * in latitude (degrees) along a great circle, hence delta_lat = delta_lon 
 * Roger uses deg2km = 111.320 
 * I think the numbers for start_lat and start_lon are ignored. Need to
 * check this out with Roger. */


  /* for 100km span at Earth's radius over x_dim pixels, need */
  options->delta_lon = (options->width_km/(options->x_dim-1))/R_EARTH; /* rad */
  options->delta_lat = (options->width_km/(options->y_dim-1))/R_EARTH; /* rad */
  options->delta_alt = options->height_km/(options->z_dim-1);
  
  /* Memory allocation */
  *iono_grid = (double *)calloc(options->t_dim*options->x_dim*options->y_dim*options->z_dim, sizeof(double));
  
  printf("Grid params\n");
  printf("Lon:  dim: %d,\tstart: %g,\tsize: %g (rad),\tdelta: %g\n",options->x_dim,options->start_lon,
            (options->delta_lon)*(options->x_dim-1),options->delta_lon);
  printf("Lat:  dim: %d,\tstart: %g,\tsize: %g (rad),\tdelta: %g\n",options->y_dim,options->start_lat,
            (options->delta_lat)*(options->y_dim-1),options->delta_lat);
  printf("Alt:  dim: %d,\tstart: %g,\tsize: %g (km),\tdelta: %g\n",options->z_dim,options->start_alt,
            options->delta_alt*(options->z_dim-1),options->delta_alt);
  printf("Time: dim: %d,\tstart: %g,\tsize: %g (s),\tdelta: %g\n",options->t_dim,options->start_time,
            options->delta_time*(options->t_dim-1),options->delta_time);
}

/* A function to write out a NetCDF header file describing
 * an ionosphere in the same format as Gary Bust's Python 
 * code. 
 *
 * Expects the information to be written out in the header
 * to be supplied via the argument list.
 *
 * Divya, 16 Sept, 04 (created)
*/
int write_NC_header(const char *outfile, int *outfile_id, int t_dim, int x_dim, int y_dim, int z_dim,
                    int *time_id, int *x_id, int *y_id, int *z_id, int *grid_id,
                    double start_time, double delta_time, double start_lat, double delta_lat,
                    double start_lon, double delta_lon, double start_alt, double delta_alt)
{
  int status;
  int grid_dim_id[4];
  char *title = "Test ionosphere";
  char *units = "electrons per meter cubed";

  /* create netCDF dataset: enter define mode */
  status = nc_create(outfile, NC_CLOBBER, outfile_id);
  if (status != NC_NOERR) 
  {
    fprintf(stderr, "### Error in nc_create, ret = %d\n", status);
    return status;
  }

/* define dimensions: from name and length */
  status = nc_def_dim(*outfile_id, "t", t_dim, time_id) 
         | nc_def_dim(*outfile_id, "x", x_dim, x_id)
         | nc_def_dim(*outfile_id, "y", y_dim, y_id)
         | nc_def_dim(*outfile_id, "z", z_dim, z_id);
  if (status != NC_NOERR) 
  {
    fprintf(stderr, "### Error in nc_def_dim, ret = %d\n", status);
    return -1;
  }

  grid_dim_id[0] = *time_id;
  grid_dim_id[1] = *x_id;
  grid_dim_id[2] = *y_id;
  grid_dim_id[3] = *z_id;
  status = nc_def_var(*outfile_id, "ionosphere_grid", NC_FLOAT, 4, grid_dim_id, grid_id)
         | nc_put_att_text(*outfile_id, *grid_id, "units", strlen(units), units);
  if (status != NC_NOERR) 
  {
    fprintf(stderr, "### Error in nc_def_var/nc_put_att_text, ret = %d\n", status);
    return -1;
  }

/* put attribute: assign attribute values */
  status = nc_put_att_text(*outfile_id, NC_GLOBAL, "title", strlen(title), title) \
         | nc_put_att_double(*outfile_id, NC_GLOBAL, "start_time", NC_DOUBLE, 1, &start_time)
         | nc_put_att_double(*outfile_id, NC_GLOBAL, "start_lat",  NC_DOUBLE, 1, &start_lat)  
         | nc_put_att_double(*outfile_id, NC_GLOBAL, "start_lon",  NC_DOUBLE, 1, &start_lon)  
         | nc_put_att_double(*outfile_id, NC_GLOBAL, "start_alt",  NC_DOUBLE, 1, &start_alt) 
         | nc_put_att_double(*outfile_id, NC_GLOBAL, "delta_time", NC_DOUBLE, 1, &delta_time) 
         | nc_put_att_double(*outfile_id, NC_GLOBAL, "delta_lat",  NC_DOUBLE, 1, &delta_lat)  
         | nc_put_att_double(*outfile_id, NC_GLOBAL, "delta_lon",  NC_DOUBLE, 1, &delta_lon) 
         | nc_put_att_double(*outfile_id, NC_GLOBAL, "delta_alt",  NC_DOUBLE, 1, &delta_alt);
  if (status != NC_NOERR) 
  {
    fprintf(stderr, "### Error in defining attributes, ret = %d\n", status);
    return -1;
  }

/* end definitions: leave define mode */
  status = nc_enddef(*outfile_id);
  if (status != NC_NOERR) 
  {
    fprintf(stderr, "### Error in nc_enddef, ret = %d\n", status);
    return -1;
  }
  
  return status;
}

void parse_commandline(int argc, char *argv[], char *optstring, ionogen_options_t *options) {
    int opt_char;
  /* parse command line args */
  while((opt_char = getopt(argc,argv,optstring)) != -1) {
    switch (opt_char) {
    case 'd': g_debug=1;
      break;
    case 'o': options->outfilename = optarg;
      break;
    case 't': options->sim_type = optarg[0];
      break;
    case 'x': options->x_dim = atoi(optarg);
      break;
    case 'z': options->z_dim = atoi(optarg);
      break;
    case 'w': options->width_km = atof(optarg);
      break;
    case 'h': options->height_km = atof(optarg);
      break;
    case 'I': options->inner_scale = atof(optarg);
      break;
    case 'O': options->outer_scale = atof(optarg);
      break;
    default:
      fprintf(stderr,"unknown option '%c'\n",opt_char);
      exit(1);
    }
  }
}

/******************************
usage 
****************************/
void usage(int argc, char *argv[],ionogen_options_t *options) {
  fprintf(stderr,"Usage:\n%s [options]\n",argv[0]);
  fprintf(stderr,"\t-o <output_filename> \n");
  fprintf(stderr,"\t-t <R,T,S,P,2> R=ramp (default), T=turbulent, P=parabolic, S=sine, 2=2D turb\n");
  fprintf(stderr,"\t-x dimension of lat & long. Must be power of 2. Default %d \n",options->x_dim);
  fprintf(stderr,"\t-z z (height) dimension. Must be power of 2. Default %d \n",options->z_dim);
  fprintf(stderr,"\t-w width (km) of lat & long. Default %g \n",options->width_km);
  fprintf(stderr,"\t-h height (km) (top-bottom) of ionosphere. Default %g \n",options->height_km);
  fprintf(stderr,"\t-I inner scale length (km) of turbulence. Default %g \n",options->inner_scale);
  fprintf(stderr,"\t-O outer scale length (km) of turbulence. Default %g \n",options->outer_scale);
  fprintf(stderr,"\t-d turn debugging output on \n");
  exit(0);
}

