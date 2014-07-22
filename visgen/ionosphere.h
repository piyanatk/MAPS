/************************************************************************
    public structure definitions and function prototypes for ionosphere.c
   
    Jan, 2008. Randall Wayth.
    Original contributors:
    Roger Cappallo, Colin Lonsdale, Shep Doeleman

************************************************************************/

#ifndef IONOSPHERE_H
#define IONOSPHERE_H

#include <fitsio.h>


#define IONO_MAX_FNAME 256
#define deg_to_km  111.317         // number of km in 1 deg of latitude at r =6378 km
#define rad_to_km  6378.0          //   "    "  "   " " rad  "    "      "    "
enum coordinates {T, X, Y, Z};

/* public structure definitions */
struct iono_screen {
  int     num_stations;    /* number of stations to compute a phase screen for */
  int     size_x,size_y;   /* Size of function in both x, y. */
  float   **phase_delay;   /* Radians at 1 MHz. Array size of num_stations, pointing to 2D
			      arrays size_x*size_y pixels. These pixel sky locations are set by
			      the colvolution beam grid, or the location of OOB sources. */
  int     *valid;          /* array the size of num_staions. Boolean flag to indicate if the
			      phase screen calc is valid for the station */
};


/* Define turbulence cube structure */
struct tcube_str {
  int nx, ny, nz;            /* Number of elements in each direction */
  double sx, sy, sz;          /* Step size in each direction [km] */
  double vx, vy, vz;          /* Translation velocity in each direction [km/s] */
  float  ***cube;             /* The 3D cube of density perturbations */
  double i_scale, o_scale;    /* Inner and outer scale of turbulence */
  double rms;                 /* rms amplitude of the turbulence */
};

/* Define a path */
struct path_str {
  
  double x0, y0, z0; /* starting position in km relative to 
                     the center of the array */
  double nx, ny, nz; /* unit vector for step */
  double step;    /* length of each step in km */
  double steps;      /* total number of steps to take */
  
};


struct ion_model {
  char iono_fname[IONO_MAX_FNAME];     // name of input ionsphere model file
  char settings_fname[IONO_MAX_FNAME];     // name of ionsphere settings file
  fitsfile *fpfit;                // pointer to open fits file of global ionosphere
  int file_id;                    // file identifier used by netcdf
  int grid_id;                    // identifies the current dataset for netcdf
  int  global_type;               // IRI in fits=1 or netcdf=0 input?
  int  zlog;                      // Z sampling is (0=linear, 1=log)
  int  do_turbulence;             // Do=1 or do not=0 include turbulence
  int  do_tid;                    // Do=1 or do not=0 include travelling ionospheric discontinuity
  int  do_simple;                 // Do=1 or do not=0 do simple linear density gradient
  int do_force;                   // Do=1 or do not=0 force time in UT
  double force_ut;                // Forced time in UT seconds
  double tid_amp;                 // Amplitude of TID as fraction at 300km
  double tid_amin;                // TID minimum altitude [km]
  double tid_amax;                // TID maximum altitude [km]
  double tid_vel;                 // TID longitudinal drift speed [km/s]
  double tid_len;                 // TID length scale [km]
  double simple_grad;             // Electron gradient e/cm^3/km
  double simple_dc;               // Electron DC value e/cm^3
  double simple_zmin;             // Zmin for simple [km]
  double simple_zmax;             // Zmax for simple [km]

  size_t dim[4];                  // dimension for each axis (nx, ny, nz, nt)
  double origin[4];               // ind. var. values for 1st element of grid
  double delta[4];                //  "    "   element spacings for all 4D
  int ntslice;                    // index of first of 4 time slices currently in RAM
  float  *tslice[4];              // ptrs to 4 sequential 3D time slices in RAM
  double t;                       // time to which 4D->3D interpolation is done
  float  *slice;                  // field of 3D values interpolated to time t
  double center_lat;              // Center latitude of global model in radians
  double center_lon;              // Center longitude of global model in radians
  struct tcube_str *tcube;
};


/* public function prototypes */
int path (struct ion_model *isphere, double enh[3], double az, double el,double *integral);
int open_ionosphere (struct ion_model *isphere);
int init_ionoscreen(const int num_stations, 
		    const int gridsize_x, const int gridsize_y, struct iono_screen *ionoscreen);
int compute_ionoscreen(struct stationspec *st,
		       int                 st_index,
		       struct beamgrid    *bgrid,
		       double             gha,
		       double             fovra,
		       double             padfact,
		       struct ion_model   *isphere,
		       int                 center_ion,
		       struct iono_screen *ionoscreen);
void ionosphere_set_debug(int level);
#endif

