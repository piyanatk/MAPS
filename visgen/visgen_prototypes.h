#ifndef VISGEN_PROTOTYPES_H
#define VISGEN_PROTOTYPES_H

// Collect all function prototypes for visgen main program into one place

int compute_gha (struct observing_param *obsparam, double time, double *gha);
int compute_uvw (struct stationspec *st1, struct stationspec *st2, double gha,
		 struct observing_param *obsparam, double *u_len, 
		 double *v_len, double *w_len);
int get_site (char *sitename, struct array_specification *arrayspec);
int integrate (struct uvpatch *patch, struct stationspec *stn1, 
	       struct stationspec *stn2, double time, double freq, 
	       int scan_no, int intg_no, struct observing_param *pobs, 
	       complex *visibility, double *umean, double *vmean, int n_pols);
int ion_density (struct ion_model *isphere, double x[3], double *value);
/* file_offset() modified by L. Benkevitch.                          * 
 * the 'FILE *outfp' argument omitted.                               *
 * int file_offset (int numprocs, int myid, off_t *refpoint,         *
 *                  struct observing_param *obsparam,                *
 *		    int scan, int nst, FILE *outfp, int *intg_start, *
 *                  int *intg_num, int npol, off_t *write_pos);      */
int file_offset (int numprocs, int myid, off_t *refpoint, 
		 struct observing_param *obsparam, int scan, int nst, 
		 int *intg_start, int *intg_num, int npol, off_t *write_pos);
int fill_beamgrid (struct beamgrid *beamgrid, double cellsize,
		   struct observing_param *obsparam);
int interp_t (struct ion_model *isphere,double t);
int open_gridfile (char *filename[4], double fov, struct uvgrid_parms *uvgrid);
/* parse_cmdline () modified by J. Steeger and L. Benkevitch         *
 * Argument 'FILE **asciifp' added for visibility ascii output       */
int parse_cmdline (int argc, char **argv, struct vg_param *vgparam, 
		   FILE **outfp, FILE **asciifp);
int patch_size (struct stationspec *st1, struct stationspec *st2, double start, 
		double freq, struct observing_param *obsparam, int scan_index, 
		int integration_numner, double cellsize,
		struct conv_fn *beamconvl,  double *udim, double *vdim);
int read_array_spec (char *filename, char * dirname, 
		     struct array_specification *arrayspec);
int read_obs_spec (char *filename, struct observing_param *obsparam,
		   struct vg_param *vgparam, double array_lon_radian);
int read_rfimodel (char *filename, struct rfi_model *rfimodel);
int write_main_header (struct observing_param *obsparam, 
		       struct vg_param *vgparam,
		       struct array_specification *arrayspec, FILE *outfp);
int write_stnfile (struct array_specification *arrayspec, char *simname);
int write_timeblock_header (struct observing_param *obsparam, int scan_no,
			    struct vg_param *vgparam, FILE *outfp);

int check_baseline(struct vg_param vgparam);
#endif
