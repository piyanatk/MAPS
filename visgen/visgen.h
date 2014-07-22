#define VISGEN_SIZE_NAME 200

struct vg_param
{
  int   n_pol_products;
  int   n_pol_inputs;
  char  simname[VISGEN_SIZE_NAME];
  char  obsfilename[VISGEN_SIZE_NAME];
  char *gridfile[4];
  char  arrayfilename[VISGEN_SIZE_NAME];
  int   do_oob;
  char  oobfilename[VISGEN_SIZE_NAME];
  int   do_ionomodel;
  char  ionomodel_filename[VISGEN_SIZE_NAME];
  char  ionosettings_filename[VISGEN_SIZE_NAME];
  int   do_rfimodel;
  char  rfimodel_filename[VISGEN_SIZE_NAME];
  int   do_noise;
  int   do_vis;
  char  site[VISGEN_SIZE_NAME];
  double reftime;
  int   center_ion;
  float global_sefd; /* global station SEFD from command line */
  char  *module_debug_flags;
  char  layoutdirname[512];
  int   vlbi; /* flag; [-1,0,1]: option [-l,none,-v] * 
	       *added by L. Benkevitch, 09-Dec-2010  */
};
