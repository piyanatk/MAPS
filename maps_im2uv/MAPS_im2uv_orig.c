/**********************************
MAPS_im2uv: a program to read a FITS sky brightness image file,
            reorder the data, fourier transform it and write it
            out in the binary format expected by VISGEN.
            If an appropriate source file is given, a point source
            is added to the uv grid for each source in the file.

author: Randall Wayth. Oct, 2006
Modiications (Leonid Benkevitch):
2010-May-15    Added east-west flipping of brightness image
2010-Jun-01    Added "long" cmd option style (L. Benkevitch)
2010-Jun-02    Provided ability to read multidimensional images (L. Benkevitch)
2011-Aug-25    Added saving central cut (crop) of the result (L. Benkevitch)
***********************************/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <fftw3.h>
#include <getopt.h>
#include <fitsio.h>
#include <string.h>

/* function prototypes */
//static void usage(char argv[]);
static void im2uv_help();

/* global vars */
int debug=1;
//static int verbose_flag;

int main(int argc, char *argv[]) {

  fitsfile *infile=NULL;
  FILE *outfile=NULL, *sourcefile=NULL;
  //  FILE *tempfile;
  char *infilename=NULL;
  char optstring[] = "a:i:o:s:f:h:l:n:dhtqmp:c:";
  char *outfilename=NULL, *sourcefilename=NULL;
  char line[256], str[64];
  char units[64];               // Added by benkev
  int result=0, naxis=0, anynull=0, i, j, n_threads=1;
  int ind=0, iind=0, jind=0, ipix=0, jpix=0, status=0;
  int  src=0, num_src=0, nkeys=0;
  int imgsize[2];
  long naxes[16], arr_size=0;
  double *data=NULL, *mirror=NULL, normalizer=1.0, ptsrc_normalizer=1.0;
  double  temp_r, phase;
  float za, *az=NULL, *el=NULL, *src_amp=NULL;
  float  param1, param2, fd, freq=0.0, ref_freq, spec_index;
  float lst=0.0, lat=0.0, ra=0.0, dec=0.0, ra0=0.0, dec0=0.0;
  float dl=0.0, dm=0.0, du=0.0, dv=0.0, u=0.0, v=0.0;
  float *l=NULL, *m=NULL;
  fftw_plan theplan;
  time_t thetime;
  // added by benkev
  double *pbeg, *pmid, *pend;
  long halfaxis0, pincr;
  short normalizer_specified = 0; // Flag indicating if -n is in cmd line
  short ewtranspose_specified = 0; // Flag indicating if -t is in cmd line
  double const norm_Jy_rad2 = 0.0304617;  // For flux in Jy/steradian
  char *strJy_rad2[] = {"Jy/steradian", "Jy/sr", "Jy steradian^-1", 
			"Jy steradian^(-1)", "Jy/rad^2", "Jy/rad2", 
			"Jy rad^-2", "Jy rad^(-2)"};
  int optc;                              // Command line option character
  int nxdimfromcmd = 0; // Number of dimensions from --axis option
  int nextradim = 0; // Number of fitsfile dimensions above two for image
  char *axestr = NULL; // Subscript string from cmd line
  long cropx = 0, cropy = 0; // Crop dimensions; 0 - no cropping
  char *cropstr = NULL; // Crop dimensions string from cmd line
  int iax[16] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1}; // First slice by default 
  long imgstart = 1;  // Where in the fits file data the image plane located
  char *brkt;
  short queryfitsheader = 0; // Print the dimensions of input fitsfile
  long padzpix = 0; /* Number of pixels from all sides to pad br. image */
  char ctype[16], ctype_num[16] ;
  double *datx = NULL; /* Temporary data container before zero-padding */
  long n_cols, n_rows; /* Dimensions of image plane */ 
  long icol1, icol2, irow1, irow2; // Cropping limits

  normalizer_specified = 0;   // Assume no normalizer (-n) in cmd line 
  ewtranspose_specified = 0;  // Assume no east-west flipping (-t) in cmd line
  if (argc == 1) im2uv_help();

  /* 
   * Parse command line arguments 
   */
  while (1) {
    static struct option long_options[] =
      {
	/* These options set a flag. */
	//{"verbose", no_argument,       &verbose_flag, 1},
	//{"brief",   no_argument,       &verbose_flag, 0},
	/* These options don't set a flag.
	   We distinguish them by their indices. */
	{"inputfilename",     required_argument, 0, 'i'},
	{"outputfilename",    required_argument, 0, 'o'},
	{"sourcefilename",    required_argument, 0, 's'},
	{"frequencyofsource", required_argument, 0, 'f'},
	{"lstofsource",       required_argument, 0, 'h'},
	{"latitude",          required_argument, 0, 'l'},
	{"normalizer",        required_argument, 0, 'n'},
	{"debug",             no_argument,       0, 'd'},
	{"help",              no_argument,       0, 'm'},
	{"axes",              required_argument, 0, 'a'},
	{"query",             no_argument,       0, 'q'},
	{"ewtranspose",       no_argument,       0, 't'},
	{"padzeropixels",     required_argument, 0, 'p'},
	{"crop",              required_argument, 0, 'c'},
	{0, 0, 0, 0}
      };
    /* getopt_long stores the option index here. */
    int option_index = 0;
     
    optc = getopt_long (argc, argv, optstring,	long_options, &option_index);
     
    /* Detect the end of the options. */
    if (optc == -1)
      break;

    switch (optc)
      {
      case 'i': infilename = optarg;
	break;
      case 'o': outfilename = optarg;
	break;
      case 's': sourcefilename = optarg;
	break;
      case 'f': freq = atof(optarg);
	break;
      case 'h': lst = atof(optarg);
	break;
      case 'l': lat = atof(optarg);
	break;
      case 'd': debug = 1;
	break;
      case 'q': queryfitsheader = 1;
	break;
      case 'a': axestr = optarg;
	break;
      case 'c': cropstr = optarg;
	break;
      case 'm': 
	im2uv_help();
	exit(0);
	break;
      case 't': ewtranspose_specified = 1; // added by benkev
	break;
      case 'p': padzpix = atoi(optarg) ; // added by benkev
	printf("padzpix = %ld\n", padzpix);
	break;
      case 'n': 
	normalizer = atof(optarg);
	normalizer_specified = 1;  // added by benkev
	break;
      default:
	abort ();
      }
  }



  if (infilename == NULL){
    fprintf(stderr, "FATAL: Name for input file must be supplied\n");
    exit(1);
  }
  if ((outfilename == NULL) && (!queryfitsheader)){
    fprintf(stderr, "FATAL: Names for output files must be supplied\n");
    exit(1);
  }

  if(debug) printf("infilename: %s\n", infilename);
  if(debug) printf("outfilename: %s\n", outfilename);
  if(debug) printf("sourcefilename: %s\n", sourcefilename);
  if(debug) printf("number of threads: %d\n", n_threads);

  /*
   * Process the --query option:
   * Print the data cube dimensions in the fitsfile 
   * as table:
   * CTYPE#: <name>, <dim>
   */
  if (queryfitsheader) {
    fits_open_image(&infile,infilename,READONLY, &result);
    if (result !=0) {
      fprintf(stderr, "FATAL: Failed to open file <%s> in FITSio. " \
	      "Error code: %d\n", infilename, result);
      exit(1);
    }
    /* get number of image dimensions */
    fits_get_img_dim(infile, &naxis, &result);
    if (result != 0) {
      fprintf(stderr,"FATAL: Failed to get image naxis with error code %d\n",
	      result);
      exit(1);
    }
    /*
     * Read up to 16 datacube dimensions into naxes[16] array
     */
    fits_get_img_size(infile, 16, naxes, &result);
    if (result != 0) {
      fprintf(stderr,"FATAL: Failed to get axis sizes with error code %d\n",
	      result);
      exit(1);
    }
    printf("Axes of the fits file data cube:\n");
    for (i = 0; i < naxis; i++) {
      sprintf(ctype_num, "CTYPE%d", i+1);
      fits_read_key(infile, TSTRING, ctype_num, ctype,  NULL, &status);
      if (status !=0) break;
      printf("%s: %s,   \t%ld\n", ctype_num, ctype, naxes[i]);
    }
    printf("\n");
    fits_close_file(infile,&result);
    if (result !=0) 
      fprintf(stderr,"WARNING: fits_close_file() failed with code %d. " \
	      "Do we care?\n", result);
    exit(0);
  }


  /*
   * Parse the axes subscript string
   */
  if (axestr != NULL) {
    char *tok = NULL;
    int j;
    tok = strtok_r(axestr, "[,]", &brkt); 
    nxdimfromcmd = 0;
    while (tok) {
      j = 2 + (nxdimfromcmd++);
      if (j > 15) {
	fprintf(stderr, "FATAL: Number of fitsfile datacube indices in "\
		"--axes option exceeds 14\n");
	exit(1);
      }
      iax[j] = atoi(tok); // Fill in axes indices,iax[], starting from the 2nd
      tok = strtok_r(NULL, "[,]", &brkt);
      
    }
  }
  printf("nxdimfromcmd = %d\n", nxdimfromcmd);
  for (i = 0; i < nxdimfromcmd; i++) {
    //printf("iax[%d] = %d\n", i, iax[2+i]);
    if (iax[2+i] <= 0) {
      printf("FATAL: wrong subscripts in --axes; they must be positive " \
	     "integers.\nExample: --axes [5,1,7]\n");
      exit(1);
    } 
  }


  /*
   * Parse the crop dimensions string
   */
  if (cropstr != NULL) {
    char *tok = NULL;
    cropx = 0;
    cropy = 0;
    tok = strtok_r(cropstr, "[,]", &brkt); 
    cropx = atoi(tok);
    if (*brkt != 0) {
      tok = strtok_r(NULL, "[,]", &brkt);
      cropy = atoi(tok);
    }
    else
      cropy = cropx;
    if (*brkt != 0) {
      fprintf(stderr, "FATAL: extra characters in '-c' option.\n");
      exit(14);
    }
    printf("cropx = %ld, cropy = %ld\n", cropx, cropy);
    if (cropx <= 0 || cropy <= 0) {
      printf("FATAL: wrong --crop dimensions; must be positive "	\
	     "integers.\nExample: -c1024   -c 980,1100   -c[512,256]\n");
      exit(1);
    }
    /* /\* Check if crop parities are the same as those of the source image *\/ */
    /* i = naxes[0]%2; */
    /* j = naxes[1]%2; */
    /* if (cropx%2 != i || cropy%2 != j) { */
    /*   printf("FATAL: wrong --crop dimensions; must be the same parity " \ */
    /* 	     "with the fits image dimensions.\n" \ */
    /* 	     "Query them using MAPS_im2uv.exe -i file.fits -q\n"); */
    /*   printf("i = %d, j = %d, cropx mod 2 = %ld, cropy mod 2 = %ld\n", i, j,  */
    /* 	     cropx%2, cropy%2); */
    /*   printf("cropx  = %ld, cropy  %ld\n", cropx, cropy); */
    /*   printf("naxes[0]  = %ld, naxes[1] = %ld\n", naxes[0], naxes[1]); */
    /*   exit(1); */
    /* } */
  }


  /* set up for multi-threads. This must be first FFTW call of any kind */
  /*  result = fftw_init_threads();
  if (result==0) {
    fprintf(stderr,"FATAL: fftw thread init failed.\n");
    exit(1);
  }
  fftw_plan_with_nthreads(n_threads);*/
  result=0;

  /* create output file now before doing any work in case it fails */
  if ((outfile = fopen(outfilename, "w")) == NULL) {
    fprintf(stderr,"FATAL: failed to open outfile file name <%s>.\n",
	    outfilename);
    exit(1);
  }

  /* check source file now before doing any work in case it doesn't exist */
  if ( sourcefilename != NULL ) {

    if ( (sourcefile=fopen(sourcefilename,"r")) == NULL ) {
      fprintf(stderr,"FATAL: Failed to open source file <%s>.\n",
	      sourcefilename);
      exit(1);
    }
    if ( freq == 0.0 ) {
      fprintf(stderr,"FATAL: The observing frequency (in MHz) is " \
	      "required when adding point sources.\n");
      exit(1);
    }
    if ( lst == 0.0 ) {
      fprintf(stderr,"FATAL: The observing lst (in rad) is " \
	      "required when adding point sources.\n");
      exit(1);
    }
    if ( lat == 0.0 ) {
      fprintf(stderr,"FATAL: The observing latitude (in rad) is " \
	      "required when adding point sources.\n");
      exit(1);
    }

    while ( fgets (line, 255, sourcefile) != NULL ) {
      if (line[0] == '#') continue;
      if (strlen (line) == 0) continue;
      status = sscanf( line, "%f %f %s %f %f %f", 
		       &param1, &param2, str, &fd, &ref_freq, &spec_index );
      if (status != 6) {
        fprintf(stderr,"FATAL: Invalid line in %s: '%s'", 
		sourcefilename, line );
        exit(1);
      }

      az      = realloc( az,      (src+1)*sizeof(float) );
      el      = realloc( el,      (src+1)*sizeof(float) );
      src_amp = realloc( src_amp, (src+1)*sizeof(float) );

      az[src]      = param1;
      el[src]      = param2;
      src_amp[src] = fd * pow( freq/ref_freq, spec_index );

      src++;

    }
    num_src = src;

    printf( "adding %d point sources to the visibilities, " \
	    "which will take around %.0f minutes with -O4\n",
            num_src, (2.9*num_src)/60.0 );

  }

  /* 
   * Read image file 
   */
  /* open it first */
  fits_open_image(&infile,infilename,READONLY, &result);
  if (result !=0) {
    fprintf(stderr, "FATAL: Failed to open file <%s> in FITSio. " \
	    "error code: %d\n", infilename, result);
    exit(1);
  }
  /* get image dimensions */
  fits_get_img_dim(infile, &naxis, &result);
  if (result != 0) {
    fprintf(stderr,"FATAL: Failed to get image naxis with error code %d\n",
	    result);
    exit(1);
  }
  if (debug) printf("Dimensionality of datacube in fitsfile: naxis = %d\n", 
		    naxis);
  /*
   * Check if number of extra dimension subscripts in cmdline
   * --axes option plus two image dimensions do not exceed the 
   * dimensionality (naxis) of the data hipercube in fits file 
   */
  // added by benkev, 2010-Jun-01

  nextradim = naxis - 2;  // Datacube in fitsfile dimensions above 2 of image
  if (nxdimfromcmd > nextradim) {
    fprintf(stderr, "FATAL: Number of subscripts in --axes option, %d, " \
	    "plus two image dimensions\nexceed actual data dimensionality, " \
	    "%d, in input fits file\n",
	    nxdimfromcmd, naxis);
    exit(1);    
  }
  
  //if (naxis!=2) {
  //  fprintf(stderr,"input image is NOT 2 dimensional. wrong wrong wrong.\n");
  //  exit(1);
  //}
  fits_get_img_size(infile, 16, naxes, &result);
  if (result != 0) {
    fprintf(stderr,"FATAL: Failed to get axis sizes with error code %d\n",
	    result);
    exit(1);
  }

  if (naxes[0] != naxes[1]) {
    fprintf(stderr,"FATAL: The input image MUST be square. " \
	    "Input was: %ld x %ld. exiting\n",naxes[0],naxes[1]);
    exit(1);
  }

  /* 
   * The image dimensions are (ncolumns x nrows) = (naxes[0] x naxes[1]).
   * In view of possible zero-padding of the brightness image,
   * two other variables for image dimensions are introduced,
   * n_cols and n_rows.
   * Added by L. Benkevitch, 22-Dec-2010
   */
  n_cols = naxes[0];
  n_rows = naxes[1];

  /* Currently, only even-sized crops are allowed */
  if (cropstr != NULL)  { /* No cropping requested */
    long nc, nr;
    if (cropx%2 != 0 || cropy%2 != 0) {
      printf("FATAL: wrong --crop dimensions; must be even-sized\n");
      exit(1);
    }
    /* Check crop sizes */
    nc = n_cols + 2*padzpix; /* padded with zeroz */
    nr = n_rows + 2*padzpix; /* padded with zeroz */
    if (cropx > nc || cropy > nc) {
      printf("FATAL: wrong --crop dimensions; cannot exceed " \
	     "(zero-padded) original\n");
      exit(1);
    }
  } /* if (cropstr != NULL) */
  else  { /* No cropping requested */
    cropx = n_cols + 2*padzpix; /* No cropping */ 
    cropy = n_rows + 2*padzpix; /* No cropping */ 
  } /* if (cropstr != NULL) {} else */


  // added by benkev:
  fits_read_key(infile, TSTRING, "BUNIT", units, NULL, &status);
  if (status!=0) 
    { fprintf(stderr,"FATAL: Error reading BUNIT from header.\n"); exit(1); }

  //if (debug) printf("\nnum_src = %d\n", num_src);

  //  if ( num_src > 0 ) {  /* commented out by benkev
  /* Now the parameters ra0, dec0, dl, dm are read from the fits 
   * file header unconditionally. If dl > 0, the data array read 
   * from the fits file will be right-left flipped to make RA grow 
   * from right to left according to the astronomical standards. 
   * Written in 2010 by L. Benkevitch */

  // DAM: read in data needed to generate the point sources.
  nkeys = status = 0;
  /* copy all the user keywords (not the structural keywords) */
  fits_get_hdrspace(infile, &nkeys, NULL, &status);
  if (status!=0) 
    { fprintf(stderr,"FATAL: Error reading header.\n"); exit(1); }
  fits_read_key(infile, TFLOAT, "CRVAL1", &ra0,  NULL, &status);
  if (status!=0) 
    { fprintf(stderr,"FATAL: Error reading CRVAL1 from header.\n"); exit(1); }
  fits_read_key(infile, TFLOAT, "CRVAL2", &dec0, NULL, &status);
  if (status!=0) 
    { fprintf(stderr,"FATAL: Error reading CRVAL2 from header.\n"); exit(1); }
  fits_read_key(infile, TFLOAT, "CDELT1", &dl,   NULL, &status);
  if (status!=0) 
    { fprintf(stderr,"FATAL: Error reading CDELT1 from header.\n"); exit(1); }
  fits_read_key(infile, TFLOAT, "CDELT2", &dm,   NULL, &status);
  if (status!=0) 
    { fprintf(stderr,"FATAL: Error reading CDELT2 from header.\n"); exit(1); }
  
  // added by benkev
  fits_read_key(infile, TSTRING, "BUNIT", units, NULL, &status);
  if (status!=0) 
    { fprintf(stderr,"FATAL: Error reading BUNIT from header.\n"); exit(1); }
  if (debug) printf("BUNIT (Flux units): %s\n", units);
  if (debug) printf("CRVAL1 (ra0, degrees): %f\n", ra0);
  if (debug) printf("CRVAL2 (dec0, degrees): %f\n", dec0);
  if (debug) printf("CDELT1 (dl, degrees): %f\n", dl);
  if (debug) printf("CDELT1 (dm, degrees): %f\n", dm);
  
  
  ra0  *= M_PI/180.0;
  dec0 *= M_PI/180.0;
  dl   *= M_PI/180.0;
  dm   *= M_PI/180.0;
  
  // added by benkev
  if (normalizer_specified) {
    if (debug) {
      printf("Normalizer is specified in command line: -n %f\n", 
	     normalizer); 
    }
  }
  else {
    short units_Jy_rad2 = 0;
    for (i = 0; i < 8; i++) {
      if (strcmp(units, strJy_rad2[i]) == 0) {
	normalizer = norm_Jy_rad2;
	units_Jy_rad2 = 1;
	break;
      }
    }
    if (debug) {
      if (units_Jy_rad2) 
	printf("Flux units are Jy/steradian; "	\
		"normalization constant is %f\n", normalizer); 
      else 
	printf("No flux units are specified in fits header; " \
		"normalization constant is assumed %f\n", normalizer); 
    }
  }

  ptsrc_normalizer = (n_cols*n_rows) / normalizer;

  if ( num_src > 0 ) {  // added by benkev
    
    // DAM: (n_rows-1)*dl = (n_cols-1)*dm = pi for the full sky. 
    du    = 1.0 / ( (n_rows-1)*dl );
    dv    = 1.0 / ( (n_cols-1)*dm );
    
    l = calloc( num_src, sizeof(float) );
    m = calloc( num_src, sizeof(float) );

    for( src=0; src<num_src; src++ ) {
      ra = lst - atan2(-sin(az[src])*cos(el[src]), cos(lat)*sin(el[src]) - 
		       sin(lat)*cos(el[src])*cos(az[src]));
      dec = asin(sin(lat)*sin(el[src]) + cos(lat)*cos(el[src])*cos(az[src]));
      //l[src] = cos(dec)*sin(ra-ra0);
      //m[src] = sin(dec)*cos(dec0) - cos(dec)*sin(dec0)*cos(ra-ra0);

      za = M_PI/2.0 - el[src];

      l[src] = sin(za)*sin(az[src]);
      m[src] = sin(za)*cos(az[src]);

      // this is to speed things up later.
      src_amp[src] *= ptsrc_normalizer;

    }
    if (debug) printf("Image size is %ld*%f x %ld*%f = %f x %f.\n", n_cols, 
		      dm, n_rows, dl, (n_cols-1)*dm, n_rows*dl);
  }


  time(&thetime);
  if (debug) printf("Original image size is %ld x %ld. %s",n_cols, 
		    n_rows, ctime(&thetime));

  /* make space for the data. Must allocate enough here for an in-place FFT
   * so there is extra padding in the last dimension of data. Use fftw_malloc 
   * so that things are word-aligned the way it needs it */
  arr_size = sizeof(double)*naxes[1]*2*(naxes[0]/2 + 1);
  data = fftw_malloc(arr_size);
  if (data == NULL) { 
    fprintf(stderr,"FATAL: malloc failed for array size %ld bytes\n", 
   	    arr_size);
    exit(1);
    }
  /*
   * Commented out by benkev 22-Dec-2010 (transferred closer to 
   * visibility file write):
   * mirror = fftw_malloc(sizeof(double)*n_cols*2);
   * if (mirror == NULL) {
   *   fprintf(stderr,"FATAL: malloc failed for array size %ld bytes\n", 
   *	    arr_size);
   * exit(1);
   * }
  */

  /* 
   * Read in the actual data. ********************
  */
  /*
   * Calculate the image plane position, imgstart, in fits file
   */
  // added by benkev, 2010-Jun-01
  /*
   * First, check if the subscripts in --axes do not exceed the
   * corresponding dimensions of the data hipercube. Simultaneously,
   * bring them to C index style (starting from 0, not 1).
   */
  for (i = 2; i < naxis; i++) {
    if (iax[i] > naxes[i]) {
    fprintf(stderr,"FATAL: a subscript in --axes, %d, exceeds corresponding " \
	    "dimension, %ld, of the datacube in fits file\n", 
	    iax[i], naxes[i]);
    exit(1);
    }
    else {
      iax[i]--; // C index style: starting from 0, not 1.
      //printf("iax[%d] = %d\n", i, iax[i]);
    }
  }


  if (nextradim != 0) { // I.e. the datacube has more than 2 axes
    long istart1; 
    imgstart = iax[naxis-1];
    istart1 = imgstart;
    //printf("Before loop: imgstart = %ld\n", imgstart);
    for (i = naxis-2; i > 1; i--) {  // Down to i = 2
      imgstart = iax[i] + imgstart*naxes[i];
      //printf("In     loop: imgstart = %ld = %d + %ld*%d\n", 
      //       imgstart, iax[i], istart1, naxes[i]);
      istart1 = imgstart;
    }
    //imgstart++;  // In fits indices start from 1, not 0.
    imgstart = 1 + imgstart*n_cols*n_rows;
  }
  if (debug) printf("Image in fitsfile datacube will be read " \
		    "from position %ld\n", imgstart);

  /* Note that even though the FITS images may contain single-precision   */
  /* float pixel values or even integer ones, this routine is reading     */
  /* the values into a double-precision array. Cfitsio automatically      */
  /* performs the datatype conversion in cases like this. See             */
  /* https://setisvn.ssl.berkeley.edu/trac/browser/lib/cfitsio/cookbook.c */
  //fits_read_img(infile, TDOUBLE, 1, n_cols*n_rows, NULL, data, 
  fits_read_img(infile, TDOUBLE, imgstart, n_cols*n_rows, NULL, data, 
		&anynull, &result );
  if (result !=0) {
    fprintf(stderr,"FATAL: reading FITS data failed with error code %d\n", 
	    result);
    exit(1);
  }
 
  time(&thetime);
  if (debug) printf("Success reading the data. %s", ctime(&thetime));
  fits_close_file(infile,&result);
  if (result !=0) 
    fprintf(stderr,"WARNING: fits_close_file() failed with code %d. " \
	    "Do we care?\n",
	    result);
  result=0;

  //added by benkev
  /*
   * Left-right image transposition
   *
   * Some programs (like LOsim) create brightness fits files with positive
   * CDELT1 parameter in header, which causes the image have East on the 
   * right and West on the left, which is non-standard. 
   * In this case the data array must be right-left flipped *
   * Added by L. Benkevitch, 2010.  
   *
   * 22 Feb 2011, modified by L. Benkevitch
   *   Current visgen operation right-left transposes the images generated
   * by standard sky simulation packages and by brigen. This is because
   * LOsim produced images with a non-standad RA direction (increasing 
   * from left to right), and visgen internally fixes this. 
   * Thus, I made the LR-transpose a normal mode here. The "-t" command
   * line option will leave the data as it is, but the visgen result will 
   * be transposed.
   */
  /* if (ewtranspose_specified) { */
  if (!ewtranspose_specified) {
    /* if (debug) printf("Right-left flipping the data.\n\n"); */
    halfaxis0 = n_cols/2;
    pincr = halfaxis0;
    if (n_cols%2) pincr++;  // Allow for ODD naxes
    pbeg = data;              // pointer to beginning of data row
    pmid = data + pincr;      // pointer to middle of data row
    pend = data + n_cols;   // pointer to end of data row
    for (i = 0; i < n_cols; i++) {
      while (pbeg < pmid) {
	temp_r = *pbeg;
	*(pbeg++) = *(--pend);
	*pend = temp_r;
      }
      pbeg += halfaxis0;
      pmid = pbeg + pincr;
      pend = pbeg + n_cols;
    }
  }
  else { /* Leaving the data as it is means LR-transposition! */
    if (debug) printf("Right-left flipping the data.\n\n");
  }


  /*
   * If the brightness zero-padding is requested (-p or --padzeropixels
   * command line option), we calculate the new array size, save data 
   * pointer in datx, allocate at data the new arrea, and copy from datx
   * into the middle of the new data.
   * Added by L. Benkevitch, 22-Dec-2010.  
   *
   */
  if (padzpix) {
    long offset, arr_size_dbls;
    int i = 0, j = 0;
    n_cols = n_cols + 2*padzpix;
    n_rows = n_rows + 2*padzpix;
    printf("Zero-padded image dimensions: %ld x %ld\n", n_cols, n_rows);
    datx = data; /* Save pointer to original image in datx */
    arr_size_dbls = n_rows*2*(n_cols/2+1); /* number of doubles */
    arr_size = sizeof(double)*arr_size_dbls;
    data = fftw_malloc(arr_size);
    if (data == NULL) { 
      fprintf(stderr,"FATAL: malloc failed for array size %ld bytes\n", 
	      arr_size);
      exit(1);
    }
    //printf("n_cols = %ld, n_rows = %ld\n", n_cols, n_rows);
    //printf("naxes[0] = %ld, naxes[1] = %ld\n", naxes[0], naxes[1]);
    //printf("padzpix = %ld\n", padzpix);
    for (i = 0; i < arr_size_dbls; i++) data[i] = 0.0;  /* Zeroize array */
    offset = (2*padzpix + naxes[0] + 1)*padzpix; /* top zero area */
    for (i = 0; i < naxes[1]; i++) {
      /* printf("i = %d, j = %d, offset = %ld\n", i, j, offset); */ 
      for (j = 0; j < naxes[0]; j++) { 
	data[offset+j+i*naxes[0]] = datx[j+i*naxes[0]];
      }
      offset = offset + 2*padzpix;
    }
    /* Debugging image output. Visualized by script viewim.py */
    /* FILE *fh; */
    /* size_t res, nelem; */
    /* int imgdims[2]; */
    /* imgdims[0] = naxes[0]; imgdims[1] = naxes[1]; */
    /* nelem = imgdims[0]*imgdims[1]; */
    /* printf("Original image dimensions: %ld x %ld\n",  */
    /* 	   naxes[0], naxes[1]); */
    /* fh = fopen("original_image.dat", "w"); */
    /* res = fwrite(imgdims,sizeof(int), 2, fh); /\* Header with 2 dims *\/ */
    /* res = fwrite(datx, sizeof(double), nelem, fh); /\* Write array *\/ */
    /* fclose(fh); */

    /* imgdims[0] = n_cols; imgdims[1] = n_rows; */
    /* nelem = imgdims[0]*imgdims[1]; */
    /* printf("Padded image dimensions: %ld x %ld\n",  */
    /* 	   n_cols, n_rows); */
    /* fh = fopen("padded_image.dat", "w"); */
    /* res = fwrite(imgdims,sizeof(int), 2, fh); /\* Header with 2 dims *\/ */
    /* res = fwrite(data, sizeof(double), nelem, fh); /\* Write array *\/ */
    /* fclose(fh); */
  }

  /* For FFT, must transpose and shift the image in the
   * X and Y axes by half the image (this has not been 
   * checked for odd sized images.) */
  /* transpose */
  for (i=1; i<n_cols; i++) {
    for (j=0; j<i; j++) {
      temp_r = data[(i+j*n_cols)  ];
      data[(i+j*n_cols)  ] = data[(j+i*n_cols)  ];
      data[(j+i*n_cols)  ] = temp_r;
    }
  }

  /* shift X */
  for (j=0; j<n_rows; j++) {
    int halfaxis = n_cols/2;
    for (i=0; i<halfaxis; i++) {
      temp_r = data[(i+halfaxis+j*n_cols)  ];
      data[(i+halfaxis+j*n_cols)  ] = data[(i+j*n_cols)  ];
      data[(i+j*n_cols)  ] = temp_r;
    }
  }

  /* shift Y */
  for (i=0; i<n_cols; i++) {
    int halfaxis = n_cols/2;
    for (j=0; j<n_rows/2; j++) {
      temp_r = data[(i+(halfaxis+j)*n_cols)  ];
      data[(i+(halfaxis+j)*n_cols)  ] = data[(i+j*n_cols)  ];
      data[(i+j*n_cols)  ] = temp_r;
    }
  }

  /* now re-arrange the data, so that there is the padding 
   * required by FFTW3 after each row */
  for (i=n_rows-1; i>0; i--) {
    memmove(data+i*2*(n_cols/2+1), data+i*n_cols, 
	    n_cols*sizeof(double));
  }

  //tempfile = fopen("testoutput.dat","w");
  //for (i=0;i<n_rows;i++) 
  //  fwrite(data+i*2*(n_cols/2+1),sizeof(double),n_cols,tempfile);
  //fclose(tempfile);

  /* FFT */
  /* make the cunning plan. Do an in-place transform. 
   * Use FFTW_ESTIMATE so that data
   * is not touched in planning */
  theplan = fftw_plan_dft_r2c_2d(n_cols, n_rows, data, 
				 (fftw_complex *) data, FFTW_ESTIMATE);
  if (theplan == NULL) {
    fprintf(stderr,"FATAL: fftw failed to make a plan\n");
    exit(1);
  }
  fftw_execute(theplan);   /* FFT  FFT  FFT  FFT  FFT  FFT  FFT  FFT  FFT */

  if ( num_src > 0 ) {
    time(&thetime);
    if (debug) printf("Done with FFT. Adding point sources. %s", 
		      ctime(&thetime));
  } else {
    time(&thetime);
    if (debug) printf("Done with FFT. Writing the data. %s", 
		      ctime(&thetime));
  }

  if ( num_src > 0 ) {
    // This is different to what's used below.
    arr_size = n_cols+2;

    // This is very time consuming. About 2.9 seconds per source 
    // with a 6602x6602 image 
    for(jpix=0; jpix<n_rows; jpix++) {
      j = jpix;
      jind = arr_size*jpix;
      if (jpix>=n_rows/2) j = jpix - n_rows;
      u = -2*M_PI*du*j;
      for(ipix=0; ipix<=n_cols/2; ipix++){
        i = ipix;
        ind = jind+2*ipix;
        if (ipix>=n_cols/2) i = ipix - n_cols;
        v = -2*M_PI*dv*i;
        for( src=0; src<num_src; src++ ) {
          phase = u*l[src] + v*m[src];
          data[ind    ] += src_amp[src]*cos(phase);
          data[ind + 1] += src_amp[src]*sin(phase);
        }
      }
    }

    time(&thetime);
    if (debug) printf("Done with point sources. Writing the data. %s",
		      ctime(&thetime));
  }

  /* 
   * Write the UV data 
   */
  /*
   * Determine crop boundaries for the saved image.
   * 
   * Added by L Benkevitch on 25 August 2011
   */
  /* start_col = (n_cols - cropx)/2; /\* n_cols and cropx: same parity! *\/ */
  /* end_col = start_col + cropx; */
  /* start_row = (n_rows - cropy)/2; /\* n_rows and cropy: same parity! *\/ */
  /* end_row = start_row + cropy; */
  /* wrstart = 2*start_col; /\* From what location in mirror to write *\/ */

  icol1 = cropx/2;
  icol2 = n_cols - icol1;
  irow1 = cropy/2;
  irow2 = n_rows - irow1;

  printf("icol1 = %ld, icol2 = %ld, irow1 = %ld, irow2 = %ld\n", 
	 icol1, icol2, irow1, irow2);
  printf("cropx = %ld, cropy = %ld\n", cropx, cropy);

  /* write header consisting of two ints. 
   * The first is the size of the image (number of pixels)
   * used to create the data (don't know why this is needed). 
   * The second is the size of the
   * UV data (number of pixels). */
  /* imgsize[0] = n_cols;  */
  imgsize[0] = cropx; 
  /* a time may come when this needs to pad the image either in UV 
   * or image space, in which case these will not be the same any more */
  /* imgsize[1] = n_cols; */
  imgsize[1] = cropy;

  result = fwrite(imgsize, sizeof(int), 2, outfile);

  printf("Image dimensions to be written: %ld x %ld\n", 
	 cropx, cropy);
	 /* n_cols, n_rows); */
  /*
   * The mirror allocation has been transferred here because 
   * at this point the zero-padded dimensions, n_col and
   * n_rows, are determined.
   * Added by L. Benkevitch, 22-Dec-2010
   */
  /*
   * The 'mirror' is an array containing one line of the resulting image.
   * It has n_cols (Re,Im) pairs or the image, so its structure is
   *    Re Im Re Im Re Im Re Im Re Im Re Im Re Im Re Im ... Re Im Re Im.
   * In the loop, the lines of 2D 'data' array are copied into 'mirror'
   * and the 'mirror' is saved in the output *.dat file.
   */
  mirror = fftw_malloc(sizeof(double)*n_cols*2);
  if (mirror == NULL) {
    fprintf(stderr,"FATAL: malloc failed for array size %ld bytes\n", 
	    sizeof(double)*n_cols*2);
    exit(1);
  }

  for(j = 0; j < 2*n_cols; j++) mirror[j]=0;  /* Zeroize the mirror line */

  arr_size = (n_cols/2 + 1);
  normalizer /= (n_cols*n_rows);


  /* for(j = start_row; j < end_row; j++) { */
  /*   for(i = start_col; i < end_col; i++) { */
  for(j = 0; j < n_rows; j++) {
    double conjugator = 1.0;
    long k, l;
    if (j < irow1 || j >= irow2) {
      for(i = 0; i < n_cols; i++) {

	if (i < icol1 || i >= icol2) {
	  iind = i;
	  jind = j;
	  if (i >= n_cols/2) {
	    iind = n_cols - i;
	    if (j != 0) jind = n_rows - j;
	    conjugator = -1.0;
	  }
	  k = 2*i;
	  l = 2*(jind*arr_size + iind);
	  mirror[k  ] = data[l  ]*normalizer;             /* Re */
	  mirror[k+1] = data[l+1]*conjugator*normalizer;  /* Im */
	  //printf("(%ld,%ld)%10.2e ", k, l, mirror[k  ]);
	}
      }

      //for(i = 0; i < icol1; i++) printf("[%ld]%g ", 2*i, mirror[2*i]);
      //printf("\n");
      result = fwrite(mirror, 2*sizeof(double), icol1, outfile);
      if (result != icol1) {
	fprintf(stderr, "ERROR: wrote %d dcomplex, expected to write %ld\n",
		result, icol1);
      }
      //printf("\n");
      //for(i = 0; i < icol1; i++) printf("[%ld]%g ", 2*(icol2+i), 
      //					mirror[2*(icol2+i)]);
      //printf("\n");
      result = fwrite(mirror+2*icol2, 2*sizeof(double), icol1, outfile);
      if (result != icol1) {
	fprintf(stderr, "ERROR: wrote %d dcomplex, expected to write %ld\n",
		result, icol1);
      }
      // printf("j = %d: fwrite\n", j);
    }
  }

  /* clean up */
  time(&thetime);

  if (debug) printf("done writing the data. %s", ctime(&thetime));
  fclose(outfile);

  fftw_destroy_plan(theplan);
  fftw_free(data);
  fftw_free(datx); /* Done right after allocation new padded data */
  fftw_free(mirror);
  free(l); free(m); free(az); free(el); free(src_amp);

  return 0;

}

//static void usage(char progname[]) {
//  fprintf( stderr, "usage: %s -i infile -o outfile "			
//	   "[-s srcfile -f frequency (MHz)] -t number_of_threads "
//	   "-n multiplicative_constant -d (debugging)\n",
//           progname);
//  exit(0);
//}

static void im2uv_help() {
  printf("Usage: MAPS_im2uv -i infile -o outfile [OPTIONS]...\n" \
	 "\t Convert brightness image from FITS file into visibility image " \
	 "using FFT.\n\n" \
	 "Example: MAPS_im2uv  -i cygnusXX.fits  -o cygnusXX.dat\n\n" \
	 "Options and interpretations:\n\n" \
	 "-i, --inputfilename\t input FITS file name\n" \
	 "-o, --outputfilename\t output file name (usually .dat)\n" \
	 "-s, --sourcefilename\t point source file name\n" \
	 "-f, --frequencyofsource\t frequency of source in MHz\n"    \
	 "-h, --lstofsource\t local siderial time in radians\n" \
	 "-l, --latitude\t\t observing latitude of point sources in radians\n" \
	 "-n, --normalyzer\t multiplicative constant to \n" \
	 "-d, --debug\t print work information\n"	\
	 "-m, --help\t print this screen\n" \
	 "-q, --query\t print table of all axis names and dimensions\n" \
	 "\t\t from the header of input fits file (at option -i)\n" \
	 "-a, --axes\t provide subscripts for miltidimensional fitsfiles;\n" \
	 "\t\t say, if file holds images for 12 different frequencies\n" \
	 "\t\t and 4 Stokes parameters for each, then the option\n" \
	 "\t\t --axes [5,2] asks the program to read image for 5th\n" \
	 "\t\t frequency and 2nd Stokes parameter\n" \
	 "-t, --ewtranspose\t flip image east-west before processing\n" \
	 "-p, --padzeropixels\t pad image with given number of zero pixels\n" \
	 "\t\t on each side\n\n" \
	 "More examples:\n\n"						\
	 "$ MAPS_im2uv.exe -i test_map.fits -o test_map.dat\n"		\
	 "\t Convert the brightness image 'test_map.fits' to\n" \
	 "\t visibilities (i.e., FFT it) and save\n" \
	 "\t it in the binary file 'test_map.dat'\n"	\
	 "\n"						\
	 "$ MAPS_im2uv.exe -i test4d_map.fits -q\n"	\
	 "  CTYPE1: RA---SIN,       2048\n"		\
	 "  CTYPE2: DEC--SIN,       2048\n"		\
	 "  CTYPE3: FREQ-LSR,       4\n"		\
	 "  CTYPE4: STOKES,         1\n"				\
	 "\t Print a table of all the axis names and dimensions from\n" \
	 "\t the file 'test4d_map.fits'.\n"			\
	 "\n"					\
	
	 "$ MAPS_im2uv.exe -i test4d_map.fits -o test4d_map.dat -a[3,1]\n" \
	 "\t Convert the brightness image of the 3rd frequency channel\n" \
	 "\t and first Stokes parameter from the input file \n" \
	 "\t 'test4d_map.fits' to visibilities and save the result\n" \
	 "\t in 'test4d_map.dat'.\n"				      \
	 "\n" \
	 "$ MAPS_im2uv.exe -i test4d_map.fits -o test4d_map.dat -a[3]\n" \
	 "\t Same action as in the previous example, because all the\n"	\
	 "\t dimensions omitted in the -a option are assumed to be\n"	\
	 "\t unity by default.\n"					\
	 "\n"								\
	 "$ MAPS_im2uv.exe -i modelsky.fits -o modelsky.dat -n 3.59e5\n" \
	 "\t Adjust the units of the brightness image from 'modelsky.fits',\n" \
	 "\t converting them to Jy/rad^2 via the multiplicative\n" \
	 "\t constant 3.59e5, using the -n option; then perform\n" \
	 "\t the FFT and save the resulting visibility data\n" \
	 "\t in 'modelsky.dat'.\n" \
	 "\n" \
	 "$ MAPS_im2uv.exe -i modelsky.fits -o modelsky.dat -p 512\n" \
	 "\t Pad the image in 'modelsky.fits' with the zero margins\n" \
	 "\t 512 pixels wide on top, bottom, left, and right sides,\n" \
	 "\t then perform the FFT and save the resulting visibility data\n" \
	 "\t in 'modelsky.dat'.\n"
	 "\n"								\
	 "MAPS_im2uv.exe -i modelsky.fits -o modelsky.dat -c 800\n"
	 "\t Make FFT of the brightness image of original size, then save\n" \
	 "\t the 800x800 central crop of resulting visibility data\n" \
	 "\t in 'modelsky.dat'.\n" \
	 "\n"						\
	 "MAPS_im2uv.exe -i sky.fits -o sky.dat -p 512 -c[300,200]\n" \
	 "\t Pad the image in 'sky.fits' with the zero margins of\n" \
	 "\t 512 pixels wide on top, bottom, left, and right sides,\n"	\
	 "\t perform the FFT, and  save the 300x200 central crop\n" \
	 "\t of resulting visibility data in sky.dat.\n");  
  exit(0);
}

