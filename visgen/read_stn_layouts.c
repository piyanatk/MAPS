/************************************************************************/
/*                                                                      */
/* This routine reads and stores station layout descriptions            */
/* specified in text files.  The resulting structure array elements     */
/* are then pointed to in the array specification structure, allowing   */
/* convenient and relatively flexible specification of array and        */
/* station design.                                                      */
/*                                                                      */
/*      Inputs:     layoutdir       Directory containing files to read  */
/*                                                                      */
/*      Output:     layouts         Structure array filled in           */
/*                  return value    Number of layouts, <=0 if error     */
/*                                                                      */
/* Created 6 Feb 2002 by CJL                                            */
/*                                                                      */
/************************************************************************/
#include <string.h>
#include <stdio.h>
#include <sys/types.h>
#include <dirent.h>
#include <stdlib.h>
#include "array_specification.h"

#define FALSE 0
#define TRUE 1

#define LAYOUT_SUFFIX ".layout"
#define LINE_SIZE 256

/* private function prototypes */
int GetAntennaParams(struct antenna *ant, char *params);

struct station_layout *
get_stn_layout( char *layoutdir, char *layoutname) {
    char fullname[FILENAME_MAX];
    char line[LINE_SIZE];
    int n,  nant, ant_type,last_num_receptors=0;
    float north, east, height, gain,phase=0;
    struct station_layout *layout=NULL;
    struct dirent *ds=NULL;
    DIR *dp=NULL;
    FILE *fp=NULL;

    /* Open directory */
    if ((dp = opendir (layoutdir)) == NULL) {
	    msg ("Could not open station layout directory '%s'", 2, layoutdir);
	    return NULL;
    }
    
    /* Examine all files in directory */
    while ((ds = readdir (dp)) != NULL) {
        
        /* Ignore self and parent directories */
        if ((strcmp (ds->d_name, ".") == 0) || 
	    (strcmp (ds->d_name, "..") == 0))  continue;

        /* Look for names of form xyz.layout */
        if( strlen(ds->d_name) != strlen(LAYOUT_SUFFIX)+strlen(layoutname) ||
            strcmp(ds->d_name+strlen(ds->d_name)-strlen(LAYOUT_SUFFIX),
		   LAYOUT_SUFFIX)!=0 ||
            strncmp(layoutname,ds->d_name,strlen(layoutname)) != 0) {
            //msg("get_stn_layout: skipping file %s",-1,ds->d_name);
            continue;
        }

         /* Got the layout we were looking for - read and parse it */
        sprintf (fullname, "%s/%s", layoutdir, ds->d_name);
        if ((fp = fopen (fullname, "r")) == NULL){
	        msg ("Could not open file '%s', ignoring", 2, ds->d_name);
	        continue;
	    }
	    msg("Reading station layout file %s", 0, ds->d_name);

	    /* make space for the layout. no arbitrary maximum size */
	    layout = calloc(1,sizeof(struct station_layout));
	    if (layout==NULL) {
	        fprintf(stderr,"FATAL: alloc for new station layout failed\n");
	        exit(1);
	    }
	    layout->ants=NULL;
        strncpy(layout->layout_name,layoutname,MAX_LNAME_SIZE);

	    /* Loop through lines */
        nant = 0;
        last_num_receptors=0;
        while (fgets (line, LINE_SIZE, fp) != NULL) {
	  // Skip comment and blank lines
            if (line[0] == '*') continue;
	        if (line[0] == '#') continue;
            rm_whitespace (line);
            if (strlen (line) == 0) continue;
        
            /* Backward compatibility: station layout files 
	     * have an optional first line which re-states the
             * name of the layout. This is totally redundant, 
	     * so we just ignore this line if it exists. */
            if (strncmp (line, "NAME", 4) == 0) {
                continue;
            }
                             // Everything else should be 3 f.p. numbers
                             // and an antenna type ID number
	        north = east = height = phase = 0.0;
	        gain =1.0;
	        ant_type=ANT_ISOTROPIC;

	        /* NOTE: if this changes so that more or less 
		 * than 6 params are read, then you MUST
	         * also change the test below for GetAntennaParams 
		 * and the hard-coded number of
	         *  items to skip in GetAntennaParams */
            n = sscanf (line, "%f %f %f %d %f %f", &north, &east, &height, 
			&ant_type, &gain, &phase);
            if (n < 3) {
	      msg ("Invalid antenna line '%s' in file '%s'", 
		   2, line, fullname);
	      return NULL;
	    }
	    msg("antenna %d. type %d, (n,e,h): (%g,%g,%g), " \
		"gain: %g, phase: %g", -1, nant, ant_type, north, east, 
		height,gain,phase);
	    
	    /* make space for an antenna struct. no arbitrary maximum */
	    layout->ants = realloc(layout->ants,
				   sizeof(struct antenna)*(nant+1));
	    if (layout->ants ==NULL) {
	      fprintf(stderr,
		      "FATAL: failed to realloc space for anntena struct\n");
	      exit(1);
	    }
	    
        layout->ants[nant].north  = north;
        layout->ants[nant].east   = east;
        layout->ants[nant].height = height;
        layout->ants[nant].type   = ant_type;
        layout->ants[nant].gain[0]   = gain;
        layout->ants[nant].phase[0]  = phase;
        layout->ants[nant].params = NULL;
        layout->ants[nant].gain[1]   = 0.0;
        layout->ants[nant].phase[1]  = 0.0;
        layout->ants[nant].pol_type[0] = POL_I;
        layout->ants[nant].pol_type[1] = POL_I;
        layout->ants[nant].num_receptors = 1;
	    
	    /* special case for dual-receptor antennas */
	    if (ant_type == ANT_CROSSED_DIPOLES_ON_GROUNDPLANE ||
		    ant_type == ANT_CROSSED_DIPOLES_HORIZONTAL  ||
		    ant_type == ANT_IDEAL_PARABOLOID_DUAL_LINEAR ||
		    ant_type == ANT_IDEAL_PARABOLOID_DUAL_CIRCULAR) {
	      layout->ants[nant].num_receptors = 2;
	    }
	    
	    /* sanity check: all or none should be dual-pol receptors */
	    if (last_num_receptors != 0 && 
		last_num_receptors!=layout->ants[nant].num_receptors) {
	      msg("ERROR: mixed receptor types for layout name %s. exiting.",
		  2, fullname);
	      return NULL;
	    }
	    
	    last_num_receptors = layout->ants[nant].num_receptors;
	    
	    /* read the antenna-specific parameters if there 
	     * are more than 3 values in the line*/
	    if (n ==6) {
	      if (GetAntennaParams(layout->ants+nant,line) != 0) {
            msg("ERROR: failed to get antenna params. exiting.",2);
            return NULL;
	      }
	    }
	    
            nant++;
	}
	/* close the open file */
	fclose(fp);
        layout->nant = nant;
        
    }
    /* close the open dir */
    closedir(dp);
    return layout;
}


/********************
 Extract out the remainder of the antenna-specific parameters
 for this antenna type. All of the rest params are in the
 string params.
********************/
int GetAntennaParams(struct antenna *ant, char *param_str) {
  int result =0,n_extracted=0,n_params=0;

  switch (ant->type) {
  case ANT_SHORT_DIPOLE_ON_GROUNDPLANE:
    n_params = 2; /* PA, height above groundplane */
    if ((ant->params = calloc(n_params,sizeof(float))) ==NULL) return -1;
    n_extracted = sscanf(param_str,"%*s %*s %*s %*s %*s %*s %f %f",
			 ant->params,ant->params+1);
    msg("antenna extra params: %f %f",-1,ant->params[0],ant->params[1]);
    break;
  case ANT_GAUSSIAN:
  case ANT_IDEAL_PARABOLOID:
    n_params = 1;/* diameter (meters) */
    if ((ant->params = calloc(n_params,sizeof(float))) ==NULL) return -1;
    n_extracted = sscanf(param_str,"%*s %*s %*s %*s %*s %*s %f",ant->params);
    msg("antenna extra params: %f",-1,ant->params[0]);
    break;
  case ANT_CROSSED_DIPOLES_ON_GROUNDPLANE:
    /* special case: there are dual receptors on this antenna */
    n_params = 6; /* PA, height, gain2, phase2, PA2, height2 */
    if ((ant->params = calloc(n_params,sizeof(float))) ==NULL) return -1;
    n_extracted = sscanf(param_str,
			 "%*s %*s %*s %*s %*s %*s %f %f %f %f %f %f",
			 ant->params, ant->params+1, ant->gain+1, ant->phase+1,
			 ant->params+2,ant->params+3);
    msg("antenna extra params: %f %f %f %f %f %f",-1,
	ant->params[0],ant->params[1],ant->gain[1], ant->phase[1],
	ant->params[2],ant->params[3]);
    ant->pol_type[0] = POL_X;
    ant->pol_type[1] = POL_Y;
    break;
  case ANT_CROSSED_DIPOLES_HORIZONTAL:
    /* special case: there are dual receptors on this antenna */
    n_params = 4; /* PA, gain2, phase2, PA2 */
    if ((ant->params = calloc(n_params,sizeof(float))) ==NULL) return -1;
    n_extracted = sscanf(param_str,"%*s %*s %*s %*s %*s %*s %f %f %f %f",
			 ant->params,ant->gain+1, ant->phase+1,ant->params+1);
    msg("antenna extra params: %f %f %f %f",-1,
	ant->params[0],ant->gain[1], ant->phase[1],ant->params[1]);
    ant->pol_type[0] = POL_X;
    ant->pol_type[1] = POL_Y;
    break;
  case ANT_IDEAL_PARABOLOID_DUAL_LINEAR:
    /* special case: there are dual receptors on this antenna */
    n_params = 6; /* PA1, diam1, gain2, phase2, PA2, diam2 */
    if ((ant->params = calloc(n_params,sizeof(float))) ==NULL) return -1;
    n_extracted = sscanf(param_str,
			 "%*s %*s %*s %*s %*s %*s %f %f %f %f %f %f",
			 ant->params, ant->params+1, ant->gain+1, 
			 ant->phase+1, ant->params+2, ant->params+3);
    msg("antenna extra params: %f %f %f %f %f %f",-1,
	ant->params[0], ant->params[1], ant->gain[1], 
	ant->phase[1], ant->params[2], ant->params[3]);
    ant->pol_type[0] = POL_X;
    ant->pol_type[1] = POL_Y;
    break;
  case ANT_IDEAL_PARABOLOID_DUAL_CIRCULAR:
    /* special case: there are dual receptors on this antenna */
    n_params = 6; /* PA1, diam1, gain2, phase2, PA2, diam2 */
    if ((ant->params = calloc(n_params,sizeof(float))) ==NULL) return -1;
    n_extracted = sscanf(param_str,
			 "%*s %*s %*s %*s %*s %*s %f %f %f %f %f %f",
			 ant->params,ant->params+1,ant->gain+1, ant->phase+1,
			 ant->params+2,ant->params+3);
    msg("antenna extra params: %f %f %f %f %f %f",-1,
	ant->params[0],ant->params[1],ant->gain[1], ant->phase[1],
	ant->params[2],ant->params[3]);
    ant->pol_type[0] = POL_R;
    ant->pol_type[1] = POL_L;
    break;
  case ANT_ISOTROPIC:
    break;
  default:
    msg("GetAntennaParams: unknown antenna type %d.",2,ant->type);
    result = -1;
  }

  /* check that we got what we expected */
  if (n_params != n_extracted) {
    msg("ERROR: GetAntennaParams: expected %d params. "\
	"Extracted %d. String was: <%s>",
	2,n_params,n_extracted,param_str);
    result = -1;
  }
  return result;
}
