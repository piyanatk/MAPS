/************************************************************************/
/*                                                                      */
/* Mechanism for modifying the behaviour of various programs in a       */
/* uniform way.  In each program, one uses an extern declaration to     */
/* access the various environment variables the user may have set.      */
/* This routine should be called early by any program whose behaviour   */
/* is under environment variable control.  Any variable which has       */
/* potential or actual application to more than one program should be   */
/* processed here.                                                      */
/*                                                                      */
/*      Inputs:                                                         */
/*                                                                      */
/*      Output:                                                         */
/*                                                                      */
/* Modelled after Mk4 routine of same name, 7 Jan 2002 by CJL           */
/*                                                                      */
/* Being modified for MPI compatability, 25 Feb 2003 by RB		*/
/* 									*/
/************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define TRUE 1
#define FALSE 0
#define MAX_LINE 256

                                        /* Declare these global, to be */
                                        /* initialized here, but used */
                                        /* throughout program */
                                        
/* FIXME: this entire design for the environment variables and input/output directories needs a
   major cleanup */

char control_dir[FILENAME_MAX];
char output_dir[FILENAME_MAX];
char skymodel_dir[FILENAME_MAX];
char array_dir[FILENAME_MAX];
char rfi_dir[FILENAME_MAX];
char ionosphere_dir[FILENAME_MAX];
char tmpdir[FILENAME_MAX];
char textdir[FILENAME_MAX];
char visdir[FILENAME_MAX];
char layoutdir[FILENAME_MAX];

void
environment()
    {                                   /* Default values */
    FILE *envfile=NULL;
    char *dummy=NULL,line[MAX_LINE],envname[MAX_LINE],simenv[MAX_LINE],*ptr=NULL;
    int nline,flen;


    if ((dummy = getenv ("CONTROLDIR")) != NULL) strcpy (control_dir, dummy);
    else strcpy (control_dir, ".");

    if ((dummy = getenv ("OUTPUTDIR")) != NULL) strcpy (output_dir, dummy);
    else strcpy (output_dir, ".");

    if ((dummy = getenv ("SKYMODELDIR")) != NULL) strcpy (skymodel_dir, dummy);
    else strcpy (skymodel_dir, ".");

    if ((dummy = getenv ("ARRAYDIR")) != NULL) strcpy (array_dir, dummy);
    else strcpy (array_dir, ".");

    if ((dummy = getenv ("RFIDIR")) != NULL) strcpy (rfi_dir, dummy);
    else strcpy (rfi_dir, ".");

    if ((dummy = getenv ("IONOSPHEREDIR")) != NULL) strcpy (ionosphere_dir, dummy);
    else strcpy (ionosphere_dir, ".");

    if ((dummy = getenv ("TMPDIR")) != NULL) strcpy (tmpdir, dummy);
    else strcpy (tmpdir, ".");

    if ((dummy = getenv ("TEXTDIR")) != NULL) strcpy (textdir, dummy);
    else strcpy (textdir, ".");

    if ((dummy = getenv ("VISDIR")) != NULL) strcpy (visdir, dummy);
    else strcpy (visdir, ".");

    if ((dummy = getenv ("LAYOUTDIR")) != NULL) strcpy (layoutdir, dummy);
    else strcpy (layoutdir, ".");

/* this is what you'll do if master is excluded from processing */
/*					RB, 25 Feb 2003		*/

  if ((envfile = fopen("sim_setup.env", "r")) == NULL) {
      printf("can't open the env file \n");
      fflush(stdout);
      return;
  }
	
  while (TRUE)
  { if ((fgets (line, MAX_LINE-1, envfile)) == NULL) break; 
  	  // printf("line %d is %s \n", nline, line); fflush(stdout); 
	  nline++; 

	  if ((ptr = strchr(line, '=')) == NULL)
	  { ptr = line + strlen (line);
	  *(ptr+1) = '\0';
	  }

	  flen = ptr - line;
	  strncpy(envname, line, flen);
	  envname[flen] = '\0';

	  strcpy(simenv, ptr + 1);
	  // *(ptr+1) = '\0';
      flen = strlen(simenv);
      simenv[flen-1] = '\0';

  	  // printf("env in line %d is %s %s \n", nline, envname, simenv); fflush(stdout); 
  	  // printf("slen %d flen %d names %s %s \n", slen, flen, envname, simenv); fflush(stdout); 

    if ( strcmp(envname, "VISDIR") == 0) strcpy (visdir, simenv);
    if ( strcmp(envname, "TEXTDIR") == 0) strcpy (textdir, simenv);
    if ( strcmp(envname, "LAYOUTDIR") == 0) strcpy (layoutdir, simenv);
    if ( strcmp(envname, "CONTROLDIR") == 0) strcpy (control_dir, simenv);
    if ( strcmp(envname, "OUTPUTDIR") == 0) strcpy (output_dir, simenv);
    if ( strcmp(envname, "SKYMODELDIR") == 0) strcpy (skymodel_dir, simenv);
    if ( strcmp(envname, "ARRAYDIR") == 0) strcpy (array_dir, simenv);
    if ( strcmp(envname, "RFIDIR") == 0) strcpy (rfi_dir, simenv);
    if ( strcmp(envname, "IONOSPHEREDIR") == 0) strcpy (ionosphere_dir, simenv);
    if ( strcmp(envname, "TMPDIR") == 0) strcpy (tmpdir, simenv);
  } // close while loop
  fclose(envfile);

  //printf("read %d lines from the env file \n", nline); fflush(stdout); 
  // printf("env VISDIR is %s %d \n", visdir, strlen(visdir)); fflush(stdout); 
  // printf("filename is %s/dummy.vis ", visdir); fflush(stdout); 

    return;
}
