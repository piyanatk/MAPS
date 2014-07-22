/************************************************************************/
/*                                                                      */
/* This routine prompts the user for confirmation, typically before     */
/* performing an operation that has been determined to be questionable. */
/*                                                                      */
/*      Inputs:         string          Message string                  */
/*                                                                      */
/*      Output:         return value    TRUE or FALSE                   */
/*                                                                      */
/* Stolen from Mk4 software, unmodified, 7 Jan 2002 by CJL              */
/*                                                                      */
/************************************************************************/

#include <stdio.h>

#define FALSE 0
#define TRUE 1

int
confirm ( /* string) */
char *string)
    {
    char buf[100];
    extern char progname[];             /* *progname no good here */
    static int interactive = TRUE;

    if (strcmp (string, "OFF") == 0) 
        {
        interactive = FALSE;
        return (TRUE);
        }
    if (strcmp (string, "ON")  == 0) 
        {
        interactive = TRUE;
        return (TRUE);
        }
                                        /* Batch mode, always do it */
    if (! interactive) return (TRUE);

    printf("%s: %s (y/n) : ", progname, string);
    while(TRUE) 
        {
        fgets (buf, 99, stdin);
        if(buf[0] == 'Y' || buf[0] == 'y') return(TRUE);
        else if(buf[0] == 'N' || buf[0] == 'n') return(FALSE);
        else printf("Answer 'y' or 'n' : ");
        }
    }
