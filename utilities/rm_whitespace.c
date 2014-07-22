/************************************************************************/
/*                                                                      */
/* Simple routine to trim whitespace from control file lines            */
/*                                                                      */
/*      Inputs:         line        Input character string              */
/*                                                                      */
/*      Output:         line        Duly sanitized                      */
/*                                                                      */
/* Created 1 Feb 2002 by CJL                                            */
/*                                                                      */
/************************************************************************/
#include <stdio.h>

#define TRUE 1
#define FALSE 0

void
rm_whitespace (/* line) */
char *line)
    {
    int nc, last_white;
    char *ptr;

    ptr = line;
                                        /* Trim leading whitespace */
    last_white = TRUE;
    nc = 0;
    while (*ptr != '\0')
        {
        switch (*ptr)
            {
                                        /* Compress all whitespace blocks to */
                                        /* single space character */
            case ' ':
            case '\t':
            case '\n':
            case '\r':
                if (! last_white)
                    {
                    line[nc++] = ' ';
                    last_white = TRUE;
                    }
                break;
                                        /* Copy all others as is */
            default:
                line[nc++] = *ptr;
                last_white = FALSE;
                break;
            }
        ptr++;
        }
                                        /* Strip trailing blank and terminate */
    if (nc > 0)
        if (line[nc-1] == ' ') nc--;
    line[nc] = '\0';
    }
