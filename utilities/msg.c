/*************************************************************/
/*                                                           */
/* Trivial routine to send messages to stderr, under control */
/* of "importance level" argument.  This allows for user     */
/* selectable verbosity levels for programs.                 */
/*                                                           */
/* Needs msglev and progname to be defined elsewhere         */
/*                                                           */
/* Stolen from Mk4, 7 Jan 2002 by CJL                        */
/*                                                           */
/*************************************************************/

#include <stdio.h>
#include <stdarg.h>

FILE *msg_fp=NULL;

void msg (char *string, int level, ...) 
{
  extern int msglev;
  extern char progname[];
  va_list args;

  if (msg_fp==NULL) msg_fp=stdout;
 
  if (level >= msglev)
    {
      /* printf("progname = %s, ptr = %p\n", progname, progname); */
      va_start (args,level);
      fprintf (msg_fp, "%s: ", progname);
      vfprintf (msg_fp, string, args);
      putc ('\n', msg_fp);
      va_end (args);
      fflush (msg_fp);
    }
}
