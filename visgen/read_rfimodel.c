/************************************************************************/
/*                                                                      */
/* This routine reads a RFI model specification from a file, and places */
/* the information into an rfi_model struct for later use               */
/*                                                                      */
/*      Inputs:     filename        Contains RFI information            */
/*                                                                      */
/*      Output:     rfimodel        Struct containing the info          */
/*                  return value    0=OK, else bad                      */
/*                                                                      */
/* Created <date> by <author>                                             */
/*                                                                      */
/************************************************************************/
#include "rfi_model.h"
#include "utility.h"

int
read_rfimodel (/* filename, rfimodel) */
char *filename,
struct rfi_model *rfimodel)
    {
    msg ("read_rfimodel stubbed", 2);
    return (0);
    }
