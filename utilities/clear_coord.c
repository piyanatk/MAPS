#include "visibility.h"

void
clear_coord (/* crd) */
struct coord *crd)
    {
    crd->ra_hrs = -1;
    crd->ra_mins = 0;
    crd->ra_secs = 0.0;
    crd->dec_degs = -100;
    crd->dec_mins = 0;
    crd->dec_secs = 0.0;
    }
