#include "visibility.h"

void
init_vb (struct visblock *vb)
    {
    vb->intg_no = -1;
    vb->station1 = -1;
    vb->station2 = -1;
    vb->u = 0.0;
    vb->v = 0.0;
    vb->w = 0.0;
    vb->nfreq = 0;
    }
