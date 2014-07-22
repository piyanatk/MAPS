/************************************************************************/
/*                                                                      */
/* Effective area of a station.  The type of antenna, and the           */
/* effective area formula to use, are determined by the observing       */
/* frequency, and the antenna type field in the stsp struct             */
/*                                                                      */
/*      Inputs:     stsp        Station struct                          */
/*                  freq        Observing frequency of this visibility  */
/*                                                                      */
/*      Output:     aeff        Area in square meters                   */
/*                  return value  0=OK, else bad                        */
/*                                                                      */
/* Created 27 June 2002 by CJL                                          */
/*                                                                      */
/************************************************************************/
#include <stdio.h>
#include "array_specification.h"

int
compute_aeff (/* stsp, freq, aeff)*/
struct stationspec *stsp,
double freq,
double *aeff)
    {
    int i, f, below, above;
    double ll, lh, hl, hh, lm, hm, value, totvalue;
                                        // Examine all antennas contributing
                                        // for this observing frequency    
    for (i=0; i<stsp->nant; i++)
        {
                                        // Point to information on this antenna type
        antenna = stsp->layout->antenna + i;
        ainf = ant_types + antenna->type;
                                        // Skip if freq not in antenna range
        if ((freq < ainf->minfreq) || (freq > ainf->maxfreq)) continue;
                                        // Find entry in elevation lookup table
                                        // zeroth entry corresponds to elevation=0,
                                        // increments are 5 degrees, 19 entries
        elev = stsp->elevation * 57.29578;
        if (elev < 0.0 || elev > 90.0)
            {
            msg ("Elevation out of range in compute_aeff", 2);
            return (-1);
            }
        low_el = floor (elev / 5.0);
        if (low_el == 18) low_el--;     // Handle special case of 90 degrees
        high_el = low_el + 1;
                                        // Linear interpolator
        eldiff1 = elev - (low_el * 5.0);
        eldiff2 = (high_el * 5.0) - elev;
                                        // Find bracketing frequency entries
        if (ainf->nfreq <= 0)
            {
            msg ("Missing elevation tables for antenna type %d", 2, antenna->type);
            return (-1);
            }
        below = above = -1;
        for (f=0; f<ainf->nfreq; f++)
            {
            if (ainf->freq[f].freq < freq) below = f;
            if (ainf->freq[f].freq > freq)
                {
                above = f;
                break;
                }
            }
        if (below == -1) below = above;
        if (above == -1) above = below;
                                        // Linear interpolator
        if (above == below) fdiff1 = fdiff2 = 1.0;
        else
            {
            fdiff1 = freq - ainf->freq[below].freq;
            fdiff2 = ainf->freq[above].freq - freq;
            }
                                        // Compute value
        ll = ainf->freq[below].el[low_el];
        lh = ainf->freq[below].el[high_el];
        hl = ainf->freq[above].el[low_el];
        hh = ainf->freq[above].el[high_el];
                                        // Elevation interpolation
        lm = (ll * eldiff2 + lh * eldiff1) / (eldiff1 + eldiff2);
        hm = (hl * eldiff2 + hh * eldiff1) / (eldiff1 + eldiff2);
                                        // Frequency interpolation
        value = (lm * fdiff2 + hm * fdiff1) / (fdiff1 + fdiff2);
                                        // Add together all antenna values
        totvalue += value;
        nant_used++;
        }

    meanvalue = totvalue / (double)nant_used;
