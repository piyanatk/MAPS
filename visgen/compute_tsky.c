/************************************************************************/
/*                                                                      */
/* Figures out the sky temperature contibution to system temperature,   */
/* based on the frequency, the observing parameters, and a map of the   */
/* galactic background radiation.                                       */
/* Currently a semi-stubbed version which assumes a uniform galactic    */
/* background                                                           */
/*                                                                      */
/*      Inputs:     freq            Frequency in MHz                    */
/*                  obsparam        Currently unused                    */
/*                                                                      */
/*      Output:     return value    Sky temperature (Kelvins)           */
/*                                                                      */
/* Created 9 September 2002 by CJL                                      */
/*                                                                      */
/************************************************************************/
#include <stdio.h>
#include <math.h>
#include "observing_param.h"

double
compute_tsky (// freq, obsparam)
double freq,
struct observing_param *obsparam)
    {
    double freq_ratio, tsky;
                                        // Assume 100K at 300 MHz, and simply
                                        // apply known galactic spectral index
    freq_ratio = freq / 300.0;
    tsky = 100.0 * pow (freq_ratio, -2.55);

    return (tsky);
    }
