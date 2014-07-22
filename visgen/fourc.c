/* fourc - slightly modified version of fourn, which was originally
 * written by Norm Brenner and found its way into Numerical Recipes.
 *
 * Replaces data by its ndim dimensional discrete Fourier transform,
 * if isign =1, and by its inverse FFT if isign = -1. nn [0..ndim-1]
 * is a vector containing the size of the input data array in each
 * dimension, all of which must be powers of two. Data is a complex
 * array with length equal to the product of the dimensional lengths.
 *
 *
 * Modified to use the complex data type, to be zero relative,
 * and beautified the code.                      rjc 2002.7.23   */

#include <stdio.h>
#include <math.h>
#include "comp_func.h"

void fourc (complex data[],         /* input and output array */
            int nn[],               /* lengths in each dimension */
            int ndim,               /* number of dimensions */
            int isign)              /* +1 forward, -1 reverse transform */
    {
    int idim, i1, i2, i3, 
        i2rev, i3rev, 
        ip1, ip2, ip3, ifp1, ifp2,
        ibit, k1, k2, n, nprev, nrem, ntot;
    
    double theta,
           wtemp;

    complex temp,
            w,
            wp;

                            /* compute total number of data points as the
                             * product of the lengths in each dimension */
    for (ntot=1,idim=0; idim<ndim; idim++)
        ntot *=  nn[idim];

    nprev = 1;
                                    /* loop over the dimensions */
    for (idim=ndim-1; idim>=0; idim--)
        {
        n = nn[idim];
        nrem = ntot / (n * nprev);
        ip1 = nprev;
        ip2 = ip1 * n;
        ip3 = ip2 * nrem;
        i2rev = 0;
                                    /* do bit reversal */
        for (i2=0; i2<ip2; i2+=ip1)
            {
            if (i2 < i2rev) 
                {
                for (i1=i2; i1<i2+ip1; i1++)
                    {
                    for (i3=i1; i3<=ip3; i3+=ip2)
                        {
                        i3rev = i2rev + i3 - i2;
                        temp = data[i3]; 
                        data[i3] = data[i3rev];
                        data[i3rev] = temp;
                        }
                    }
                }
            ibit = ip2 / 2;
            while (ibit >=  ip1 && i2rev >= ibit)
                {
                i2rev -=  ibit;
                ibit /= 2;
                }
            i2rev +=  ibit;
            }
                                    /* Danielson-Lanczos section - does
                                     * actual calculations */
        ifp1 = ip1;
        while (ifp1 < ip2)
            {
            ifp2 = ifp1 * 2;
            
            theta = isign * 6.28318530717959 / (ifp2 / ip1);
            wtemp = sin (0.5 * theta);
#ifdef _COMPLEX_H
            wp = -2.0 * wtemp * wtemp + I*sin(theta);
#else
            wp.re  =  -2.0 * wtemp * wtemp;
            wp.im = sin (theta);
#endif
            
            w = c_make(1.0,0.0);
            for (i3=0; i3<ifp1; i3+=ip1)
                {
                for (i1=i3; i1<i3+ip1; i1++)
                    {
                    for (i2=i1; i2<ip3; i2+=ifp2)
                        {
                        k1 = i2;
                        k2 = k1 + ifp1;
#ifdef _COMPLEX_H
                        temp = c_mult (w, data[k2]);
                        data[k2] = c_sub (data[k1], temp);
                        data[k1] = c_add (data[k1], temp);
#else
                        /* explicitly write out complex arithmetic for speed's sake */
                        temp.re = w.re*data[k2].re - w.im*data[k2].im;
                        temp.im = w.im*data[k2].re + w.re*data[k2].im;
                        data[k2].re = data[k1].re - temp.re;
                        data[k2].im = data[k1].im - temp.im;
                        data[k1].re = data[k1].re + temp.re;
                        data[k1].im = data[k1].im + temp.im;
#endif
                        }
                    }
		/* trigonometric recursion */
                temp = w;
                w = c_mult (w, wp);
                w = c_add (w, temp);
                }
            ifp1 = ifp2;
            }
        nprev *=  n;
        }
    }
