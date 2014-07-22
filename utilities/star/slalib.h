/***************************************************************************
 *
 *   Copyright (C) 2005 by Willem van Straten
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/
#ifndef __SLA_DECL_H
#define __SLA_DECL_H

/*
  These C wrappers are required in order that the PSRCHIVE software may
  be linked against either the C or Fortran version of SLALIB.

  Original wrappers by W van Straten
  Extracted from the PSRCHIVE library and added to MAPS by S. Ord 
 */
#ifdef __cplusplus
extern "C" {
#endif

double slaEQEQX(double mjd);

void slaMAP (double RA_mean, double DEC_mean, double proper_motion_RA,
		double proper_motion_DEC, double parallax, double rv, 
		double epoch, double mjd, double *ra, double *dec);

void slaCaldj (int year, int month, int day, double *mjd, int *status);
void slaAltaz ( double ha, double dec, double phi,
                double *az, double *azd, double *azdd,
                double *el, double *eld, double *eldd,
                double *pa, double *pad, double *padd );

void slaCs2c ( float a, float b, float v[3] ); /* spherical to cartesian - single precision */

void slaDcs2c ( double a, double b, double v[3] ); /* spherical to cartesian - double precision  */

void slaDE2H ( double ha, double dec, double phi,double *az,double *el);

void slaDmxv ( double dm[3][3], double va[3], double vb[3] );

double slaDsep ( double a1, double b1, double a2, double b2 );

double slaDtt ( double dju );

double slaDvdv ( double va[3], double vb[3] );

void slaE2H ( float ha, float dec, float phi,float *az,float *el);

double slaEpj ( double date );

void slaEqgal ( double dr, double dd, double *dl, double *db );

void slaEvp ( double date, double deqx,
              double dvb[3], double dpb[3],
              double dvh[3], double dph[3] );

void slaGaleq ( double dl, double db, double *dr, double *dd );

void slaGeoc ( double lat_radian, double ht_meters, double *r_au, double *z_au );

double slaGmst ( double ut1 );

double slaPa ( double ha, double dec, double phi );

void slaPrec ( double ep0, double ep1, double rmatp[3][3] );

#ifdef __cplusplus
           }
#endif

#endif

