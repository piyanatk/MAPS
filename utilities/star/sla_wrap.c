/***************************************************************************
 *
 *   Copyright (C) 2005 by Willem van Straten
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

/*
  This file provides C wrappers to the Fortran SLALIB bindings.
*/

/* *********************************************************************** */

/* Define to a macro mangling the given C identifier (in lower and upper
   case), which must not contain underscores, for linking with Fortran. */

   /* NOTE: all fortran subroutine names are made lower case by GNU C. */

// #define F77_FUNC_(name,NAME) name ## _

/* As F77_FUNC, but for C identifiers containing underscores. */
#define F77_FUNC_(name,NAME) name ## __ 

#define SLA_eqeqx F77_FUNC_(sla_eqeqx,slaEQEQX)

double SLA_eqeqx (double *);

double slaEQEQX (double mjd) {

  return SLA_eqeqx(&mjd);
}

#define SLA_map F77_FUNC_(sla_map,slaMAP)

void SLA_map (double *,double *,double *,double *,
		double *, double *, double *, double *,
		double *, double *);

void slaMAP (double rm, double dm, double pr, double pd,
		double px, double rv, double eq, double mjd, 
		double *ra, double *da) {

  SLA_map(&rm,&dm,&pr,&pd,&px,&rv,&eq,&mjd,ra,da);

}
#define SLA_caldj F77_FUNC_(sla_caldj,slaCald)

void SLA_caldj (int *, int *, int *, double *, int *);

void slaCaldj (int year, int month, int day, double *mjd, int *status) {
  SLA_caldj(&year,&month,&day,mjd,status);
}

#define SLA_altaz F77_FUNC_(sla_altaz,slaAltaz)

void SLA_altaz (double*, double*, double*, double*, double*, double*,
		double*, double*, double*, double*, double*, double*);

void slaAltaz (double ha, double dec, double phi,
                double *az, double *azd, double *azdd,
                double *el, double *eld, double *eldd,
                double *pa, double *pad, double *padd)
{
  SLA_altaz (&ha, &dec, &phi,
	     az, azd, azdd,
	     el, eld, eldd,
	     pa, pad, padd);
}


/* *********************************************************************** */

#define SLA_dcs2c F77_FUNC_(sla_dcs2c,slaDcs2c)

void SLA_dcs2c (double* ra, double* dec, double* v);

void slaDcs2c (double a, double b, double v[3])
{
  SLA_dcs2c (&a, &b, v);
}

/* *********************************************************************** */

#define SLA_de2h F77_FUNC_(sla_de2h,slaDE2H)
// fortran fake function prototype
void sla_de2h (double *ha,double *dec, double *lat, double *az, double *el);
// C wrapper function
void slaDE2H (double ha,double dec, double lat, double *az, double *el)
{
  SLA_de2h (&ha,&dec,&lat,az,el);
}


/* *********************************************************************** */

#define SLA_dmxv F77_FUNC_(sla_dmxv,slaDmxv)

void SLA_dmxv (double* dm, double* va, double* vb);

void slaDmxv (double dm[3][3], double va[3], double vb[3])
{
  SLA_dmxv (&(dm[0][0]), va, vb);
}

/* *********************************************************************** */

#define SLA_dsep F77_FUNC_(sla_dsep,slaDsep)

double SLA_dsep (double *, double *, double *, double*);

double slaDsep (double a1, double b1, double a2, double b2)
{
  return SLA_dsep (&a1, &b1, &a2, &b2);
}

/* *********************************************************************** */

#define SLA_dtt F77_FUNC_(sla_dtt,slaDtt)

double SLA_dtt (double* mjd);

double slaDtt (double dju)
{
  return SLA_dtt (&dju);
}

/* *********************************************************************** */

#define SLA_dvdv F77_FUNC_(sla_dvdv,slaDvdv)

double SLA_dvdv (double* va, double* vb);

double slaDvdv (double va[3], double vb[3])
{
  return SLA_dvdv (va, vb);
}

/* *********************************************************************** */

#define SLA_e2h F77_FUNC_(sla_e2h,slaE2H)
// fortran fake function prototype
void sla_e2h (float *ha,float *dec, float *lat, float *az, float *el);
// C wrapper function
void slaE2H (float ha,float dec, float lat, float *az, float *el)
{
  SLA_e2h (&ha,&dec,&lat,az,el);
}

/* *********************************************************************** */

#define SLA_epj F77_FUNC_(sla_epj,slaEpj)

double SLA_epj (double* mjd);

double slaEpj (double date)
{
  return SLA_epj (&date);
}

/* *********************************************************************** */

#define SLA_eqgal F77_FUNC_(sla_eqgal,slaEqgal)

void SLA_eqgal (double *, double *, double *, double *);

void slaEqgal (double dr, double dd, double *dl, double *db)
{
  SLA_eqgal (&dr, &dd, dl, db);
}

/* *********************************************************************** */

#define SLA_evp F77_FUNC_(sla_evp,slaEvp)

void SLA_evp (double* tdb, double* ep,
              double* dvb, double* dpb,
              double* dvh, double* dph);

void slaEvp (double date, double deqx,
             double dvb[3], double dpb[3],
             double dvh[3], double dph[3])
{
  SLA_evp (&date, &deqx, dvb, dpb, dvh, dph);
}

/* *********************************************************************** */

#define SLA_galeq F77_FUNC_(sla_galeq,slaGaleq)

void SLA_galeq (double *, double *, double *, double *);

void slaGaleq (double dl, double db, double *dr, double *dd)
{
  SLA_galeq (&dl, &db, dr, dd);
}

/* *********************************************************************** */

#define SLA_geoc F77_FUNC_(sla_geoc,slaGeoc)

void SLA_geoc (double *, double *, double *, double *);

void slaGeoc (double dl, double db, double *dr, double *dd)
{
  SLA_geoc (&dl, &db, dr, dd);
}

/* *********************************************************************** */

#define SLA_gmst F77_FUNC_(sla_gmst,slaGmst)

double SLA_gmst (double* mjd);

double slaGmst (double ut1)
{
  return SLA_gmst (&ut1);
}

/* *********************************************************************** */

#define SLA_pa F77_FUNC_(sla_pa,slaPa)

double SLA_pa (double* HA, double* DEC, double* PHI);

double slaPa (double ha, double dec, double phi)
{
  return SLA_pa (&ha, &dec, &phi);
}

/* *********************************************************************** */

#define SLA_prec F77_FUNC_(sla_prec,slaPrec)

void SLA_prec (double* ep0, double* ep1, double* rmatp);

void slaPrec (double ep0, double ep1, double rmatp[3][3])
{
  SLA_prec (&ep0, &ep1, &(rmatp[0][0]));
}

