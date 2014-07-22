// DataStructures.h: Defines special data types (structs) used in LOsim

// Author:	Peter Sherwood	sherwood@computer.org	(617) 244-0836

// 8/22/01	add Stokes parameters and spectral density to Source

#ifndef dataStructuresHaveBeenDefined
#include <sys/stat.h>

struct Complex
{
	double re;
	double im;
};

// File status
struct statL
{
	struct stat s;
	long long st_sizeL;
};

// All sources are Gaussian and elliptical
struct Source
{
	double fluxDensity;	// total flux
	double Q,U,V;		// Stokes parameters
	double specIx;		// spectral index
	double xOffset;
	double yOffset;
	double majorAxis;
	double minorAxis;
	double positionAngle;
	char id[10+1];
};

// A station consists of 1-MAX_ANTENNAEPERSTATION or more dipole antennae
struct Station
{
	double x,y,z;				// Cartesian coords relative to center of earth, in m
	double longi,lat;			// -pi<longitude<=pi (+=E)  -pi/2<latitude<=pi/2 (+=N)
	int nAnt;					// # antennae for this station
	unsigned short iAnt;		// the index of the 1st of the nAnt antennae in the antenna array
	char id[20+1];
};

struct Antenna
{
	struct Complex wt;			// weight
	double x,y,z;				// Cartesian coords relative to center of earth, in m
	char id[10+1];
};

#define dataStructuresHaveBeenDefined 1
#endif
