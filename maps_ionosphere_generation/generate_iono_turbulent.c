#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fourc.c"
#include "iono-gen.h"

/*
 * A series of functions that generate an ionosphere with
 * electron content corresponding to the Kolmogorov Spectra
 *
 * Sean Ting 08 July 05
 */

static const double maxHeight = 1250;
static const double minHeight = 65;
static const double eVeloc = .000073;
static const double spacing = .5;

/*
 * The envelope function is multiplied with the general function for
 * the Kolmogorov Spectra in k-space to reduce the effect of waves
 * that are both too small and too large where these sizes are specified
 * by the o_scale (large wave limit) and i_scale (small wave limit)
 */

double envelope_function(double k, double o_scale, double i_scale)
{
  if(k < o_scale) return exp(k-o_scale);
  if(k > i_scale) return 0; 
  return 1; 
} 

void read_profile(int *eProfile, const char *filename)
{
  FILE *file;
  int i;

  if((file = fopen(filename, "r")) == NULL) {
    printf("Can't open %s\n", filename);
    exit(1);
  }
  printf("read_profile: there are %d steps\n",(int)((maxHeight - minHeight)/spacing + 1));
  for(i = 0; i < (int) ((maxHeight - minHeight)/spacing + 1); i++) {
    fscanf(file, "%d", &eProfile[i]);
    //printf("Got density %d\n",eProfile[i]);
  }
  fclose(file);
}

/*
 * A function that takes the result of the complex fourier transform in
 * k-space and converts it to an ionosphere with real electron content values. 
 */

 void make_real(double* newArray, complex *array, int t_dim, int x_dim, int y_dim, int z_dim, const char *profile, double start_alt, 
                double delta_alt, double cellspertime)
{
  int t, i, j, k, index,h_index;
  int *eProfile;

  eProfile = malloc(sizeof(int) * ((int) (maxHeight - minHeight)/spacing + 1));
  read_profile(eProfile, profile); 
  for(t = 0; t < t_dim; t++)
    for(i = 0; i < x_dim/*/2*/; i++) 
      for(j = 0; j < y_dim; j++)
        for(k = 0; k < z_dim; k++) {
            index = t*x_dim*y_dim*z_dim/*/2*/ + i*y_dim*z_dim + j*z_dim + k;
            h_index = (((start_alt-minHeight)/spacing + 1) + k*delta_alt/spacing);
            if(t==0 && i==0 && j==0) printf("k: %d. h_index: %d, val: %d\n",k,h_index,eProfile[h_index]);
	        //newArray[index] = (array[(i)*y_dim*z_dim + j*z_dim + k].re+1)*1e6*eProfile[h_index];
	        /* crude adjustment for reducing the dim of z to few pixels */
	        newArray[index] = (c_real(array[(i)*y_dim*z_dim + j*z_dim + k])+1)*1e6*eProfile[h_index]*400.0/(maxHeight - minHeight);        
        }
}


void make_real2D(double* newArray, complex *array, int t_dim, int x_dim, int y_dim, int z_dim,
                double start_alt,double delta_alt, double cellspertime)
{
  int t, i, j, k, index;
  const float tec=1e15;

  for(t = 0; t < t_dim; t++)
    for(i = 0; i < x_dim/*/2*/; i++) 
      for(j = 0; j < y_dim; j++)
        for(k = 0; k < z_dim; k++) {
            index = t*x_dim*y_dim*z_dim/*/2*/ + i*y_dim*z_dim + j*z_dim + k;
	        newArray[index] = (c_real(array[(i)*y_dim*z_dim + j*z_dim + k])+1)*tec;        
        }
}  

int TranslateCoord(int coord, int dim_length)
{
  coord -= floor(dim_length/2);
  if(coord < 0) coord += dim_length;
  return coord;
}

int TranslateCoord2(int coord, int dim_length)
{
  coord = dim_length - coord;
  if(coord == dim_length) coord = 0;
  return coord;
}

/*
 * Translates the values in the passed complex array so that the origin, which is considered
 * to lie at the element with coordinates (x_dim/2, y_dim/2, z_dim/2) lies at the origin.
 * This is done so that the array can be passed to the fourc function.
 */

void translate_array(int x_dim, int y_dim, int z_dim, complex **temp_iono)
{
  complex *temp, *new_iono = malloc(sizeof(complex)*x_dim*y_dim*z_dim);
  int i, j, k;
  int x, y, z;
  int index, nIndex;

  for(i = 0; i < x_dim; i++) 
    for(j = 0; j < y_dim; j++) 
      for(k = 0; k < z_dim; k++) {
        index = i*y_dim*z_dim + j*z_dim + k;
        nIndex = index;
        x = floor(((double) nIndex)/(y_dim*z_dim));  nIndex -= x*y_dim*z_dim;
        y = nIndex/z_dim; nIndex -= y*z_dim;
        z = nIndex;
        x = TranslateCoord(x, x_dim); y = TranslateCoord(y, y_dim); z = TranslateCoord(z, z_dim);
        nIndex = x*y_dim*z_dim + y*z_dim + z;
        new_iono[nIndex] = (*temp_iono)[index];
      }
  temp = *temp_iono;
  *temp_iono = new_iono;
  free(temp);
}

void print_darray(int x_dim, int y_dim, int z_dim, double *array)
{
  int i, j, k, index;
  int numNeg = 0, numPos = 0;
  double min = 0, max = 0;
  double maxNeg = 0, minPos = 0, sumPos = 0, sumNeg = 0; 
  for(i = 0; i < x_dim; i++)
    for(j = 0; j < y_dim; j++) 
      for(k = 0; k < z_dim; k++) {
        index = i*y_dim*z_dim + j*z_dim + k;
        if(array[index] < 0) {
          numNeg++;
          sumNeg += array[index];
          if(array[index] < min) min = array[index];
          if(array[index] > maxNeg || maxNeg == 0) maxNeg = array[index];
	} else {
          numPos++;
          sumPos += array[index];
          if(array[index] > max) max = array[index];
          if(array[index] < minPos || minPos == 0) minPos = array[index];
	}
        /*if(k == 0) printf("%lf at %d %d %d\n", array[index], i, j, k);*/
      }
  printf("\n");
  printf("Min: %f Max: %f Sum: %f\n", min, max, sumNeg + sumPos);
  printf("avgNeg: %f avgPos: %f rNeg: %f rPos: %f\n", sumNeg/numNeg, sumPos/numPos, min/maxNeg, max/minPos);
}

void print_carray(int x_dim, int y_dim, int z_dim, complex *array)
{
 int i, j, k, index;
 for(i = 0; i < x_dim; i++) {
    for(j = 0; j < y_dim; j++) { 
      for(k = 0; k < z_dim; k++) {
        index = i*y_dim*z_dim + j*z_dim+k;
        /*printf("Density is %f + %fi at %d %d %d\n", array[index].re, array[index].im, i, j, k);*/
        printf("(%f + %fi) ", c_real(array[index]), c_imag(array[index])); 
      }
    }
  printf("\n");
 }
 printf("\n");
}

/*
 * Generates a turbulent ionosphere by first creating a complex array representing the desired k-space distribution
 * of the ionpshere.  Values in the array are chosen by fitting the amplitude to the Kolmogorov Spectra value (c*k^-11/3) 
 * where k is the distance from the origin (the center of the array) with the envelope function and taking a random phase.  
 * The array is further forced to be Hermetian so that the output after a Fourier transformation is real.  The imaginary values
 * are then dropped and the array is copied into the array representing the ionosphere.
 */


void generate_iono_turbulent(int t_dim, int x_dim, int y_dim, int z_dim, double c, double o_scale, double i_scale, const char* profile, double delta_lat, double delta_lon, double start_alt, double delta_alt, double delta_time, double* iono_grid)
{
  complex *temp_iono;
  int i, j, k, index, tIndex;
  int dims[3], idum[1];
  double dist;
  double amplitude, phase;
  double cellspertime = eVeloc*delta_time/delta_lon;

  temp_iono = calloc(x_dim*y_dim*z_dim, sizeof(complex)); 
  printf("Please enter in a negative integer to be used as a seed for the random number generator: ");
  scanf("%d", idum);
  srandom(idum[0]);

  dims[0] = x_dim; dims[1] = y_dim; dims[2] = z_dim;
  for(i = 0; i < x_dim; i++) 
    for(j = 0; j < y_dim; j++)
      for(k = 0; k < z_dim; k++) {
        index = i*y_dim*z_dim + j*z_dim + k;
        tIndex = TranslateCoord2(i, x_dim)*y_dim*z_dim + TranslateCoord2(j, y_dim)*z_dim + TranslateCoord2(k, z_dim);
        if(index >= x_dim*y_dim*z_dim/2 + y_dim*z_dim/2 + z_dim/2) break;
        dist = pow(pow(x_dim/2-i, 2) + delta_lat*x_dim/(delta_lon*y_dim)*pow(y_dim/2-j, 2) + delta_lat*x_dim/(delta_alt*z_dim/R_EARTH)*pow(z_dim/2-k, 2), .5);  
        amplitude = envelope_function(dist, o_scale, i_scale)*c*(pow(dist, -11/3));
        phase = 2.0*M_PI*((double)random()/(double)RAND_MAX); 
        temp_iono[index] = s_mult(c_exp(phase),amplitude);
        temp_iono[tIndex] = c_conj(temp_iono[index]);
      }
  translate_array(x_dim, y_dim, z_dim, &temp_iono);
  fourc(temp_iono,dims, 3, -1);
  translate_array(x_dim, y_dim, z_dim, &temp_iono);
  make_real(iono_grid, temp_iono, t_dim, x_dim, y_dim, z_dim, profile, start_alt, delta_alt, cellspertime);
  free(temp_iono);
}
  
  
void generate_iono_turbulent2D(int t_dim, int x_dim, int y_dim, int z_dim, double c, double o_scale, double i_scale, double delta_lat, double delta_lon, double start_alt, double delta_alt, double delta_time, double* iono_grid)
{
  complex *temp_iono;
  int i, j, k, index, tIndex;
  int dims[3], idum[1];
  double dist;
  double amplitude, phase;
  double cellspertime = eVeloc*delta_time/delta_lon;

  temp_iono = calloc(x_dim*y_dim*z_dim, sizeof(complex)); 
  printf("Please enter in a negative integer to be used as a seed for the random number generator: ");
  scanf("%d", idum);
  srandom(idum[0]);

  dims[0] = x_dim; dims[1] = y_dim; dims[2] = z_dim;
  for(i = 0; i < x_dim; i++) 
    for(j = 0; j < y_dim; j++)
      for(k = 0; k < z_dim; k++) {
        index = i*y_dim*z_dim + j*z_dim + k;
        tIndex = TranslateCoord2(i, x_dim)*y_dim*z_dim + TranslateCoord2(j, y_dim)*z_dim + TranslateCoord2(k, z_dim);
        if(index >= x_dim*y_dim*z_dim/2 + y_dim*z_dim/2 + z_dim/2) break;
        dist = pow(pow(x_dim/2-i, 2) + delta_lat*x_dim/(delta_lon*y_dim)*pow(y_dim/2-j, 2) + delta_lat*x_dim/(delta_alt*z_dim/R_EARTH)*pow(z_dim/2-k, 2), .5);  
        amplitude = envelope_function(dist, o_scale, i_scale)*c*(pow(dist, -11/3));
        phase = 2.0*M_PI*((double)random()/(double)RAND_MAX); 
        temp_iono[index] = s_mult(c_exp(phase),amplitude);
        temp_iono[tIndex] = c_conj(temp_iono[index]);
      }
  translate_array(x_dim, y_dim, z_dim, &temp_iono);
  fourc(temp_iono,dims, 3, -1);
  translate_array(x_dim, y_dim, z_dim, &temp_iono);
  make_real2D(iono_grid, temp_iono, t_dim, x_dim, y_dim, z_dim, start_alt, delta_alt, cellspertime);
  free(temp_iono);
}


