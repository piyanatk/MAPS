#include <stdio.h>

int main()
{
  int i;
  double curED;
  FILE *file, *outfile;

  file = fopen("HeightProfile1.txt", "r");
  outfile = fopen("MATLAB_profile.txt", "w");

  for(i = 0; i < 70; i++) fscanf(file, "%lf", &curED);
  fprintf(outfile, "[");
  for(i = 0; i < 1000; i++) fprintf(outfile, "%lf, ", 100 + (double) i/2);
  fprintf(outfile, "];\n");
  fprintf(outfile, "[");
  for(i = 0; i < 1000; i++) {
    fscanf(file, "%lf", &curED);
    fprintf(outfile, "%lf, ", curED);
  }
  fprintf(outfile, "];\n");
  fclose(file);
  fclose(outfile);
  return 0;
}
