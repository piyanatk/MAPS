#include <stdio.h>

void parseHProfile(double maxHeight, const char *filename, const char *newfilename)
{
  FILE *file, *newfile;
  double height, dgarbage[2];
  int density, garbage[11];
  char buffer[1024];

  if((file = fopen(filename, "r")) == NULL) {
    printf("Can't open %s\n", filename);
    exit(1);
  }
  newfile = fopen(newfilename, "w");
  fscanf(file, "%lf %d %lf %d %d %d %d %d %d %d %d %d %d %lf %d", &height, &density, &garbage[0], &dgarbage[0], &garbage[1], &garbage[2], &garbage[3], &garbage[4], 
         &garbage[5], &garbage[6], &garbage[7], &garbage[8], &garbage[9], &dgarbage[1], &garbage[10]);
  while(height < maxHeight) { 
    fprintf(newfile, "%d\n", density);
    fscanf(file, "%lf %d %lf %d %d %d %d %d %d %d %d %d %d %lf %d", &height, &density, &garbage[0], &dgarbage[0], &garbage[1], &garbage[2], &garbage[3], &garbage[4], 
         &garbage[5], &garbage[6], &garbage[7], &garbage[8], &garbage[9], &dgarbage[1], &garbage[10]);
  }
  fprintf(newfile, "%d\n", density);
  fclose(file);
  fclose(newfile);
}

int main()
{
  parseHProfile(1250, "h6.txt", "HeightProfile6.txt");
  return 0;
}


  
