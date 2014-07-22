#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include "uvfits.h"
#include <string.h>

/*
 * Reads a uvfits file into a C data structure.
 *
 * For testing purposes, it also writes that data to a file and reads it in again,
 * to make sure that the reading and writing process does not lose any information.
 */
int main(int argc, char *argv[]){
    int status=0;
    uvdata *data, *data2;
    
    if (argc<2) {
        fprintf(stderr,"Usage: %s filename\n",argv[0]);
        exit(0);
    }
    
    uvfitsSetDebugLevel(1);
    status = readUVFITS(argv[1],&data);
    if (status != 0) {
        fprintf(stderr,"readUVFITS returned status %d\n",status);
        exit(1);
    }
    
    
    printAntennaData(data,stdout);
    // TODO: print FQ data here also
    printUVData(data,stdout);
    
    
    
    // Print to a file so can diff it
    FILE *fp = fopen("final.uvfits", "w");
    printUVData(data, fp);
    
    // Write from one data structure, read into another
    writeUVFITS("temp.uvfits", data);
    readUVFITS("temp.uvfits", &data2);
    
    // Use comparison function
    if(!compareUVFITS(data, data2)) {
        printf("read/written data is the SAME\n");
    }
    
    freeUVFITSdata(data);
    freeUVFITSdata(data2);
    
    return 0;
}
