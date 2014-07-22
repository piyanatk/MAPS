/** OBSOLETE. Do not use. ***/
// Functions for large file access
// The function xL has the same prototype and semantics as the function x (x=fopen, fclose, etc)

#include "LargeFileLimits.h"
#include "LargeFiles.h"
#include <string.h>
#include <errno.h>

struct SegmentEntry
{
//	long length;	 // length of this segment (usually 2GB)
	FILE *fileptr;
};

static int hasBeenInitialized=0;
static char userFilename[FILES_MAX][FILENAME_MAX]; // name of the (simulated) file known to user
static char userMode[FILES_MAX][10]; // mode specified by user on fopenL()
static long long filePosition[FILES_MAX],
  fileLength[FILES_MAX]; // -1 for fileLength means this file is not in use
static struct SegmentEntry segmentEntry[FILES_MAX][LARGEFILE_MAX_SEGMENTS];

// Local function prototypes
FILE *fopenSegment(FileHandle fh, int i);
int FileHandleCheck(FileHandle fh);

// Get information on a simulated file
int statL(const char *filenameL, struct statL *buf)
{
	int iSeg,er;
	long long fileLength=0;
	char filename[FILENAME_MAX];
	struct stat fileInfo;	
	
	for (iSeg=0; iSeg<LARGEFILE_MAX_SEGMENTS; iSeg++)
	{
		sprintf(filename,"%s_%02d",filenameL,iSeg); // name of this segment is <username>_nn
		er=0;
		if(stat(filename,&fileInfo)<0) // get info about a single segment
		{
			er=errno;
			if(er==ENOENT) er=(iSeg==0); // past end segment-not an error; else file does not exist-pass on the error
			break;
		}
		if(iSeg==0) buf->s=fileInfo; // all info except size comes from 1st segment
		fileLength+=fileInfo.st_size;
	}
	buf->st_sizeL=fileLength;
	return er;
}

// Open a simulated file
// If successful, return the "file handle"; if not, return -1
// LIMITATION: mode must not be append
FileHandle fopenL(const char *filename, const char *mode)
{
	int i,nr,iSeg; FileHandle fh; long long lo,hi,t; char dummy[1];
	if(!hasBeenInitialized) // make all the files not in use the first time through
	{
		for (fh=0; fh<FILES_MAX; fh++) fileLength[fh]=-1;
		hasBeenInitialized=1;
	}
	// Assign a file handle
	for (fh=0; fh<FILES_MAX; fh++) if (fileLength[fh]<0) break;
	if(fh>=FILES_MAX) return -1; // no more files can be opened
	strcpy(userFilename[fh],filename); strcpy(userMode[fh],mode); // save the file info for later
	filePosition[fh]=fileLength[fh]=0;
	for (i=0; i<LARGEFILE_MAX_SEGMENTS; i++) segmentEntry[fh][i].fileptr=NULL; // mark all segments not in use
	
	// For read and append, we need to get file length; do so by attempting to open all possible segments
	// For write, do nothing; segment 0 will be opened when data is written
	if(userMode[fh][0]=='r' || userMode[fh][0]=='a')
	{
		// Locate the highest file segment and get its exact length
		for (iSeg=LARGEFILE_MAX_SEGMENTS-1; iSeg>=0; --iSeg) if(fopenSegment(fh,iSeg)!=NULL) break;
		if(iSeg>=0)
		{
			// The length is determined by binary search
			// However, unlike the standard binary search, we have no way of knowing when we've
			// found the search object (i e, read the file's last byte), so we depend on
			// the fact that after each ineration, either lo increases, or hi decreases. This last
			// is true because if lo<hi, (lo+hi)/2<hi due to integer division truncating.
			// lo <= file length <=hi
			lo=0; hi=LARGEFILE_MAX_SEGMENTSIZE;
			while(lo<hi)
			{
				t=(lo+hi)/2;
				fseek(segmentEntry[fh][iSeg].fileptr,t,SEEK_SET);
				nr=fread(dummy,1,1,segmentEntry[fh][iSeg].fileptr); // attempt a single-byte read
				if(nr) lo=t+1; // the file is at least t+1 bytes long (position t<=>length t+1)
				else hi=t; // the file is no longer than 't' bytes
			}
			fileLength[fh]=hi+(iSeg*LARGEFILE_MAX_SEGMENTSIZE);
			// Open the remaining segments that exist
			for (; iSeg>=0; --iSeg) fopenSegment(fh,iSeg);
			if(userMode[fh][0]=='a') filePosition[fh]=fileLength[fh];
		}
	}
	return fh;
}
	
FILE *fopenSegment(FileHandle fh, int i)
{
	FILE *fp;
	char segmentFilename[FILENAME_MAX],tmp[4];
	sprintf(tmp,"_%02d",i); // the suffix for the segment file
	strcpy(segmentFilename,userFilename[fh]); strcat(segmentFilename,tmp);
	fp=fopen(segmentFilename,userMode[fh]);
	segmentEntry[fh][i].fileptr=fp;
	return fp;
}

extern int fcloseL(FileHandle fh)
{
	int i,error=0;
	if((error=FileHandleCheck(fh))) return error;
	for (i=0; i<LARGEFILE_MAX_SEGMENTS; i++)
	{
		if(segmentEntry[fh][i].fileptr!=NULL) error|=fclose(segmentEntry[fh][i].fileptr);
	}
	fileLength[fh]=-1; // file no longer in use
	if(error) return 0;
	else return EOF;
}

// The following functions use 'long long' to represent the position, assuming that type is 64 bits
int fseekL(FileHandle fh, long long offset, int wherefrom)
{
	long long newPosition; int error;
	if((error=FileHandleCheck(fh))) return error;
	switch(wherefrom)
	{
		case(SEEK_SET):	if(offset<0) return 3; // file positions are >=0
				newPosition=offset; break;
		case(SEEK_CUR): newPosition=filePosition[fh]+offset; break;
		case(SEEK_END): newPosition=fileLength[fh]; break;
		default:	return 1; // error
	}
	
	// Check for a valid file position (can be beoyond previous EOF)
	if((newPosition)>LARGEFILE_MAX_FILESIZE) return 2; // can't make files this long
	if((newPosition)<0) return 4; // past start of file
	
	filePosition[fh]=newPosition;
	return 0;
}

long long ftellL(FileHandle fh)
{
	if(FileHandleCheck(fh)) return -1LL;
	return filePosition[fh];
}

size_t freadL (void *ptr, size_t element_size, size_t count, FileHandle fh)
{
	long long fPos,remSeg; size_t nRead=0; int iSeg,dummy;
	long long readSize=element_size*count; // # bytes user wants read
	size_t partSize; // # bytes to be read in this file segment
	if(FileHandleCheck(fh)) return 0;
	
	while(readSize>0)
	{
		iSeg=filePosition[fh]/LARGEFILE_MAX_SEGMENTSIZE;
		fPos=filePosition[fh]%LARGEFILE_MAX_SEGMENTSIZE;
		if(iSeg>=LARGEFILE_MAX_SEGMENTS) break; // error-tried to read file > maximum size
		if(segmentEntry[fh][iSeg].fileptr==NULL) break; // this segment not yet created
		if(filePosition[fh]>=fileLength[fh]) break; // trying to read past EOF
		if(iSeg<(fileLength[fh]/LARGEFILE_MAX_SEGMENTSIZE)) remSeg=LARGEFILE_MAX_SEGMENTSIZE-fPos;
		else remSeg=fileLength[fh]-filePosition[fh]; // we're in last segment, which may be partly full
		partSize=readSize; if(partSize>remSeg) partSize=remSeg; // bytes to read from this segment
		
		// If we read from a segment that hasn't been written, we read 0s with no error
		// (this simulates a "sparse file")
		if(segmentEntry[fh][iSeg].fileptr==NULL)
		{
			nRead+=partSize;
			memset(ptr,0,partSize);
 		}
		else
		{
			char errMsg[200]; int err,eof;
			fseek(segmentEntry[fh][iSeg].fileptr,fPos,SEEK_SET);
			nRead+=fread(ptr,1,partSize,segmentEntry[fh][iSeg].fileptr);
			err=ferror(segmentEntry[fh][iSeg].fileptr); eof=feof(segmentEntry[fh][iSeg].fileptr);
			if(err) {err=errno; strerror_r(err,errMsg,200);} // for debugging only
			dummy=err+0; // a dummy line to be able to set a breakpoint
		}
		filePosition[fh]=fPos+nRead+(iSeg*LARGEFILE_MAX_SEGMENTSIZE);
		readSize-=partSize; ptr=(void *)(((char *)ptr)+partSize); // advance for next part
	}
	return nRead/element_size;
}

// Write to the current file position (which may be beyond the physical EOF)
size_t fwriteL(void *ptr, size_t element_size, size_t count, FileHandle fh)
{
	long long fPos,remSeg; size_t nWritten=0; int iSeg;
	long long writeSize=element_size*count; // # bytes user wants written
	long partSize; // # bytes to be written in this file segment
	if(FileHandleCheck(fh)) return 0;
	fPos=filePosition[fh];
	
	while(writeSize>0)
	{
		iSeg=filePosition[fh]/LARGEFILE_MAX_SEGMENTSIZE;
		fPos=filePosition[fh]%LARGEFILE_MAX_SEGMENTSIZE;
		if(iSeg>=LARGEFILE_MAX_SEGMENTS) break; // error-tried to write file > maximum size
		if(segmentEntry[fh][iSeg].fileptr==NULL) fopenSegment(fh,iSeg); // this segment not yet created
		// Writing into an existing file (possibly extending it)
		fseek(segmentEntry[fh][iSeg].fileptr,fPos,SEEK_SET);
		remSeg=LARGEFILE_MAX_SEGMENTSIZE-fPos;
		partSize=writeSize; if(partSize>remSeg) partSize=remSeg;
		nWritten+=fwrite(ptr,1,partSize,segmentEntry[fh][iSeg].fileptr);
		filePosition[fh]=fPos+nWritten+(iSeg*LARGEFILE_MAX_SEGMENTSIZE);
		if(filePosition[fh]>fileLength[fh]) fileLength[fh]=filePosition[fh]; // file is extended
		writeSize-=partSize; ptr=(void *)(((char *)ptr)+partSize); // advance for next part
	}
	return nWritten/element_size;
}

int feofL(FileHandle fh)
{
	if(FileHandleCheck(fh)) return 0;
	return filePosition[fh]==fileLength[fh];
}

// Local functions
// Validate a file handle
int FileHandleCheck(FileHandle fh)
{
	if(fh<0 || fh>FILES_MAX) return 1; // illegal file handle
	if(fileLength[fh]<0) return 2; // file was not open
	return 0;
}
