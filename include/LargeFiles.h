// Function prototypes for large file access
// The function xL has the same prototype and semantics as the function x (x=fopen, fclose, etc)

#include <stdio.h>
#include <sys/stat.h>
#include "DataStructures.h"

typedef int FileHandle;
FileHandle fopenL(const char *filename, const char *mode);
int fcloseL(FileHandle);
int statL(const char *filenameL, struct statL *buf)
;
// The following functions use 'long long' to represent the position, assuming that type is 64 bits
int fseekL(FileHandle, long long offset, int wherefrom);
long long ftellL(FileHandle);

size_t freadL (void *ptr, size_t element_size, size_t count, FileHandle);
size_t fwriteL(void *ptr, size_t element_size, size_t count, FileHandle);
int feofL(FileHandle);
