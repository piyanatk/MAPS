// Size parameters for large file access

#define FILES_MAX 10 // number of large files that can be open at once

//Characteristics of each file
#define LARGEFILE_MAX_SEGMENTSIZE 2147483648LL //2^31 // 1073741824LL // 2^30     2147483648LL // 2^31 (ext2fs doesn't support 2^31!)
#define LARGEFILE_MAX_SEGMENTS 100
#define LARGEFILE_MAX_FILESIZE (LARGEFILE_MAX_SEGMENTSIZE*LARGEFILE_MAX_SEGMENTS)

