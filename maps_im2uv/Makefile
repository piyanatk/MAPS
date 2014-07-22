INCS = -I$(CFITSIO_INC)
LIBS = -L$(CFITSIO_LIB)
CFLAGS= -Wall -g -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE

maps_im2uv: maps_im2uv.c
#	cc $(CFLAGS) -o maps_im2uv maps_im2uv.c $(INCS) $(LIBS) -lfftw3_threads -lfftw3 -lcfitsio -lm -lpthread 
	cc $(CFLAGS) -o maps_im2uv maps_im2uv.c $(INCS) $(LIBS) -lfftw3 -lcfitsio -lm

install: maps_im2uv
	mv maps_im2uv $(SIM_BIN)

clean:
	rm -f *.o maps_im2uv
