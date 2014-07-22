# Makefile for maps2uvfits. Created Jan 2007. RBW.
INCS=-I$(SIM_INC) -I$(CFITSIO_INC)
CFLAGS=-g -Wall $(INCS) -L$(CFITSIO_LIB) -D_FILE_OFFSET_BITS=64

maps2uvfits: maps2uvfits.c uvfits.c
	cc $(CFLAGS) -o maps2uvfits maps2uvfits.c uvfits.c -lcfitsio -lm

install: maps2uvfits
	mv maps2uvfits $(SIM_BIN)

clean:
	rm -f *.o test_readuvfits test_uvfits maps2uvfits

test_readuvfits: uvfits.c test_readuvfits.c
	cc $(CFLAGS) -o test_readuvfits test_readuvfits.c uvfits.c -lcfitsio -lm
#	mv test_readuvfits $(SIM_BIN)

test_uvfits: uvfits.c test_uvfits.c
	cc $(CFLAGS) -o test_uvfits test_uvfits.c uvfits.c -lcfitsio -lm

