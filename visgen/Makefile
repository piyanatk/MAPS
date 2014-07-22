# Makefile for visgen                   CJL Jan 31 2002
#

LIBS=$(SIM_LIB)/libutil.a
CFLAGS=-g -O -Wall -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE  \
	-I/usr/include/netcdf-3/ -L/usr/lib/netcdf-3/ -L/usr/lib64/netcdf-3/ \
	-I$(CFITSIO_INC) -L$(CFITSIO_LIB) -L$(SIM_LIB)


.c.o:
#	mpicc $(CFLAGS) -c -I$(SIM_INC) $<
	gcc $(CFLAGS) -c -I$(SIM_INC) $<

SOURCES=\
compute_beamconvl.c  integrate.c      read_array_spec.c  write_main_header.c\
ionosphere.c         write_stnfile.c  \
compute_uvw.c        oob.c         write_timeblock_header.c\
convolve.c           open_gridfile.c  read_rfimodel.c    read_obs_spec.c\
get_patch.c          parse_cmdline.c  sizes.c            compute_gha.c\
get_site.c           patch_size.c\
read_stn_layouts.c interp_vis.c\
fourc.c          makebeam.c         add_noise.c \
fill_beamgrid.c      file_offset.c    cspline.c \
pol_response.c   coord_trans.c check_baseline.c \
get_patch_dummy.c visgen.c get_array_border.c \
convex_hull_2d.c list_proc.c

OBJECTS=\
compute_beamconvl.o  integrate.o      read_array_spec.o  write_main_header.o\
ionosphere.o  write_stnfile.o\
compute_uvw.o        oob.o         write_timeblock_header.o\
convolve.o	         open_gridfile.o  read_rfimodel.o    read_obs_spec.o\
get_patch.o          parse_cmdline.o  sizes.o            compute_gha.o\
get_site.o           patch_size.o\
read_stn_layouts.o interp_vis.o\
fourc.o makebeam.o add_noise.o  \
fill_beamgrid.o      file_offset.o    cspline.o \
pol_response.o   coord_trans.o check_baseline.o \
get_patch_dummy.o visgen.o get_array_border.o \
convex_hull_2d.o list_proc.o


visgen: $(OBJECTS) $(LIBS)
	mpicc $(CFLAGS) $(OBJECTS) -o visgen $(LIBS) -lnetcdf -lnovas_c -lcfitsio -lm
#	cc -o visgen $(OBJECTS) $(LIBS) -lm

install: visgen
	chmod 0775 visgen
	chmod a+x  scripts/*.py
	cp scripts/*.py $(SIM_BIN)
	mv visgen $(SIM_BIN)

clean:
	rm -f *.o
	rm -rf $(SIM_BIN)/plot_hull.py


