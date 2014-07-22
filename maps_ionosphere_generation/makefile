CC = gcc
FLAGS = -Wall -g
NETCDF_INCLUDE = -I$(SIM)/include -I/usr/include/netcdf-3
NETCDF_LIB = -L/usr/lib/netcdf-3 -L$(SIM_LIB) -lnetcdf -lutil

OBJS = generate_iono_ramp.o generate_iono_parabolic.o iono-gen.h generate_iono_turbulent.o generate_sine.o

.c.o:
	$(CC) $(FLAGS) $(NETCDF_INCLUDE) -c $< 

generate_ionoscreen: generate_ionoscreen.c $(OBJS)
	$(CC) generate_ionoscreen.c -o generate_ionoscreen $(OBJS) $(FLAGS) $(NETCDF_INCLUDE) -lm $(NETCDF_LIB) 

install: generate_ionoscreen
	chmod 0775 generate_ionoscreen
	mv generate_ionoscreen $(SIM_BIN)

clean:
	rm -f *.o generate_ionoscreen
