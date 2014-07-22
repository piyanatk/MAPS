FC = g77

LIBRARY  = libstar.a

SOURCES = $(wildcard *.f) 

# C wrappers to the Fortran SLALIB bindings
HEADER = slalib.h
SOURCES += sla_wrap.c
BASENAMES = $(basename $(SOURCES))
OBJS = $(addsuffix .o, $(BASENAMES));

all: $(LIBRARY) $(HEADER)

%.o:%.f
	$(FC) $(FFLAGS) -c $<
%.o:%.c
	$(CC) $(CFLAGS) -c $<

$(LIBRARY): $(OBJS) $(HEADERS)
	$(AR) $(ARFLAGS) $(LIBRARY) $(OBJS)

install: all
	cp $(HEADER) $(SIM_INC)
	cp $(LIBRARY) $(SIM_LIB)

clean:
	rm -f *.o
	rm -f libstar.a
	rm -f $(SIM_LIB)/$(LIBRARY) 
	rm -f $(SIM_INC)/$(HEADER)

