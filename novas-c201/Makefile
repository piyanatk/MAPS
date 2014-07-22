# Makefile for NOVAS-C source code. Makes a library and copies the
# library and header files into appropriate MAPS dirs
CFLAGS = -O
OBJS = novascon.o novas.o solsys3.o readeph0.o
NOVASLIB = libnovas_c.a

install: lib
	mv $(NOVASLIB) $(SIM_LIB)
	cp *.h $(SIM_INC)

lib: $(OBJS)
	$(AR) $(ARFLAGS) $(NOVASLIB) $(OBJS)

check: $(OBJS)
	$(CC) -o checkout checkout-st.c $(OBJS) -lm

clean:
	rm -f *.o checkout $(NOVASLIB)
	rm -f $(SIM_INC)/novas.h
	rm -f $(SIM_INC)/novascon.h
	rm -f $(SIM_INC)/solarsystem.h
