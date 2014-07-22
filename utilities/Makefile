# Makefile for libutil                   CJL Jan 31 2002
#

RESULT=libutil.a

.c.o:
	cc -g -O -D_FILE_OFFSET_BITS=64 -c -I$(SIM_INC) $(RTS_INCS) $<
# 	cc -c $(SIM_CFLAGC) -I$(SIM_INC) $<

SOURCES=\
clear_date.c  double_to_time.c  report_times.c\
account.c     confirm.c     environment.c     msg.c     time_to_double.c\
rm_whitespace.c clear_coord.c comp_func.c

OBJECTS=\
clear_date.o  double_to_time.o   report_times.o\
account.o     confirm.o     environment.o     msg.o     time_to_double.o\
rm_whitespace.o clear_coord.o comp_func.o

# Jan 2009. SLALIB now obsolete and unsupported. Use novas instead
STAR=star/libstar.a

libutil.a: $(OBJECTS)
	ar cvr libutil.a $(OBJECTS)
	ranlib libutil.a

$(STAR):
	cd star;$(MAKE) install;
	
install: libutil.a 
	chmod 0664 $(RESULT)
# 	chgrp vlbi $(RESULT)
	cp $(RESULT) $(SIM_LIB)

clean:
	rm -f *.o libutil.a
	cd star;$(MAKE) clean

