all: semiconductor fvsg nano nano_test

clean:
	/bin/rm -f *.o

#include /opt/phg/phg-0.9.3/share/phg/Makefile.inc
PHG_MAKEFILE_INC = /share/soft/phg/phg-0.9.4-mvapich2-20190318/share/phg/Makefile.inc
include ${PHG_MAKEFILE_INC}

semiconductor.o: semiconductor.c

fvsg.o:	fvsg.c	

#nano.o:	/opt/phg/phg-0.9.3/lib/libphg${LIB_SUFFIX} nano.c PNP_analytic.c

nano.o: PNP_coefficient.h nano.c PNP_analytic.c 

nano_test.o: PNP_coefficient.h nano_test.c PNP_analytic.c 

runmp:
	mpirun -np 4 ./semiconductor ./fvsg ./nano
