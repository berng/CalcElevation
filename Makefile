CC = gcc -DKERFIX
CFLAGS = -Iinclude/superdarn -Iinclude/base -Iinclude/general -I/usr/include -Iusr/include/superdarn -I./src
LDFLAGS = -lm -Llib -lz -lfit.1 -ldmap.1 -laacgm.1 -lradar.1 -lrcnv.1 -lcfit.1 -lraw.1 -lrscan.1 -lrtime.1
FFLAGS = -freal-4-real-8
#LDFLAGS = -lm -Llib -lz

all:	recalc_elevation_DAVIT_linear


recalc_elevation_DAVIT_linear:	recalc_elevation_DAVIT_linear.c
	$(CC) -o recalc_elevation_DAVIT_linear recalc_elevation_DAVIT_linear.c $(CFLAGS) $(LDFLAGS) 

distclean:
	rm -f *.o

clean:
	rm -f recalc_elevation_DAVIT_linear
