#Makefile for stereotomography

CFLAGS=                 \
	-I../include	\
	-I$(CWPROOT)/include        \
        -I./include     \
	-I$(SUITROOT)/include	\
	-L../lib       \
	-L/usr/local/gfortran/lib	\
	-L$(SUITROOT)/lib	\
	-L$(CWPROOT)/lib       \
	-L./lib
LFLAGS=                 \
	./lib/libstereo.a	\
	-lspqr -lsuitesparseconfig -lcholmod -lamd -lcolamd	\
	-lm -llapack -lblas -ltbb	\
	-lwx -lsub -lsu -lpar -lcwp -lstdc++

r:
	cd ./common; $(MAKE)
	cd ./sub; $(MAKE)
	cd ./src; $(MAKE)
	mpicc -std=c99 -O0 	\
	./main.c  \
	$(CFLAGS)	\
	-o ../bin/stereo2d	\
	$(LFLAGS) -g #-Wall
