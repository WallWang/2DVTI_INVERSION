LIB = ../lib/libstereo.a

ARCH =		\
	$(LIB)(spqr_solve_matrix.o)	\
	$(LIB)(sparse_matrix.o)	\
	$(LIB)(regu.o)	\
	$(LIB)(lsqr.o)	\
	$(LIB)(bsplines.o)	\
	$(LIB)(forward_modeling.o)	\
	$(LIB)(frechet.o)	\
	$(LIB)(data_fitting_err.o)	\
	$(LIB)(update_model.o)	\
	$(LIB)(stereo_subroutines.o)	

CC=mpicc
CFLAGS=	\
	-I$(CWPROOT)/include	\
	-I$(MPIROOT)/include    \
	-I../../include        \
        -I$(SUITROOT)/include      \
	-I../include -g

INSTALL: $(LIB)

$(LIB): $(ARCH)

clean:
	@-rm -f $(LIB)
