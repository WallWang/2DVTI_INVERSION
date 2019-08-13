include ../Makefile.config

LIB = ../lib/libsub.a

ARCH =		\
	$(LIB)(RaySub.o)

CFLAGS=-I$(CWPROOT)/include -I../include -L../lib -lwyx -O0 -g -Wall

INSTALL: $(LIB)

$(LIB): $(ARCH)

clean:
	@-rm -f $(LIB)
