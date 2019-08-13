include ../Makefile.config

LIB = ../lib/libwx.a

ARCH =		\
	$(LIB)(iofile.o)	\
	$(LIB)(wx.o)

CFLAGS=	\
	-I$(CWPROOT)/include	\
	-I../include	\
	-O0 -lstdc++

INSTALL: $(LIB)

$(LIB): $(ARCH)

clean:
	@-rm -f $(LIB)
