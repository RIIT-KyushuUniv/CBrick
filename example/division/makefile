SHELL=/bin/bash

CSRCS = div.cpp
FSRCS =
SRCS  = $(CSRCS) $(FSRCS)

.SUFFIXES: .o .f90 .cpp
COBJS = $(CSRCS:.cpp=.o)
FOBJS = $(FSRCS:.f90=.o)
OBJS  = $(COBJS) $(FOBJS)

TARGET = div

CC = icpc
#CFLAGS = -O3 -D_REAL_IS+_DOUBLE_
CFLAGS = -O3 -I../../src -DDebug
CLIBS = -L../../src -lCBrick
#LDFLAGS=
LDFLAGS

F90 =
F90FLAGS = -O3 -ffree-form -fopenmp -Wall
#F90FLAGS = -O3 -ffree-form -fdefault-real-8
F90LIBS =

RM = rm

$(TARGET):		$(OBJS)
	$(CC) $(CFLAGS) -o $(@) $(OBJS) $(F90LIBS) $(CLIBS) $(LDFLAGS)

.cpp.o:
	$(CC) $(CFLAGS) -c $<

.f90.o:
	$(F90) $(F90FLAGS) -c $<

clean:
	$(RM) $(OBJS) $(TARGET)

depend:
	makedepend -Y *.c *.h
# DO NOT DELETE
