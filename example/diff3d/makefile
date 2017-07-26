SHELL=/bin/bash

CSRCS = utility.c diff3dp.c table.c comm.c
FSRCS = fcore.f90
SRCS  = $(CSRCS) $(FSRCS)

.SUFFIXES: .o .f90 .c
COBJS = $(CSRCS:.c=.o)
FOBJS = $(FSRCS:.f90=.o)
OBJS  = $(COBJS) $(FOBJS)

TARGET = diff3dhp

CC = mpicc
#CFLAGS = -O3 -D_REAL_IS+_DOUBLE_
CFLAGS = -O3 -fopenmp
CLIBS =
#LDFLAGS=
LDFLAGS= -L/usr/local/gfortran/lib -lgfortran

F90 = mpif90
F90FLAGS = -O3 -ffree-form -fopenmp -Wall
#F90FLAGS = -O3 -ffree-form -fdefault-real-8
F90LIBS =

RM = rm

$(TARGET):		$(OBJS)
	$(CC) $(CFLAGS) -o $(@) $(OBJS) $(F90LIBS) $(CLIBS) $(LDFLAGS)

.c.o:
	$(CC) $(CFLAGS) -c $<

.f90.o:
	$(F90) $(F90FLAGS) -c $<

clean:
	$(RM) $(OBJS) $(TARGET)

depend:
	makedepend -Y *.c *.h
# DO NOT DELETE
