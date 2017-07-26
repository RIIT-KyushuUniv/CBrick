SHELL=/bin/bash

TARGET = libCBrick.a
SRCS = CB_SubDomain.cpp

.SUFFIXES: .o .cpp
OBJS = $(SRCS:.cpp=.o)
all : $(TARGET)

#MPIINC  = -I/usr/local/mpich/include
#INCPATH = -I../../include -I../../include/sphmpi $(MPIINC)

CXX      = icpc
CXXFLAGS = -Wall -g -DDebug $(INCPATH)
RANLIB   = echo
MV       = \mv

$(TARGET): $(OBJS)
	ar -rv $(TARGET) $(OBJS)
	$(RANLIB) $(TARGET)

.cpp.o:
	$(CXX) $(CXXFLAGS) -c $<

clean:
	$(RM) *.o

depend:
	makedepend -I. $(INCPATH) -Y *.cpp

# DO NOT DELETE

CB_DomainManager.o: CB_DomainManager.h CB_Define.h
CB_SubDomain.o: CB_SubDomain.h CB_Define.h
