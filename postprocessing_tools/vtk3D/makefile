CPP=	g++
OPTFLAGS=  -O2

INC_HDF5 = -I/usr/local/hdf5/include
LIB_HDF5 = -L/usr/local/hdf5/lib
INC = -I../../include

HDF5LIBS =  -lhdf5 -lhdf5_hl  -lz

all: conv

conv: convHDF5unico.cpp
	${CPP}	${OPTFLAGS} -o  convHDF5unico ${INC_HDF5} ${INC} convHDF5unico.cpp  ${LIB_HDF5} ${HDF5LIBS}

clean:
	rm convHDF5unico
