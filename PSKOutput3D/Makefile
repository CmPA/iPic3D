
HDF5DIR = /home/network/ec/hdf5-1.6.5
HDF5INCDIR = $(HDF5DIR)/include
HDF5LIBDIR = $(HDF5DIR)/lib
HDF5LIBFLAGS =  -L$(HDF5LIBDIR) -lhdf5 -lhdf5_hl 



CPLFLAGS=-v -I$(HDF5INCDIR) -L$(HDF5LIBDIR) 
CFLAGS=-v -I$(HDF5INCDIR)

CPL=g++
CC=gcc


.cc.o:
	$(CPL) $(CPLFLAGS)  $(PIC) -c $*.cc
.c.o:
	$(CC) $(CFLAGS)  $(PIC) -c $*.c

test_pskoutput: test_pskoutput.o 
	$(CPL) -o test_pskoutput test_pskoutput.o $(HDF5LIBFLAGS)

testhdf5: testhdf5.o
	$(CPL) -o testhdf5 testhdf5.o $(HDF5LIBFLAGS)

testhdf5.o: testhdf5.cc PSKException.h  PSKhdf5adaptor.h  PSKOutput.h

test_pskoutput.o: test_pskoutput.cc PSKException.h \
                  PSKhdf5adaptor.h  PSKOutput.h
