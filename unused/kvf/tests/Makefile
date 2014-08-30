
CPL=g++

CPLFLAGS=-v -I"../src"

.cc.o:
	$(CPL)  $(CPLFLAGS)  -c $*.cc

tKVFDataSource: tKVFDataSource.o
	$(CPL) -v -o tKVFDataSource tKVFDataSource.o  -L../src -lkvf



