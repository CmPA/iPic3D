#makefile for PARSEK PROJECT
CPP=	mpicxx
OPTFLAGS= -DMPICH_IGNORE_CXX_SEEK -O2

INC_KVF = -I ${PWD}/kvf/src
INC_HDF5 = 
INC_MPI = 

LIB_KVF = -L${PWD}/kvf/src
LIB_HDF5 =  
LIB_MPI = 

HDF5LIBS = -lhdf5 -lhdf5_hl 
KVFLIB = -lkvf 
MPELIB = -lmpe

ipic3D: iPIC3D.cpp  iPIC3D.o Particles3Dcomm.o Particles3D.o
	${CPP}	${OPTFLAGS} -o  iPIC3D ${INC_HDF5} ${INC_MPI} ${INC_KVF}\
	iPIC3D.cpp iPIC3D.o Particles3Dcomm.o Particles3D.o ${LIB_KVF} ${LIB_HDF5} ${LIB_KVF} ${LIB_MPI} \
${HDF5LIBS}	${KVFLIB}	${MPELIB}

iPIC3D.o:	iPIC3D.cpp
	${CPP} ${OPTFLAGS} ${INC_HDF5} ${INC_MPI} ${INC_KVF} -c iPIC3D.cpp 

Particles3Dcomm.o:	./particles/Particles3Dcomm.cpp
	${CPP}	${OPTFLAGS}	${INC_HDF5} ${INC_KVF} ${LIB_KVF}  -c ./particles/Particles3Dcomm.cpp

Particles3D.o:	./particles/Particles3D.cpp 
	${CPP}	${OPTFLAGS} ${INC_HDF5} ${INC_KVF} ${LIB_KVF}  -c ./particles/Particles3D.cpp

clean:
	rm -rf *.o iPIC3D 
