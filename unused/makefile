#DYLD_LIBRARY_PATH = /users/cpa/volshevs/local/hdf5/lib

CPP = mpicxx
OPTFLAGS=  -O2 #-DMPICH_IGNORE_CXX_SEEK

HDF5LIBS = -lhdf5 -lhdf5_hl 

ipic3D: iPIC3D.cpp Particles3Dcomm.o Particles3D.o ConfigFile.o
	${CPP} ${OPTFLAGS} -o iPIC3D \
	iPIC3D.cpp Particles3Dcomm.o Particles3D.o ConfigFile.o \
        ${HDF5LIBS}

iPIC3D.o: iPIC3D.cpp
	${CPP} ${OPTFLAGS} -c iPIC3D.cpp 

ConfigFile.o: ./ConfigFile/src/ConfigFile.cpp
	${CPP} ${OPTFLAGS} -c ./ConfigFile/src/ConfigFile.cpp

Particles3Dcomm.o: ./particles/Particles3Dcomm.cpp
	${CPP} ${OPTFLAGS} -c ./particles/Particles3Dcomm.cpp

Particles3D.o: ./particles/Particles3D.cpp 
	${CPP} ${OPTFLAGS} -c ./particles/Particles3D.cpp

clean:
	rm -rf *.o iPIC3D

run:
	mpiexec -n 2 ./iPIC3D inputfiles/GEM.inp

tags: retags

retags:
	rm -f tags
	find . -name '*.h' -or -name '*.cpp' | xargs ctags --extra=+q
	find . -name '*.h' -or -name '*.cpp' | xargs makefiletags >> tags
	LC_ALL=C sort -u tags -o tags

