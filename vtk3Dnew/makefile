#makefile for PARSEK PROJECT
CPP=	mpicxx
OPTFLAGS=  -O2 

SYSTEM = mac

ifeq ($(SYSTEM),pleiades)
CPP=icc
INC_HDF5 = -I/nasa/hdf5/1.8.3/intel/mpt/include
LIB_HDF5 = -L/nasa/hdf5/1.8.3/intel/mpt/lib
endif

ifeq ($(SYSTEM),mms)
# UBUNTU 
# include files
INC_HDF5 = -I/shared/hdf5-1.8.9/hdf5/include
# libs
LIB_HDF5 = -L/shared/hdf5-1.8.9/hdf5/lib
endif
ifeq ($(SYSTEM),ubuntu-home)
# UBUNTU 
# include files
INC_HDF5 = -I/home/gianni/Downloads/hdf5-1.6.10/hdf5/include
# libs
LIB_HDF5 = -L/home/gianni/Downloads/hdf5-1.6.10/hdf5/lib
endif

ifeq ($(SYSTEM),ubuntu)
# UBUNTU 
# include files
INC_HDF5 = -I/home/gianni/computer-stuff/hdf5-1.6.10/hdf5/include
# libs
LIB_HDF5 = -L/home/gianni/computer-stuff/hdf5-1.6.10/hdf5/lib
endif

ifeq ($(SYSTEM),mac)
#MAC
# include files
INC_HDF5 = -I/Users/gianni/Downloads/hdf5serial1.8.11/hdf5/include
# libs
LIB_HDF5 = -L/Users/gianni/Downloads/hdf5serial1.8.11/hdf5/lib 
endif

ifeq ($(SYSTEM),mac-uff)
#MAC
# include files
INC_HDF5 = -I/Users/giovannilapenta/Downloads/hdf5-1.8.14/hdf5/include
# libs
LIB_HDF5 = -L/Users/giovannilapenta/Downloads/hdf5-1.8.14/hdf5/lib
endif

ifeq ($(SYSTEM),curie)
# CURIE
# include files
INC_HDF5 = -I/ccc/cont005/home/ra0747/lapentg/hdf5-1.6.10/hdf5/include
# libs
LIB_HDF5 = -L/ccc/cont005/home/ra0747/lapentg/hdf5-1.6.10/hdf5/lib
MPELIB = -L/opt/mpi/bullxmpi/1.1.14.3/lib/
endif

ifeq ($(SYSTEM),buteo)
# CURIE
# include files
INC_HDF5 = -I/Users/Gianni/Downloads/hdf5-1.8.9/hdf5/include

# libs
LIB_HDF5 = -L/Users/Gianni/Downloads/hdf5-1.8.9/hdf5/lib
endif



HDF5LIBS =  -lhdf5 -lhdf5_hl  -llapack -lblas

xyzGDA_par: xyzGDA_par.cpp
	${CPP}  ${OPTFLAGS} -o  xyzGDA_par ${INC_HDF5} \
	xyzGDA_par.cpp  ${LIB_HDF5}  \
${HDF5LIBS}

xyzHDF5_par: xyzHDF5_par.cpp
	${CPP}  ${OPTFLAGS} -o  xyzHDF5_par ${INC_HDF5} \
	xyzHDF5_par.cpp  ${LIB_HDF5}  \
${HDF5LIBS}

read_ASCII: read_ASCII.cpp 
	${CPP}	${OPTFLAGS} -o  read_ASCII ${INC_HDF5} \
	read_ASCII.cpp  ${LIB_HDF5}  \
${HDF5LIBS}	

convHDF5unico: convHDF5unico.cpp 
	${CPP}	${OPTFLAGS} -o  convHDF5unico ${INC_HDF5} \
	convHDF5unico.cpp  ${LIB_HDF5}  \
${HDF5LIBS}	

ytHDF5: ytHDF5.cpp 
	${CPP}	${OPTFLAGS} -o  ytHDF5 ${INC_HDF5} \
	ytHDF5.cpp  ${LIB_HDF5}  \
${HDF5LIBS}	

xtHDF5: xtHDF5.cpp 
	${CPP}	${OPTFLAGS} -o  xtHDF5 ${INC_HDF5} \
	xtHDF5.cpp  ${LIB_HDF5}  \
${HDF5LIBS}	

xtGDA: xtGDA.cpp 
	${CPP}	${OPTFLAGS} -o  xtHDF5 ${INC_HDF5} \
	xtGDA.cpp  ${LIB_HDF5}  \
${HDF5LIBS}	

xtGDA_pa: xtGDA_par.cpp 
	${CPP}	${OPTFLAGS} -o  xtHDF5_par ${INC_HDF5} \
	xtGDA_par.cpp  ${LIB_HDF5}  \
${HDF5LIBS}	

readVTK: read_vtk.cpp 
	${CPP}	${OPTFLAGS} -o  readVTK  ${INC_HDF5} \
	read_vtk.cpp  ${LIB_HDF5}  \
${HDF5LIBS}	

reconHDF5: reconHDF5.cpp 
	${CPP}	${OPTFLAGS} -o  reconHDF5  ${INC_HDF5} \
	reconHDF5.cpp  ${LIB_HDF5}  \
${HDF5LIBS}	

xyzHDF5: xyzHDF5.cpp 
	${CPP}	${OPTFLAGS} -o  xyzHDF5 ${INC_HDF5} \
	xyzHDF5.cpp  ${LIB_HDF5}  \
${HDF5LIBS}	

xyzHDF5_3sp: xyzHDF5_3sp.cpp 
	${CPP}	${OPTFLAGS} -o  xyzHDF5_3sp ${INC_HDF5} \
	xyzHDF5_3sp.cpp  ${LIB_HDF5}  \
${HDF5LIBS}	

AVGxyHDF5: AVGxyHDF5.cpp 
	${CPP}	${OPTFLAGS} -o  AVGxyHDF5 ${INC_HDF5} \
	AVGxyHDF5.cpp  ${LIB_HDF5}  \
${HDF5LIBS}	

AVGxyPHDF5: AVGxyPHDF5.cpp 
	${CPP}	${OPTFLAGS} -o  AVGxyPHDF5 ${INC_HDF5} \
	AVGxyPHDF5.cpp  ${LIB_HDF5}  \
${HDF5LIBS}	

xyzGDA: xyzGDA.cpp
	${CPP}  ${OPTFLAGS} -o  xyzGDA ${INC_HDF5} \
	xyzGDA.cpp  ${LIB_HDF5}  \
${HDF5LIBS}

cut_xzHDF5: cut_xzHDF5.cpp 
	${CPP}	${OPTFLAGS} -o  cut_xzHDF5 ${INC_HDF5} \
	cut_xzHDF5.cpp  ${LIB_HDF5}  \
${HDF5LIBS}	

partHDF5: partHDF5.cpp 
	${CPP}	${OPTFLAGS} -o  partHDF5 ${INC_HDF5} \
	partHDF5.cpp  ${LIB_HDF5}  \
${HDF5LIBS}	

presHDF5: presHDF5.cpp
	${CPP}  ${OPTFLAGS} -o  presHDF5 ${INC_HDF5} \
	presHDF5.cpp  ${LIB_HDF5}  \
${HDF5LIBS}

zenitani: zenitani.cpp
	${CPP}  ${OPTFLAGS} -o  zenitani ${INC_HDF5} \
	zenitani.cpp  ${LIB_HDF5}  \
${HDF5LIBS}

clean:
	rm -rf convHDF5unico
