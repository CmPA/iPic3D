##
## USER-DEFINED SECTION
##

IPIC_HOME = $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
IPIC_EXE = $(IPIC_HOME)/iPic3D
IPIC_LIB = $(IPIC_HOME)/libiPic3Dlib.a

##
## SYSTEM-DEPENDENT SECTION
##

#  -- Modify these variables according to your system.
#  -- Possible "IPIC_FLAGS" are:
#     -DPARALLEL_IO  Use parallel HDF5
#     -DUSEH5HUT     Use H5hut (must include also -DPARALLEL_IO)
#     -DBATSRUS      Coupling with BATS-R-US

CXX         = mpicxx

# Uncomment this line if you want to use the parallel HFD5 library:
#IPIC_FLAGS  += -DPARALLEL_IO

# If you have installed HDF5 manually in a non-default location, please
# uncomment this block and set the location in "HDF5_HOME", if necessary:
#HDF5_HOME   += /usr/local/hdf5/1.8.11-par
#HDF5_INCS   += -I$(HDF5_HOME)/include
#HDF5_LIBS   += -L$(HDF5_HOME)/lib

# Uncomment these lines if you want to use the H5HUT library:
#IPIC_FLAGS  += -DUSEH5HUT
#H5HUT_LIBS  += -lH5hut

# If you have installed H5HUT manually in a non-default location, please
# uncomment this block and set the location in "H5HUT_HOME", if necessary:
#H5HUT_HOME  += /usr/local/H5hut/1.99.12
#H5HUT_INCS  += -I$(H5HUT_HOME)/include
#H5HUT_LIBS  += -L$(H5HUT_HOME)/lib

##
## GENERIC SECTION (do not edit, unless you know what you do!)
##

# Some HDF5 library is mandatory, be it the sequential or parallel one:
HDF5_LIBS   += -lhdf5 -lhdf5_hl

# Optimization level:
OPTIM       += -O3
OPTIM       += -march=native

# Debugging options:
DEBUG       += -g3 -ggdb

# Include directories:
INC_DIRS    += -I./include
INC_DIRS    += $(HDF5_INCS)
INC_DIRS    += $(H5HUT_INCS)

# Additional libararies:
LD_LIBS     += -ldl
LD_LIBS     += $(HDF5_LIBS)
LD_LIBS     += $(H5HUT_LIBS)

.SUFFIXES:
.SUFFIXES: .cpp .o .h

SRC = \
	$(IPIC_HOME)/grids/Grid3DCU.cpp \
	$(IPIC_HOME)/fields/BCStructure.cpp \
	$(IPIC_HOME)/fields/EMfields3D.cpp \
	$(IPIC_HOME)/inputoutput/phdf5.cpp \
	$(IPIC_HOME)/inputoutput/Restart3D.cpp \
	$(IPIC_HOME)/inputoutput/ParallelIO.cpp \
	$(IPIC_HOME)/inputoutput/Collective.cpp \
	$(IPIC_HOME)/performances/Timing.cpp \
	$(IPIC_HOME)/PSKOutput3D/PSKhdf5adaptor.cpp \
	$(IPIC_HOME)/bc/BcParticles.cpp \
	$(IPIC_HOME)/bc/BcFields3D.cpp \
	$(IPIC_HOME)/mathlib/EllipticF.cpp \
	$(IPIC_HOME)/solvers/CG.cpp \
	$(IPIC_HOME)/solvers/GMRES.cpp \
	$(IPIC_HOME)/ConfigFile/src/ConfigFile.cpp \
	$(IPIC_HOME)/main/iPic3Dlib.cpp \
	$(IPIC_HOME)/particles/Particles3Dcomm.cpp \
	$(IPIC_HOME)/particles/Particles3D.cpp \
	$(IPIC_HOME)/communication/ComNodes3D.cpp \
	$(IPIC_HOME)/communication/ComParser3D.cpp \
	$(IPIC_HOME)/communication/ComInterpNodes3D.cpp \
	$(IPIC_HOME)/communication/ComParticles3D.cpp \
	$(IPIC_HOME)/communication/ComBasic3D.cpp \
	$(IPIC_HOME)/communication/VCtopology3D.cpp \
	$(IPIC_HOME)/utility/debug.cpp \
	$(IPIC_HOME)/utility/asserts.cpp

ALL_OBJS = $(subst .cpp,.o,$(SRC))

all : io lib main

io :
	CXX=$(CXX) HDF5_HOME=$(HDF5_HOME) H5HUT_HOME=$(H5HUT_HOME) IPIC_FLAGS="$(IPIC_FLAGS)" $(MAKE) -C $(IPIC_HOME)/H5hut-io

lib : $(ALL_OBJS)
	$(AR) sr $(IPIC_LIB) $(ALL_OBJS)
	ranlib $(IPIC_LIB)

main : lib iPic3D.o
	$(CXX) iPic3D.o -o $(IPIC_EXE) $(IPIC_LIB) $(LD_LIBS)

clean : cleanio
	$(RM) $(ALL_OBJS)
	$(RM) $(IPIC_LIB)
	$(RM) $(IPIC_EXE).o
	$(RM) $(IPIC_EXE)

cleanio :
	$(MAKE) -C $(IPIC_HOME)/H5hut-io clean

%.o : %.cpp
	$(CXX) $(OPTIM) $(DEBUG) $(IPIC_FLAGS) $(INC_DIRS) -c $< -o $@

