IPIC_HOME   = $(dir $(abspath $(lastword $(MAKEFILE_LIST))))

## SECTION THAT DEPENDS ON THE SYSTEM
#  -- Modify this variables according
#  -- to your system options.
#  -- Possible flags:
#     -DPARALLEL_IO  Use parallel HDF5
#     -DUSEH5HUT Use H5hut (must include also -DPARALLEL_IO)
#     -DBATSRUS Coupling with BATS-R-US
CXX         = mpicxx
HDF5_HOME   = /usr/local/hdf5/1.8.11-par
H5HUT_HOME  = /usr/local/H5hut/1.99.12
IPIC_FLAGS  = "-DUSEH5HUT -DPARALLEL_IO"
## END OF SECTION

OPTIM       = -O3

HDF5_LIB    = $(HDF5_HOME)/lib/libhdf5_hl.a $(HDF5_HOME)/lib/libhdf5.a
H5HUT_LIB   = $(H5HUT_HOME)/lib/libH5hut.a
H5HUTIO_LIB = $(IPIC_HOME)/H5hut-io/libH5hut-io.a

INC_DIR     = ./include
INC_H5HUT   = $(H5HUT_HOME)/include
INC_H5HUTIO = $(IPIC_HOME)/H5hut-io/include
INC_HDF5    = $(HDF5_HOME)/include

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

ALLOBJ = $(subst .cpp,.o,$(SRC))

IPIC3D_EXE = $(IPIC_HOME)/iPic3D
IPIC3D_LIB = $(IPIC_HOME)/libiPic3Dlib.a
LDLIBS = $(IPIC3D_LIB) $(H5HUTIO_LIB) $(H5HUT_LIB) $(HDF5_LIB) -ldl

all : io lib main

io :
	CXX=$(CXX) HDF5_HOME=$(HDF5_HOME) H5HUT_HOME=$(H5HUT_HOME) IPIC_FLAGS=$(IPIC_FLAGS) $(MAKE) -C $(IPIC_HOME)/H5hut-io

lib : $(ALLOBJ)
	$(AR) sr $(IPIC3D_LIB) $(ALLOBJ)
	ranlib $(IPIC3D_LIB)

main : iPic3D.o
	$(CXX) $(LDFLAGS) -I$(INC_DIR) -I$(INC_HDF5) iPic3D.cpp -o $(IPIC3D_EXE) $(LDLIBS)

clean : cleanio
	$(RM) $(ALLOBJ)
	$(RM) $(IPIC3D_LIB)
	$(RM) iPic3D.o
	$(RM) $(IPIC3D_EXE)

cleanio :
	$(MAKE) -C $(IPIC_HOME)/H5hut-io clean

%.o : %.cpp
	echo " Compiling " $@
	$(CXX) $(CXXFLAGS) $(OPTIM) $(IPIC_FLAGS) -I$(INC_DIR) -I$(INC_H5HUTIO) -I$(INC_H5HUT) -I$(INC_HDF5) -c $< -o $@

