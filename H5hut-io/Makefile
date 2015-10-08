
INC_DIR   = ./include
INC_HDF5  = $(HDF5_HOME)/include
INC_H5HUT = $(H5HUT_HOME)/include

SRC      = src/H5hut-io.cpp
ALLOBJ   = $(subst .cpp,.o,$(SRC))

H5HUTIO_LIB = libH5hut-io.a

main : $(ALLOBJ)
	$(AR) sr $(H5HUTIO_LIB) $(ALLOBJ)

clean :
	$(RM) $(ALLOBJ)
	$(RM) $(H5HUTIO_LIB)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(IPIC_FLAGS) -I$(INC_DIR) -I$(INC_H5HUT) -I$(INC_HDF5) -c $< -o $@
