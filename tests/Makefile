
INCLUDE = -I../include

OBJECTS = \
  test_arrays.o \
  ../utility/asserts.o \
  debug.o

FLAGS = -O3 -DNO_MPI -fno-exceptions #-DCHECK_BOUNDS -ggdb

COMPILER = g++ #icpc # g++

test: clean_test_arrays test_arrays

test_arrays: clean_test_arrays $(OBJECTS)
	$(COMPILER) $(FLAGS) $(INCLUDE) $(OBJECTS) -o test_arrays

test_arrays.o:
	$(COMPILER) $(FLAGS) -c test_arrays.cpp $(INCLUDE) -o test_arrays.o


../utility/asserts.o:
	$(COMPILER) $(FLAGS) -c ../utility/asserts.cpp $(INCLUDE) -o ../utility/asserts.o

debug.o:
	$(COMPILER) $(FLAGS) -c ../utility/debug.cpp $(INCLUDE) -o debug.o

clean: clean_old_test_arrays clean_test_arrays
	rm -f test_arrays $(OBJECTS)

clean_test_arrays:
	rm -f test_arrays test_arrays.o

old_test_arrays: clean_old_test_arrays
	$(COMPILER) $(FLAGS) old_test_arrays.cpp  -o old_test_arrays

clean_old_test_arrays:
	rm -f old_test_arrays old_test_arrays.o
