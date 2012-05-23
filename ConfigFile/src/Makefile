# Makefile for ConfigFile class

CC = xlC
CFLAGS = -qsuppress=1501-245
#CC = g++
#CFLAGS = -Wall -ansi

all: configfile

configfile: 
	$(CC) $(CFLAGS) -c ConfigFile.cpp
tidy:
	@ rm -f *.o

clean: tidy
	@ rm -f example tester exampleparsek example.out test.out exampleparsek.out core*
