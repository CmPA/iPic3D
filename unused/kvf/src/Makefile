
#CPL=pathCC
#CC=pathcc
CPL=g++
CC=gcc

CFLAGS=-v
CPLFLAGS=-v
.cc.o:
	$(CPL) $(CFLAGS)  $(PIC) -c $*.cc
.c.o:
	$(CC) $(CPLFLAGS)  $(PIC) -c $*.c

KVF_HDRS= \
  kvfDataBuffer.h \
  kvfDataSource.h \
  KVFparser.h \
  KVFyacc.h \
  kvfNaming.h \
  kvfVar.h 

MDOT_OBJS=\
  kvfDataBuffer.o \
  kvfDataSource.o \
  kvfNaming.o \
  KVFscanner.o \
  KVFyacc.o 
  
all: libkvf.a

libkvf.a: $(MDOT_HDRS) $(MDOT_OBJS)
	ar ru libkvf.a $(MDOT_OBJS)
	ranlib libkvf.a

clean:
	'rm' -f $(MDOT_OBJS) libkvf.a

install_mdot:
	cp -f $(MDOT_HDRS) ../include
	cp libkvf.a ../lib

KVFyacc.c: KVFyacc.y KVFparser.h
	bison -d -oKVFyacc.c -pKVFyy KVFyacc.y
KVFyacc.h: KVFyacc.c
KVFyacc.o: KVFyacc.c KVFparser.h
	$(CPL) -c $(CPLFLAGS)  $(PIC) KVFyacc.c

KVFscanner.cc: KVFscanner.l KVFyacc.h kvfDataSource.h\
               kvfNaming.h kvfVar.h kvfDataBuffer.h kvfDataSource.h
	flex -oKVFscanner.cc -PKVFyy KVFscanner.l
KVFscanner.o: KVFscanner.l KVFyacc.y KVFparser.h 


#new: $(QDOS_OBJS) install libqdos.so
#	$(INSTALL)  libqdos.so $(QS_HOME)/lib 

libqdos.so: $(QDOS_OBJS) 
	$(LD) $(PIC_LD) $(PIC) -o libqdos.so $(QDOS_OBJS) 
        

kvfDataBuffer.o: kvfDataBuffer.h kvfException.h kvfVar.h kvfNaming.h
kvfDataSource.o: kvfDataSource.h kvfException.h \
                 kvfDataBuffer.h kvfVar.h kvfNaming.h

install:
	$(INSTALL) Adapter.h $(QS_HOME)/include
