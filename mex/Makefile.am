AM_CFLAGS=-DDLONG $(MEXOPT) -I$(SRC_DIR)/lib
AM_LDFLAGS=$(LDMEX) -lmwblas -lmwlapack -lm -lfftw3
if USE_LONG
AM_LDFLAGS+=$(BUILD_DIR)/lib/libaos_long.la
else
AM_LDFLAGS+=$(BUILD_DIR)/lib/libaos_int.la
endif
AM_LDFLAGS+=-Xcompiler -shared -no-install
EXTRA_DIST=null.c io.c io.h interface.h Makefile.mex
EXEEXT=.$(MEXEXT) #Automatic add suffix

.NOTPARALLEL:io.c

bin_PROGRAMS=read write mkgmex crandcreate crandn crandp clrand mkwmex proj genotfmex svdmex

#.c.o:
#	$(MEX) $(AM_CFLAGS) -c $< 

#Two special ones that does not depend on aos_long
#read$(EXEEXT):read.c io.c
#	$(MEX) $(SRC_DIR)/mex/read.c $(SRC_DIR)/mex/io.c -lz
#write$(EXEEXT):write.c io.c
#	$(MEX) $(SRC_DIR)/mex/write.c $(SRC_DIR)/mex/io.c -lz

read_SOURCES=read.c io.c
read_CFLAGS=$(MEXOPT)
read_LDFLAGS=$(LDMEX) -Xcompiler -shared

write_SOURCES=write.c io.c
write_CFLAGS=$(MEXOPT)
write_LDFLAGS=$(LDMEX) -Xcompiler -shared
