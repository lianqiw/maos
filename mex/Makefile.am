#To use this make file. Please prepare as follows
# 1) Compile maos with configure --enable-shared --enable-long --disable-icc
# 2) Copy lib/.libs/libaos_long.* to /usr/local/lib or other place in or LDPATH
# 3) make -f Makefile.full

#.NOTPARALLEL

mexopts=$(shell mex -h|grep -q largeArrayDims && echo -largeArrayDims)
mexopts+=-DDLONG -D_BSD_SOURCE -D_POSIX_C_SOURCE=200112L -I$(SRC_DIR)/lib/
AM_LDFLAGS=-DDLONG -lz -lmwblas -lmwlapack -lm -largeArrayDims
if USE_LONG
AM_LDFLAGS+=../lib/.libs/libaos_long.so 
else
AM_LDFLAGS+=../lib/.libs/libaos_int.so 
endif

MEX=@mex $(mexopts) 
CCLD=mex
EXTRA_DIST=null.c io.c io.h interface.h
EXEEXT=$(MEXEXT) #Automatic add suffix

bin_PROGRAMS=read write mkgmex crandcreate crandn crandp clrand mkwmex proj genotfmex svdmex

.c.o:
	$(MEX) $(AM_CFLAGS) -c $< 

#Two special ones that does not depend on aos_long
read$(EXEEXT):read.c io.c
	$(MEX) $(SRC_DIR)/mex/read.c $(SRC_DIR)/mex/io.c -lz
write$(EXEEXT):write.c io.c
	$(MEX) $(SRC_DIR)/mex/write.c $(SRC_DIR)/mex/io.c -lz

