AM_CFLAGS=-DDLONG $(MEXOPT) -I$(SRC_DIR)/lib
AM_LDFLAGS=-lfftw3 $(LDMEX) -lmwblas -lmwlapack -lm -static
AM_LDFLAGS+=-Xcompiler -shared -no-fast-install
#The following should be in LDADD instead of AM_LDFLAGS to recompile exes when they are updated
if USE_LONG
LDADD=$(BUILD_DIR)/lib/libaos_long.la $(BUILD_DIR)/sys/libaosys.la 
else
LDADD=$(BUILD_DIR)/lib/libaos_int.la $(BUILD_DIR)/sys/libaosys.la 
endif
EXTRA_DIST=null.c io.c io.h interface.h  Makefile.mex
EXEEXT=.$(MEXEXT) #Automatic add suffix

.NOTPARALLEL:io.c

bin_PROGRAMS=read write mkgmex crandcreate crandn crandp clrand mkwmex proj genotfmex svdmex mk2dcov servo_optim_mex servo_sim_mex psd1d mkcirmap

#.c.o:
#	$(MEX) $(AM_CFLAGS) -c $< 

#Two special ones that does not depend on aos_long
#read$(EXEEXT):read.c io.c
#	$(MEX) $(SRC_DIR)/mex/read.c $(SRC_DIR)/mex/io.c -lz
#write$(EXEEXT):write.c io.c
#	$(MEX) $(SRC_DIR)/mex/write.c $(SRC_DIR)/mex/io.c -lz

read_SOURCES=read.c io.c
read_CFLAGS=$(MEXOPT)
read_LDFLAGS=$(LDMEX) -static -Xcompiler -shared -no-fast-install
read_LDADD=
write_SOURCES=write.c io.c
write_CFLAGS=$(MEXOPT)
write_LDFLAGS=$(LDMEX) -static -Xcompiler -shared -no-fast-install
write_LDADD=
