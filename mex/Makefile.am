AM_CFLAGS=-DDLONG $(MEXOPT) -I$(SRC_DIR)/lib -I$(SRC_DIR)/math
AM_LDFLAGS= $(LDMEX) -lm -static  -lfftw3 -lmwblas -lmwlapack
AM_LDFLAGS+=-Xcompiler -shared -no-fast-install
#The following should be in LDADD instead of AM_LDFLAGS to recompile exes when they are updated
LDADD=$(BUILD_DIR)/lib/libaos.la $(BUILD_DIR)/math/libaomath.la $(BUILD_DIR)/sys/libaosys.la
EXTRA_DIST=null.c io.c io.h interface.h  Makefile.mex compile.sh
EXEEXT=.$(MEXEXT) #Automatic add suffix

.NOTPARALLEL:io.c

bin_PROGRAMS=read write mkgmex mkhmex crandcreate crandn crandp clrand mkwmex proj genotfmex svdmex mk2dcov servo_optim_mex servo_sim_mex psd1d mkcirmap sde_fit_mex sde_kalman_mex kalman_sim_mex calcenc act_extrap cn2est

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
