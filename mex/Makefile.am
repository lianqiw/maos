AM_CFLAGS=$(MEXOPT) -I$(SRC_DIR)/lib -I$(SRC_DIR)/math 
AM_LDFLAGS=$(LDMEX)
#LDFLAGS goes before LDADD. Must put maos libraries before matlab libraries.
#The following should be in LDADD instead of AM_LDFLAGS to recompile exes when they are updated
LDADD=../lib/libaos.la -lcholmod -lfftw3 $(LIBMEX) 
EXTRA_DIST=null.c io.c io.h interface.h Makefile.mex compile.sh
EXEEXT=.$(MEXEXT) #Automatic add suffix

.NOTPARALLEL:io.c

bin_PROGRAMS=read write mkgmex mkhmex crandcreate crandn crandp clrand mkwmex proj genotfmex svdmex mk2dcov servo_optim_mex servo_sim_mex psd1d mkcirmap sde_fit_mex sde_kalman_mex kalman_sim_mex calcenc act_extrap cn2est maos

#.c.o:
#	$(MEX) $(AM_CFLAGS) -c $< 

#Two special ones that does not depend on aos_long
#read$(EXEEXT):read.c io.c
#	$(MEX) $(SRC_DIR)/mex/read.c $(SRC_DIR)/mex/io.c -lz
#write$(EXEEXT):write.c io.c
#	$(MEX) $(SRC_DIR)/mex/write.c $(SRC_DIR)/mex/io.c -lz

read_SOURCES=read.c io.c
read_CFLAGS=$(MEXOPT)
read_LDADD=$(LIBMEX)
write_SOURCES=write.c io.c
write_CFLAGS=$(MEXOPT)
write_LDADD=$(LIBMEX)

maos_LDADD=../maos/libmaos.la 
if USE_CUDA
maos_LDADD+=$(LIBCUDA) #must be before LIBMEX, otherwise linking error in debug mode
endif
maos_LDADD+=-lcholmod -lfftw3 $(LIBMEX) 


BUILT_SOURCES=maos2mex.h
maos2mex.h:../maos/parms.h ../maos/types.h ../scripts/maos2mex.py
	$(SRC_DIR)/scripts/maos2mex.py $(SRC_DIR)/mex/maos2mex.h
