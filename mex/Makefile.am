AM_CFLAGS=$(MEXOPT) -I$(SRC_DIR)/lib -I$(SRC_DIR)/math 
AM_LDFLAGS=$(LDMEX) $(LDSO)
#LDFLAGS goes before LDADD. Must put maos libraries before matlab libraries.
#The following should be in LDADD instead of AM_LDFLAGS to recompile exes when they are updated
LDADD=../lib/libaos.la $(LIBMEX) 
EXTRA_DIST=io.h interface.h Makefile.mex compile.sh aolib.h 
EXEEXT=.$(MEXEXT) #Automatic add suffix

.NOTPARALLEL:io.c
if USE_MEX
bin_PROGRAMS=read write aolib maos 

#.c.o:
#	$(MEX) $(AM_CFLAGS) -c $< 

#Two special ones that does not depend on aos library
#read$(EXEEXT):read.c io.c
#	$(MEX) $(SRC_DIR)/mex/read.c $(SRC_DIR)/mex/io.c -lz
#write$(EXEEXT):write.c io.c
#	$(MEX) $(SRC_DIR)/mex/write.c $(SRC_DIR)/mex/io.c -lz

read_SOURCES=read.c io.c
read_CFLAGS=$(MEXOPT)
read_LDADD=-lmx -lmex -lz

write_SOURCES=write.c io.c
write_CFLAGS=$(MEXOPT)
write_LDADD=-lmx -lmex -lz

maos_LDADD=../maos/libmaos.la 
if USE_CUDA
maos_LDADD+=$(LIBCUDA) #must be before LIBMEX, otherwise linking error in debug mode
endif
maos_LDADD+=-lmwcholmod $(LIBMEX) 
endif

BUILT_SOURCES=maos2mex.h aolib.c
maos2mex.h:../scripts/maos2mex.py ../sys/scheduler.h ../maos/parms.h ../maos/types.h ../lib/accphi.h ../lib/cure.h \
	../lib/gensei.h ../lib/locfft.h ../lib/mkh.h ../lib/muv.h ../lib/proj.h ../lib/slaving.h \
	../lib/aos.h ../lib/draw.h ../lib/hyst.h ../lib/misc.h ../lib/mkw.h ../lib/pcg.h ../lib/psd.h ../lib/stfun.h\
	../lib/chol.h ../lib/fractal.h ../lib/kalman.h ../lib/mkdtf.h ../lib/mkz.h ../lib/petal.h ../lib/readcfg.h ../lib/turbulence.h\
	../lib/cn2est.h ../lib/genotf.h ../lib/laplacian.h ../lib/mkg.h ../lib/mtch.h ../lib/phase.h ../lib/servo.h ../lib/zernike.h ../lib/pywfs.h
	@$(SRC_DIR)/scripts/maos2mex.py $(SRC_DIR) maos2mex.h
aolib.c:../maos/parms.h ../maos/types.h ../scripts/lib2mex.py ../scripts/maos_parse.py aolib.h
	@$(SRC_DIR)/scripts/lib2mex.py $(SRC_DIR) aolib.c
nodist_aolib_SOURCES=aolib.c

CLEANFILES=aolib.c maos2mex.h
