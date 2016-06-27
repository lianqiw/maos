AM_CFLAGS=$(MEXOPT) -I$(SRC_DIR)/lib -I$(SRC_DIR)/math 
AM_LDFLAGS=$(LDMEX)
#LDFLAGS goes before LDADD. Must put maos libraries before matlab libraries.
#The following should be in LDADD instead of AM_LDFLAGS to recompile exes when they are updated
LDADD=../lib/libaos.la $(LIBMEX) 
EXTRA_DIST=io.c io.h interface.h Makefile.mex compile.sh aolib.h
EXEEXT=.$(MEXEXT) #Automatic add suffix

.NOTPARALLEL:io.c

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
maos_LDADD+=-lcholmod $(LIBMEX) 


BUILT_SOURCES=maos2mex.h aolib.c aolibmex.c
maos2mex.h:../maos/parms.h ../maos/types.h ../lib/aos.h ../scripts/maos2mex.py
	$(SRC_DIR)/scripts/maos2mex.py $(SRC_DIR) $(SRC_DIR)/mex/maos2mex.h
aolib.c:../maos/parms.h ../maos/types.h ../scripts/aolib.py aolib.h
	$(SRC_DIR)/scripts/aolib.py $(SRC_DIR) $(SRC_DIR)/mex/aolib.c
aolib_SOURCES=aolibmex.c aolib.c
