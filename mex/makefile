BUILD=${HOME}/work/aos/comp/optim
all:read.mexa64 write.mexa64 mkgmex.mexa64  proj.mexa64 cp genotfmex.mexa64
rand=crandcreate.mexa64 crandn.mexa64 crandp.mexa64 clrand.mexa64 mkw.mexa64
read.mexa64:read.c  io.c
	mex -largeArrayDims read.c io.c -lz
write.mexa64:write.c io.c
	mex -largeArrayDims write.c io.c -lz
mkg.o:../lib/mkg.c
	gcc -c -std=c99 ../lib/mkg.c -fPIC -DMATLAB_MEX_FILE=1
mkgmex.mexa64:mkgmex.c mkg.o
	mex  mkgmex.c mkg.o
cp:read.mexa64 write.mexa64 
	cp read.mexa64 write.mexa64 ~/work/matlab_routines
crandcreate.mexa64:crandcreate.c  
	mex crandcreate.c  -I../lib/ -lsupport
crandn.mexa64:crandn.c 
	mex crandn.c  -I../lib/ -lsupport 
crandp.mexa64:crandp.c 
	mex crandp.c  -I../lib/ -lsupport
clrand.mexa64:clrand.c 
	mex clrand.c  -I../lib/ -lsupport
mkw.mexa64:mkwmex.c 
	mex mkwmex.c interface.c -lsupport  
proj.mexa64:proj.c
	mex proj.c
MATHSRC=../lib/genotf.c ../lib/dmat.c ../lib/cmat.c ../lib/cmat_extra.c		\
../lib/loc.c ../lib/mathmisc.c ../lib/sys/misc.c ../lib/fft.c ../lib/dsp.c	\
../lib/csp.c ../lib/random.c ../lib/bin.c ../lib/sys/thread.c ../lib/path.c
genotfmex.mexa64:genotfmex.c interface.c $(MATHSRC)
	mex genotfmex.c interface.c $(MATHSRC)  -lz -lmwblas -lmwlapack -lfftw3 -lm 
clean:
	rm *.mexa64