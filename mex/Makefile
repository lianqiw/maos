all:read.mexa64 write.mexa64
read.mexa64:read.c io.c io.h
	mex -largeArrayDims read.c io.c -lz
write.mexa64:write.c io.c io.h
	mex -largeArrayDims write.c io.c -lz
clean:
	rm *.mexa64