SHELL := /bin/bash
all:read.mexa64 write.mexa64
mexopts=$(shell mex -h|grep -q largeArrayDims && echo -largeArrayDims)
read.mexa64:read.c io.c io.h
	mex $(mexopts) read.c io.c -lz
write.mexa64:write.c io.c io.h
	mex $(mexopts) write.c io.c -lz
clean:
	rm *.mexa64
