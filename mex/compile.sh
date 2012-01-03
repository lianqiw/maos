mexopts="$(mex -h|grep -q largeArrayDims && echo -largeArrayDims) "
mex ${mexopts} read.c io.c -lz 
mex ${mexopts} write.c io.c -lz
