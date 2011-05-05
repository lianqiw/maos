#include "../lib/aos.h"
/**
   Warps of slaving for use standalone.

   Usage: slaving loc.bin HA.bin amp.bin

   
*/
int main(int argc, char **argv){
    if(argc<3){
	info2("Usage:\n"
	      "slaving loc.bin HA.bin [amp.bin] [thres]\n"
	      "loc.bin contains the grid. \n"
	      "HA.bin  contains the influence function from loc to destination\n"
	      "amp.bin contains the amplitude weighting in the destination (optional)\n"
	      "thres   contains the threshold to treat actuator as not coupled [0 1] (optional)\n"
	      "The result will be saved in slaving.bin\n"
	      );				
	exit(0);
    }
    loc_t *aloc=locread("%s",argv[1]);
    dsp *ha=spread("%s",argv[2]);
    double thres=argc>3?strtod(argv[3],NULL):0.5;
    spcell*has=spcellnew(1,1); has->p[0]=ha;
    spcell *slave=slaving(&aloc, has, NULL, NULL, thres, 1);
    spwrite(slave->p[0],"slaving");
}
