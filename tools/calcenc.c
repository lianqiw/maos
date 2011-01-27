/**
 * Standalone to postprocess MAOS PSF output and compute enclosed energe
 */
#include <unistd.h>
#include <getopt.h>

#include "../lib/aos.h"
static void calcenc_do(file_t *fpin,   /**<Input file pointer*/
		       file_t *fpout,  /**<Output file pointer*/
		       uint32_t magic, /**<The magic number if already read*/
		       dmat *restrict dvec,     /**<Vector of the diameter*/
		       int type        /**<Type of encloded energy: 0: square, 1: circle*/
		       ){
    dmat *psf=dreaddata(fpin, magic);
    int free_dvec=0;
    if(!dvec){
	free_dvec=1;
	dvec=dlinspace(0, 1, psf->nx);
    }
    dmat *enc=denc(psf, dvec, type);
    dwritedata(fpout, enc);
    dfree(psf);
    dfree(enc);
    if(free_dvec) dfree(dvec);
}
static void calcenc(const char *fn, dmat *dvec, int type){
    char fnout[PATH_MAX];
    char *suffix=strstr(fn, ".bin");
    if(!suffix) error("%s has wrong suffix\n", suffix);
    snprintf(fnout, suffix-fn+1, fn);
    fnout[suffix-fn]='\0';
    switch(type){
    case -1:
	strcat(fnout, "_azmean.bin");
	break;
    case 0:
	strcat(fnout, "_ensquare.bin");
	break;
    case 1:
	strcat(fnout, "_encircle.bin");
	break;
    case 2:
	strcat(fnout, "_enstrip.bin");
	break;
    default:
	error("Not implemented\n");
    }
    file_t *fp=zfopen(fn, "r");
    file_t *fpout=zfopen(fnout, "w");
    info("%s --> %s\n", fn, fnout);
    uint32_t magic=read_magic(fp, NULL);
    long nx, ny;
    if(iscell(magic)){
	write_magic(magic, fpout);
	zfreadlarr(fp, 2, &nx, &ny);
	zfwritelarr(fpout, 2, &nx, &ny);
	for(long i=0; i<nx*ny; i++){
	    info("%ld of %ld\n", i, nx*ny);
	    calcenc_do(fp, fpout, 0, dvec, type);
	}
    }else{
	calcenc_do(fp, fpout, magic, dvec, type);
    }
    zfclose(fp);
    zfclose(fpout);
}
static void usage(){
    info2("Usage: calcenc [options] PSF1.bin PSF2.bin \n"
	  "Options are:\n"
	  "-n num --nthread: Maximum number of threads\n"
	  "-d num --diam   : Maximum diameter to compute enclosed energy. radius if type==-1\n"
	  "-s num --step   : Step in diamter\n"
	  "-t num --type   : type of computation:\n"
	  "\t\t\t0: ensquare energy inside a square\n"
	  "\t\t\t1: encircle energy inside a circle\n"
	  "\t\t\t2: enstriped energy inside a strip\n"
	  "\t\t\t-1: azimuthal average\n"
	  "-h         : Print this help"
	  );
}
int main(int argc, char *argv[]){
    if(argc==1){
	usage();
    }
    double dmax=0;
    double dstep=1;
    int nthread=0;
    int ipos=0;
    int type=0;//default is square
    static struct option long_options[]={
	{"nthread", 1, 0, 'n'},
	{"diam",    1, 0, 'd'},
	{"step",    1, 0, 's'},
	{"type",    1, 0, 't'},
	{"help",    0, 0, 'h'},
	{NULL,0,0,0}
    };
    int option_index = 0;
    while(1){
	int c=getopt_long(argc, argv, "n:d:s:t:h",
			  long_options, &option_index);
	if(c==-1) break;
	switch(c){
	case 'n':
	    nthread=strtol(optarg, NULL, 10);
	    break;
	case 'd':
	    dmax=strtod(optarg, NULL);
	    break;
	case 's':
	    dstep=strtod(optarg, NULL);
	    break;
	case 't':
	    type=strtol(optarg, NULL, 10);
	    break;
	case 'h':
	    usage();
	    break;
	default:
	    error("");
	}
    }
    ipos=optind;
    dmat *dvec=NULL;
    if(dmax>0){
	dvec=dlinspace(0, dstep, (long)ceil(dmax/dstep));
    }
    info("ipos=%d, argv[ipos]=%s\n", ipos, argv[ipos]);
    for(; ipos<argc; ipos++){
	calcenc(argv[ipos], dvec, type);
    }
    dfree(dvec);
}
