/**
 * Standalone to postprocess MAOS PSF output and compute enclosed energe
 */
#include <unistd.h>
#include <getopt.h>

#include "../lib/aos.h"
static void calcenc_do(file_t *fpin, file_t *fpout, uint32_t magic, dmat *rvec){
    TIC;tic;
    dmat *psf=dreaddata(fpin, magic);
    toc("read");
    long ncomp=psf->nx>psf->ny?psf->nx:psf->ny;
    long ncomp2=ncomp*2;
    cmat *psf2=cnew(ncomp2, ncomp2);
    toc("cnew");
    cfft2plan(psf2, -1);
    toc("Plan");
    cembedd(psf2, psf, 0);
    toc("embed");
    cfft2(psf2, -1);
    toc("fft");
    dmat *psf3=NULL;
    creal2d(&psf3, 0, psf2, 1);
    toc("real2d");
    dfree(psf);
    cfree(psf2);
    double dk=1./(ncomp2);
    double dfx=dk;
    double dfy=dk;
    double pi2=2*M_PI;
    dscale(psf3, pow(ncomp2,-2));
    toc("scale");
    PDMAT(psf3, ppsf);
    int free_rvec=0;
    if(!rvec){
	free_rvec=1;
	rvec=dlinspace(0, 1, ncomp);
    }
    dmat *enc=dnew(rvec->nx, 1);
    for(long iy=0; iy<ncomp2; iy++){
	double k2y=pow((iy<ncomp?iy:iy-ncomp2), 2);
	for(long ix=0; ix<ncomp2; ix++){
	    double k2x=pow((ix<ncomp?ix:ix-ncomp2),2);
	    double k2=k2x+k2y;
	    for(long ir=0; ir<rvec->nx; ir++){
		double s=j0(sqrt(k2)*pi2*rvec->p[ir]);
		enc->p[ir]+=s*ppsf[iy][ix];
	    }
	}
    }
    toc("compute");
    dfree(psf3);
    dwritedata(fpout, enc);
    dfree(enc);
    if(free_rvec) dfree(rvec);
}
static void calcenc(const char *fn, dmat *rvec){
    char fnout[PATH_MAX];
    char *suffix=strstr(fn, ".bin");
    if(!suffix) error("%s has wrong suffix\n", suffix);
    snprintf(fnout, suffix-fn, fn);
    fnout[suffix-fn]='\0';
    strcat(fnout, "_enc.bin");
    
    file_t *fp=zfopen(fn, "r");
    file_t *fpout=zfopen(fn, "w");
    uint32_t magic=read_magic(fp, NULL);
    long nx, ny;
    if(iscell(magic)){
	write_magic(magic, fpout);
	zfreadlarr(fp, 2, &nx, &ny);
	zfwritelarr(fpout, 2, &nx, &ny);
	for(long i=0; i<nx*ny; i++){
	    info("%ld of %ld\n", i, nx*ny);
	    calcenc_do(fp, fpout, 0, rvec);
	}
    }else{
	calcenc_do(fp, fpout, magic, rvec);
    }
    zfclose(fp);
    zfclose(fpout);
}
static void usage(){
    error("Usage: calcenc [options] PSF1.bin PSF2.bin \n"
	  "-n num     :Maximum number of threads\n"
	  "-r num     :Maximum radius in unit of pixel to compute\n"
	  "-h         :Print this help"
	  );
}
int main(int argc, char *argv[]){
    if(argc<1){
	usage();
    }
    double rmax=0;
    double rstep=1;
    int nthread=0;
    int ipos=0;
    static struct option long_options[]={
	{"nthread", 1, 0, 'n'},
	{"radius_max", 1, 0, 'r'},
	{"radius_step", 1, 0, 'd'},
	{"help", 0, 0, 'h'},
	{NULL,0,0,0}
    };
    int option_index = 0;
    while(1){
	int c=getopt_long(argc, argv, "n:r:",
			  long_options, &option_index);
	if(c==-1) break;
	switch(c){
	case 'n':
	    nthread=strtol(optarg, NULL, 10);
	    break;
	case 'r':
	    rmax=strtod(optarg, NULL);
	    break;
	case 'd':
	    rstep=strtod(optarg, NULL);
	    break;
	case 'h':
	    usage();
	    break;
	default:
	    error("");
	}
    }
    ipos=optind;
    dmat *radius=NULL;
    if(rmax>0){
	radius=dlinspace(0, rstep, (long)ceil(rmax/rstep)+1);
    }
    info("ipos=%d, argv[ipos]=%s\n", ipos, argv[ipos]);
    for(; ipos<argc; ipos++){
	calcenc(argv[ipos], radius);
    }
    dfree(radius);
}
